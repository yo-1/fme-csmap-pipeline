"""LAS/LAZ and irregular text -> buffered DEM blocks via external PDAL."""
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
try:
    from .input_sources import cancelled, say, crs_of, text_rows, NODATA
except ImportError:
    from input_sources import cancelled, say, crs_of, text_rows, NODATA


def executable(c):
    requested=c['pdal_path'].strip()
    path=str(Path(requested).expanduser().resolve()) if requested else shutil.which('pdal')
    if not path or not Path(path).is_file():
        raise ValueError('PDAL executable not found. Install the supplied conda environment and select its pdal.exe (not python.exe or pdal_wrench.exe).')
    if Path(path).suffix.lower() in ('.bat','.cmd'):
        raise ValueError('Select pdal.exe itself, not a shell script')
    return path


def run_process(exe,args,log,feedback=None,json_output=False):
    cancelled(feedback)
    # File-backed stdout/stderr avoids pipe deadlock and unbounded in-memory logs.
    # A conda pdal.exe may be invoked from QGIS without activating conda.
    env=os.environ.copy();parent=Path(exe).parent
    env['PATH']=str(parent)+os.pathsep+env.get('PATH','')
    sharedirs=[parent.parent/'share',parent.parent/'Library'/'share']
    for var,sub in [('PROJ_DATA','proj'),('GDAL_DATA','gdal')]:
        candidate=next((d/sub for d in sharedirs if (d/sub).is_dir()),None)
        if candidate:env[var]=str(candidate)
    if 'PROJ_DATA' in env:env['PROJ_LIB']=env['PROJ_DATA']
    env.pop('GDAL_DRIVER_PATH',None)
    with open(log,'wb') as stdout,open(str(log)+'.stderr','wb') as stderr:
        process=subprocess.Popen([exe,*args],stdout=stdout,stderr=stderr,stdin=subprocess.DEVNULL,
            shell=False,env=env,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        try:
            while process.poll() is None:
                cancelled(feedback)
                time.sleep(.15)
        except BaseException:
            process.terminate()
            try:process.wait(timeout=5)
            except subprocess.TimeoutExpired:process.kill();process.wait()
            raise
    if process.returncode:
        with open(str(log)+'.stderr','rb') as f:
            f.seek(max(0,Path(str(log)+'.stderr').stat().st_size-4000))
            error=f.read().decode('utf-8',errors='replace')
        raise RuntimeError(f'PDAL failed ({process.returncode}); {log}: {error}')
    if json_output:return json.loads(Path(log).read_text(encoding='utf-8-sig'))
    return None


def overlap(a,b):
    return a[0]<=b[2] and a[2]>=b[0] and a[1]<=b[3] and a[3]>=b[1]


def transform_bounds(bounds,source,target,osr):
    b=osr.CoordinateTransformation(source,target).TransformBounds(*bounds,41)
    if len(b)!=4 or not all(math.isfinite(v) for v in b):raise ValueError('Point bounds reprojection failed')
    return tuple(b)


def stage_unicode_source(path,work,number,feedback=None):
    """Return an ASCII-only PDAL path for a Unicode-named LAS/LAZ source."""
    source=Path(path).resolve()
    if str(source).isascii():return str(source),None,None
    stage=Path(work)/'pdal_ascii_inputs'
    if not str(stage.resolve()).isascii():
        raise ValueError('PDAL cannot read the Unicode input path and the output work path also contains non-ASCII characters. Choose an ASCII-only output folder such as C:\\gis_data\\output.')
    stage.mkdir(parents=True,exist_ok=True)
    destination=stage/f'input_{number:06d}{source.suffix.lower()}'
    if destination.exists():destination.unlink()
    try:
        os.link(source,destination);method='hardlink'
    except OSError:
        say(feedback,'Unicode path detected; copying point input to the ASCII PDAL work area: '+str(source))
        shutil.copy2(source,destination);method='copy'
    return str(destination),destination,method


def inspect_sources(c,work,exe,osr,feedback):
    target=crs_of(c['target_crs'],osr);records=[];staged=[]
    for number,path in enumerate(c['inputs']):
        say(feedback,'Selected point input: '+path)
        if c['input_type']=='text':
            cleaned=work/f'points_{number:06d}.csv';count=0
            xmin=ymin=math.inf;xmax=ymax=-math.inf
            with cleaned.open('w',encoding='ascii',newline='') as f:
                f.write('X,Y,Z\n')
                for x,y,z in text_rows(path,c,feedback):
                    f.write(f'{x:.15g},{y:.15g},{z:.15g}\n');count+=1
                    xmin=min(xmin,x);xmax=max(xmax,x);ymin=min(ymin,y);ymax=max(ymax,y)
            if not count:raise ValueError('No valid text points: '+path)
            source=crs_of(c['input_crs'],osr)
            bounds=(xmin,ymin,xmax,ymax)
            record=dict(path=path,read_path=str(cleaned),reader='readers.text',count=count)
        else:
            read_path,staged_path,stage_method=stage_unicode_source(path,work,number,feedback)
            if staged_path is not None:
                staged.append(staged_path)
                say(feedback,'PDAL ASCII input: '+read_path)
            data=run_process(exe,['info','--summary',read_path],work/f'info_{number:06d}.json',feedback,True)
            summary=data.get('summary',{})
            b=summary.get('bounds',{})
            if not all(k in b for k in ('minx','miny','maxx','maxy')):
                raise ValueError('PDAL summary has no readable bounds: '+path)
            bounds=tuple(b[k] for k in ('minx','miny','maxx','maxy'))
            sr=summary.get('srs',{})
            value=sr if isinstance(sr,str) else (sr.get('compoundwkt') or sr.get('wkt') or sr.get('horizontal') or '')
            source=crs_of(c['input_crs'] or value,osr)
            record=dict(path=path,read_path=read_path,reader='readers.las',count=int(summary.get('num_points',0)),source_srs=value)
            if stage_method:record['pdal_path_staging']=stage_method
            if record['count']<=0:raise ValueError('Empty LAS/LAZ or missing point count: '+path)
        if not all(math.isfinite(v) for v in bounds):raise ValueError('Non-finite point bounds')
        record.update(crs=source.ExportToWkt(),bounds=list(transform_bounds(bounds,source,target,osr)))
        records.append(record)
    return records,staged


def build_pipeline(records,bounds,c,destination,classified_path=None):
    """Explicit tags keep independent reader branches from inheriting each other."""
    cell=c['cell_size'];auto=c['input_type']=='lidar' and c['lidar_mode']=='auto'
    halo=math.ceil((c['point_radius']+(2*c['smrf_window'] if auto else 0))/cell)+2
    west,south,east,north=bounds
    extended=(west-halo*cell,south-halo*cell,east+halo*cell,north+halo*cell)
    crop=f'([{extended[0]},{extended[2]}],[{extended[1]},{extended[3]}])'
    pipeline=[];tags=[]
    for i,r in enumerate(records):
        tag=f'R{i}';proj=f'P{i}';cut=f'C{i}';keep=f'K{i}'
        reader=dict(type=r['reader'],filename=r['read_path'],override_srs=r['crs'],tag=tag)
        pipeline.append(reader)
        pipeline.append(dict(type='filters.reprojection',in_srs=r['crs'],out_srs=c['target_crs'],inputs=[tag],tag=proj))
        pipeline.append(dict(type='filters.crop',bounds=crop,inputs=[proj],tag=cut))
        if c['input_type']=='lidar':
            expression='Withheld == 0 && Classification != 7 && Classification != 18'
            if c['lidar_mode']=='classified':
                expression+=' && ('+' || '.join(f'Classification == {k}' for k in c['ground_classes'])+')'
            pipeline.append(dict(type='filters.expression',expression=expression,inputs=[cut],tag=keep))
            tags.append(keep)
        else:tags.append(cut)
    if auto:
        pipeline.append(dict(type='filters.merge',inputs=tags,tag='MERGED'))
        # Reset pre-existing labels: auto mode must not retain old class-2 points.
        pipeline.append(dict(type='filters.assign',value='Classification = 0',inputs=['MERGED'],tag='RESET'))
        pipeline.append(dict(type='filters.smrf',cell=c['smrf_cell'],slope=c['smrf_slope'],
            threshold=c['smrf_threshold'],scalar=c['smrf_scalar'],window=c['smrf_window'],
            returns='first,last,intermediate,only',inputs=['RESET'],tag='CLASSIFIED'))
        if classified_path:
            pipeline.append(dict(type='writers.las',filename=str(classified_path),inputs=['CLASSIFIED'],
                minor_version=4,dataformat_id=6,scale_x=.001,scale_y=.001,scale_z=.001,
                offset_x='auto',offset_y='auto',offset_z='auto'))
        pipeline.append(dict(type='filters.expression',expression='Classification == 2',inputs=['CLASSIFIED'],tag='GROUND'))
        tags=['GROUND']
    pipeline.append(dict(type='writers.gdal',filename=str(destination),inputs=tags,
        resolution=cell,origin_x=extended[0],origin_y=extended[1],
        width=round((east-west)/cell)+2*halo,height=round((north-south)/cell)+2*halo,
        radius=c['point_radius'],output_type=c['point_method'],window_size=0,nodata=NODATA,
        data_type='float',gdaldriver='GTiff',gdalopts='TILED=YES,COMPRESS=DEFLATE,BIGTIFF=IF_SAFER',
        override_srs=c['target_crs'],allow_empty=True))
    return {'pipeline':pipeline},halo


def points_to_rasters(c,work,gdal,osr,feedback=None):
    exe=executable(c)
    run_process(exe,['--version'],work/'pdal_version.txt',feedback)
    version=(work/'pdal_version.txt').read_text(errors='replace')
    match=re.search(r'\b(\d+)\.(\d+)\.(\d+)\b',version)
    if not match or tuple(map(int,match.groups()))<(2,9,0):
        raise ValueError('PDAL >=2.9 is required for point processing')
    # c['inputs'] is the validated config/CLI selection. Do not consult any
    # remembered or test path during point processing.
    records,staged=inspect_sources(c,work,exe,osr,feedback)
    (work/'point_sources.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
    cell=c['cell_size'];size=c['point_tile_size'];span=cell*size
    indices=set()
    for record in records:
        a,b,d,e=record['bounds']
        ix0,ix1=math.floor(a/span),math.floor(d/span)
        iy0,iy1=math.floor(b/span),math.floor(e/span)
        if (ix1-ix0+1)*(iy1-iy0+1)>c['max_input_files']:raise ValueError('Too many point DEM blocks')
        indices.update((ix,iy) for ix in range(ix0,ix1+1) for iy in range(iy0,iy1+1))
        if len(indices)>c['max_input_files']:raise ValueError('Too many point DEM blocks')
    xs=[i for i,j in indices];ys=[j for i,j in indices]
    mosaic_pixels=(max(xs)-min(xs)+1)*(max(ys)-min(ys)+1)*size*size
    if mosaic_pixels>c['max_pixels']:raise ValueError('Point output mosaic exceeds max_pixels; reduce region or use coarser cells')
    auto=c['input_type']=='lidar' and c['lidar_mode']=='auto'
    buffer=(math.ceil((c['point_radius']+(2*c['smrf_window'] if auto else 0))/cell)+2)*cell
    if (size+2*math.ceil(buffer/cell))**2>c['max_input_pixels']:
        raise ValueError('Buffered point block exceeds max_input_pixels')
    outputs=[];empty=0
    try:
      for i,(ix,iy) in enumerate(sorted(indices)):
        cancelled(feedback)
        bounds=(ix*span,iy*span,(ix+1)*span,(iy+1)*span)
        expanded=(bounds[0]-buffer,bounds[1]-buffer,bounds[2]+buffer,bounds[3]+buffer)
        selected=[r for r in records if overlap(r['bounds'],expanded)]
        name=f'block_{i:06d}';bufferfile=work/f'{name}_buffer.tif'
        classified=work/f'{name}_classified.laz' if auto else None
        pipeline,halo=build_pipeline(selected,bounds,c,bufferfile,classified)
        if auto:
            # Stream the buffered subset to disk before running non-streaming SMRF.
            # This bounds memory by local point count, not original file sizes.
            cutfile=work/f'{name}_subset.laz'
            merge_index=next(j for j,stage in enumerate(pipeline['pipeline']) if stage['type']=='filters.merge')
            branches=pipeline['pipeline'][:merge_index]
            inputs=pipeline['pipeline'][merge_index]['inputs']
            extract={'pipeline':branches+[dict(type='writers.las',filename=str(cutfile),inputs=inputs,
                minor_version=4,dataformat_id=6,scale_x=.001,scale_y=.001,scale_z=.001,
                offset_x='auto',offset_y='auto',offset_z='auto')]}
            extractpath=work/f'{name}_subset_pipeline.json'
            extractpath.write_text(json.dumps(extract,ensure_ascii=False,indent=2),encoding='utf-8')
            say(feedback,f'Extracting buffered point subset {i+1}/{len(indices)}')
            run_process(exe,['pipeline',str(extractpath),'--stream'],work/f'{name}_subset.log',feedback)
            summary=run_process(exe,['info','--summary',str(cutfile)],work/f'{name}_subset_info.json',feedback,True).get('summary',{})
            count=int(summary.get('num_points',0))
            if count==0:
                empty+=1;cutfile.unlink();continue
            if count>c['max_points_auto']:
                raise ValueError('Buffered subset exceeds max_points_auto. Reduce point_tile_size or increase the limit according to available memory. Subset LAZ has been retained.')
            rest=pipeline['pipeline'][merge_index+1:]
            pipeline={'pipeline':[dict(type='readers.las',filename=str(cutfile),override_srs=c['target_crs'],tag='MERGED')]+rest}
        pipelinefile=work/f'{name}_pipeline.json'
        pipelinefile.write_text(json.dumps(pipeline,ensure_ascii=False,indent=2),encoding='utf-8')
        say(feedback,f'Point DEM block {i+1}/{len(indices)}; sources={len(selected)}')
        args=['pipeline',str(pipelinefile)]
        if not auto:args.append('--stream')
        run_process(exe,args,work/f'{name}_pdal.log',feedback)
        ds=gdal.Open(str(bufferfile))
        if ds is None or ds.RasterCount!=1:raise RuntimeError('PDAL did not create a single-band DEM')
        gt=ds.GetGeoTransform();expected=(bounds[0]-halo*cell,cell,0,bounds[3]+halo*cell,0,-cell)
        if any(abs(a-b)>1e-6 for a,b in zip(gt,expected)) or ds.RasterXSize!=size+2*halo or ds.RasterYSize!=size+2*halo:
            raise RuntimeError('PDAL raster grid differs from requested grid; refusing to shift terrain')
        a=ds.GetRasterBand(1).ReadAsArray(halo,halo,size,size)
        import numpy as np
        if a is None:raise RuntimeError('PDAL DEM read failed')
        if np.any(np.isfinite(a)&(a!=NODATA)):
            dest=work/f'{name}.tif'
            output=gdal.Translate(str(dest),ds,srcWin=[halo,halo,size,size],format='GTiff',
                creationOptions=['TILED=YES','COMPRESS=DEFLATE','BIGTIFF=IF_SAFER'])
            if output is None:raise RuntimeError('Point DEM crop failed')
            output.FlushCache();output=None;outputs.append(str(dest))
        else:empty+=1
        ds=None
        # Reproducible transient raster; classified LAZ, pipelines and logs remain.
        bufferfile.unlink()
        if auto:cutfile.unlink()
    finally:
        for path in staged:
            try:path.unlink()
            except FileNotFoundError:pass
        stage_dir=work/'pdal_ascii_inputs'
        if stage_dir.is_dir():
            try:stage_dir.rmdir()
            except OSError:pass
    return outputs,dict(sources=records,pdal_version=version.strip(),blocks=len(indices),empty_blocks=empty,
        classification_mode=c['lidar_mode'] if c['input_type']=='lidar' else 'text_ground_points',
        interpolation=c['point_method'],radius_m=c['point_radius'],fallback_fill=False,
        note='Auto-classified LAZ blocks contain overlapping buffers; inspect before use. Original point files are unchanged.')
