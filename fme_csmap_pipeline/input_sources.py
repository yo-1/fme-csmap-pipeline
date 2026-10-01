"""Input discovery and conversion shared by FME and QGIS. GPL-3.0-only."""
import csv
import glob
import json
import math
from pathlib import Path
import zipfile
import numpy as np
try:
    from .gsi_dem import parse_dem
    from .forest_dem import expand_sources, classify_sources, read_lem, autodetect_xyz
except ImportError:
    from gsi_dem import parse_dem
    from forest_dem import expand_sources, classify_sources, read_lem, autodetect_xyz

NODATA = -999999.0
DEFAULTS = dict(input_type='raster', recursive=True, normalize_rasters=False,
    input_crs='', vertical_reference='', text_delimiter='comma', text_encoding='utf-8-sig',
    text_skip_rows=1, text_columns=[1,2,3], text_axis_order='east_north',
    text_mode='grid', text_grid_cell=1.0, point_method='idw', point_radius=2.0,
    lidar_mode='classified', ground_classes=[2], pdal_path='', point_tile_size=512,
    forest_axis_order='north_east', forest_grid_cell='auto',
    forest_lem_coordinate_scale='auto', forest_nodata=-9999.0,
    max_points_auto=10000000, max_input_pixels=25000000, max_input_files=100000,
    smrf_cell=1.0, smrf_slope=0.15, smrf_threshold=0.5, smrf_scalar=1.25, smrf_window=18.0)
EXTENSIONS = {'raster': {'.tif','.tiff','.img','.asc','.vrt'}, 'gsi': {'.xml','.zip'},
              'text': {'.xyz','.csv','.txt'}, 'lidar': {'.las','.laz'},
              'forest': {'.lem','.csv','.txt','.xyz','.tif','.tiff','.zip'}}


def cancelled(feedback):
    if feedback is not None and feedback.isCanceled():
        raise RuntimeError('Input conversion cancelled; partial outputs are incomplete')


def say(feedback, text):
    cancelled(feedback)
    if feedback is not None: feedback.pushInfo(text)
    else: print(text, flush=True)


def validate_input(c):
    if c['input_type'] not in EXTENSIONS:
        raise ValueError('input_type: raster, gsi, text, lidar or forest')
    for key in ('recursive','normalize_rasters'):
        if type(c[key]) is not bool: raise ValueError(key+' must be boolean')
    for key in ('text_skip_rows','point_tile_size','max_points_auto','max_input_pixels','max_input_files'):
        if type(c[key]) is not int or c[key] < (0 if key=='text_skip_rows' else 1):
            raise ValueError(key+' must be an integer in range')
    if not 64 <= c['point_tile_size'] <= 2048: raise ValueError('point_tile_size: 64..2048')
    for key in ('text_grid_cell','point_radius','smrf_cell','smrf_slope','smrf_threshold','smrf_scalar','smrf_window'):
        if isinstance(c[key],bool) or not isinstance(c[key],(int,float)) or not math.isfinite(c[key]) or c[key]<=0:
            raise ValueError(key+' must be finite and positive')
    for key,choices in [('text_delimiter',('comma','tab','space','semicolon')),
        ('text_axis_order',('east_north','north_east')),('text_mode',('grid','points')),
        ('point_method',('idw','mean','min')),('lidar_mode',('classified','ground_only','auto')),
        ('forest_axis_order',('east_north','north_east'))]:
        if c[key] not in choices: raise ValueError('Invalid '+key)
    cols=c['text_columns']
    if not isinstance(cols,list) or len(cols)!=3 or len(set(cols))!=3 or any(type(i)is not int or i<1 for i in cols):
        raise ValueError('text_columns requires three distinct 1-based column numbers')
    classes=c['ground_classes']
    if not isinstance(classes,list) or not classes or any(type(i)is not int or not 0<=i<=255 for i in classes):
        raise ValueError('ground_classes must contain integers in 0..255')
    for key in ('input_crs','vertical_reference','text_encoding','pdal_path'):
        if not isinstance(c[key],str): raise ValueError(key+' must be text')
    if c['input_type']=='text' and not c['input_crs'].strip():
        raise ValueError('Text data requires input_crs (source CRS)')
    if c['forest_lem_coordinate_scale'] != 'auto':
        try: scale=float(c['forest_lem_coordinate_scale'])
        except (TypeError,ValueError) as exc: raise ValueError('forest_lem_coordinate_scale must be auto or positive') from exc
        if not math.isfinite(scale) or scale<=0: raise ValueError('forest_lem_coordinate_scale must be auto or positive')
    if c['forest_grid_cell'] != 'auto':
        try: grid_cell=float(c['forest_grid_cell'])
        except (TypeError,ValueError) as exc: raise ValueError('forest_grid_cell must be auto or positive') from exc
        if not math.isfinite(grid_cell) or grid_cell<=0: raise ValueError('forest_grid_cell must be auto or positive')
    if c['forest_nodata'] is not None and (isinstance(c['forest_nodata'],bool) or
            not isinstance(c['forest_nodata'],(int,float)) or not math.isfinite(c['forest_nodata'])):
        raise ValueError('forest_nodata must be null or a finite number')


def discover(inputs, base_dir, c):
    if not isinstance(inputs,list) or not inputs: raise ValueError('Specify input files or folders')
    files=[]; seen=set()
    for entry in inputs:
        p=Path(entry).expanduser()
        if not p.is_absolute():p=Path(base_dir)/p
        matches=[p] if p.exists() else [Path(x) for x in sorted(glob.glob(str(p),recursive=c['recursive']))]
        if not matches:raise ValueError('No input matched: '+str(entry))
        found=0
        for match in matches:
            items=sorted(match.rglob('*') if c['recursive'] else match.glob('*')) if match.is_dir() else [match]
            for f in items:
                if not f.is_file():continue
                if f.suffix.lower() not in EXTENSIONS[c['input_type']]:
                    if not match.is_dir():raise ValueError('Input extension does not match input_type: '+str(f))
                    continue
                found+=1
                path=str(f.resolve())
                if path not in seen:files.append(path);seen.add(path)
                if len(files)>c['max_input_files']:raise ValueError('Too many input files')
        if not found:raise ValueError('No supported files in: '+str(entry))
    return files


def crs_of(value, osr):
    crs=osr.SpatialReference()
    if not value or crs.SetFromUserInput(value)!=0:raise ValueError('Invalid/missing input CRS: '+value)
    if crs.IsCompound() or crs.GetAxesCount()!=2:
        raise ValueError('Use an explicit two-dimensional horizontal CRS; heights are preserved without vertical conversion')
    crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return crs


def raster_write(path, array, gt, crs, gdal):
    h,w=array.shape
    ds=gdal.GetDriverByName('GTiff').Create(str(path),w,h,1,gdal.GDT_Float32,
        options=['TILED=YES','COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=IF_SAFER'])
    if ds is None:raise RuntimeError('Cannot create '+str(path))
    ds.SetGeoTransform(gt);ds.SetProjection(crs.ExportToWkt())
    ds.GetRasterBand(1).SetNoDataValue(NODATA)
    for row in range(0,h,512):
        if ds.GetRasterBand(1).WriteArray(array[row:row+512],0,row)!=0:raise RuntimeError('Raster write failed')
    ds.FlushCache();ds=None


def has_world_file(path):
    """Return True when a TIFF has a same-stem world file (case-insensitive)."""
    path=Path(path)
    expected={path.stem.lower()+suffix for suffix in ('.tfw','.tifw','.wld')}
    return any(p.is_file() and p.name.lower() in expected for p in path.parent.iterdir())


# 森林航空レーザ成果のテキスト形式（forest_dem.read_lem）で欠測扱いにしている値。
# TIFFでも同じ値が標高として紛れ込むと、CS立体図に大きな乱れが出る。
FOREST_NODATA_CODES = (-9999.0, -1111.0)


def forest_nodata_codes(band_nodata, forest_nodata):
    """TIFFのNoData設定とは別に、欠測として扱うべき値を返す。

    gdal.WarpのsrcNodataは1バンドに1値しか指定できないため、NoData設定
    （band_nodata）と異なる欠測値が画素値に混在していると、その値が標高
    として扱われる。ここで返す値は mask_forest_nodata_codes() で置き換える。"""
    codes={float(v) for v in FOREST_NODATA_CODES}
    if forest_nodata is not None:codes.add(float(forest_nodata))
    if band_nodata is not None:codes.discard(float(band_nodata))
    return tuple(sorted(codes))


def mask_forest_nodata_codes(path, dest, codes, band_nodata, gdal, feedback=None):
    """codesに該当する画素をNoDataに置き換えたGeoTIFFをdestに作る。

    該当画素が無ければ何も作らず (path, band_nodata, 0) を返す（大半の入力で
    コピーを作らないため）。該当があればFloat32のコピーを作り、
    (dest, 新しいNoData値, 置換件数) を返す。読み書きは512行ずつ行い、
    入力全体をメモリに載せない。"""
    if not codes:return path,band_nodata,0
    src=gdal.Open(str(path))
    if src is None:raise ValueError('Cannot open forestry TIFF: '+str(path))
    band=src.GetRasterBand(1);w,h=src.RasterXSize,src.RasterYSize
    targets=np.asarray(codes,dtype='float64')
    count=0
    for row in range(0,h,512):
        cancelled(feedback)
        strip=band.ReadAsArray(0,row,w,min(512,h-row))
        count+=int(np.count_nonzero(np.isin(strip.astype('float64'),targets)))
    if count==0:
        src=None
        return path,band_nodata,0
    fill=NODATA if band_nodata is None else float(band_nodata)
    ds=gdal.GetDriverByName('GTiff').Create(str(dest),w,h,1,gdal.GDT_Float32,
        options=['TILED=YES','COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=IF_SAFER'])
    if ds is None:raise RuntimeError('Cannot create '+str(dest))
    ds.SetGeoTransform(src.GetGeoTransform())
    if src.GetProjection():ds.SetProjection(src.GetProjection())
    out=ds.GetRasterBand(1);out.SetNoDataValue(fill)
    for row in range(0,h,512):
        cancelled(feedback)
        strip=band.ReadAsArray(0,row,w,min(512,h-row)).astype('float32')
        strip[np.isin(strip.astype('float64'),targets)]=fill
        if out.WriteArray(strip,0,row)!=0:raise RuntimeError('Raster write failed')
    ds.FlushCache();ds=src=None
    return str(dest),fill,count



def normalize(path, dest, c, gdal, osr, feedback=None, override_nodata=True):
    """Create a disk-backed warped VRT; province-scale data stay virtual."""
    cancelled(feedback)
    src=gdal.Open(str(path))
    if src is None or src.RasterCount!=1:raise ValueError('Expected single-band elevation raster: '+str(path))
    band=src.GetRasterBand(1)
    if override_nodata and c['source_nodata'] is None and band.GetNoDataValue() is None:
        raise ValueError('Source raster has no NoData definition; specify source_nodata: '+str(path))
    if band.GetScale() not in (None,1) or band.GetOffset() not in (None,0):
        raise ValueError('Apply raster scale/offset to elevation values before import: '+str(path))
    source=crs_of(c['input_crs'] or src.GetProjection(),osr)
    kwargs=dict(format='VRT',srcSRS=source.ExportToWkt(),dstSRS=c['target_crs'],
        xRes=c['cell_size'],yRes=c['cell_size'],targetAlignedPixels=True,
        outputType=gdal.GDT_Float32,resampleAlg='bilinear',dstNodata=NODATA,errorThreshold=0)
    if override_nodata and c['source_nodata'] is not None:kwargs['srcNodata']=c['source_nodata']
    dest=Path(dest).with_suffix('.vrt')
    vrt=gdal.Warp(str(dest),src,**kwargs)
    if vrt is None:raise RuntimeError('Input reprojection failed: '+str(path))
    vrt.FlushCache();vrt=src=None
    cancelled(feedback)
    return str(dest)


def text_rows(path,c,feedback=None):
    separators={'comma':',','tab':'\t','semicolon':';'}
    with open(path,encoding=c['text_encoding'],newline='') as f:
        for _ in range(c['text_skip_rows']):next(f,None)
        rows=(line.split() for line in f) if c['text_delimiter']=='space' else csv.reader(f,delimiter=separators[c['text_delimiter']])
        for line,row in enumerate(rows,c['text_skip_rows']+1):
            if line%100000==0:cancelled(feedback)
            if not row or not any(x.strip() for x in row):continue
            try:a,b,z=(float(row[i-1]) for i in c['text_columns'])
            except (IndexError,ValueError) as exc:raise ValueError(f'{path}:{line}: check columns, delimiter and skipped rows') from exc
            if not all(math.isfinite(v) for v in (a,b,z)):
                raise ValueError(f'{path}:{line}: non-finite point')
            if c['source_nodata'] is not None and z==c['source_nodata']:continue
            yield (b,a,z) if c['text_axis_order']=='north_east' else (a,b,z)


def text_grid(path,dest,c,gdal,osr,feedback=None):
    crs=crs_of(c['input_crs'],osr)
    if not crs.IsProjected() or abs(crs.GetLinearUnits()-1)>1e-9:
        raise ValueError('Regular text grid requires a projected metre CRS; use point interpolation for longitude/latitude')
    xmin=ymin=math.inf;xmax=ymax=-math.inf;count=0
    for x,y,z in text_rows(path,c,feedback):
        xmin=min(xmin,x);xmax=max(xmax,x);ymin=min(ymin,y);ymax=max(ymax,y);count+=1
    if not count:raise ValueError('No valid text elevations: '+str(path))
    step=c['text_grid_cell'];w=round((xmax-xmin)/step)+1;h=round((ymax-ymin)/step)+1
    backing=Path(str(dest)+'.grid.bin');a=np.memmap(backing,dtype='float32',mode='w+',shape=(h,w))
    try:
        for row0 in range(0,h,512):a[row0:row0+512].fill(NODATA)
        for x,y,z in text_rows(path,c,feedback):
            col=(x-xmin)/step;row=(ymax-y)/step;ix=round(col);iy=round(row)
            if abs(col-ix)>1e-5 or abs(row-iy)>1e-5 or not 0<=ix<w or not 0<=iy<h:
                raise ValueError('Text points are not a regular grid at text_grid_cell; select point interpolation')
            if a[iy,ix]!=NODATA and not np.isclose(a[iy,ix],z,rtol=0,atol=1e-4):raise ValueError('Conflicting elevations at duplicate text coordinates')
            a[iy,ix]=z
        a.flush();raster_write(dest,a,(xmin-step/2,step,0,ymax+step/2,0,-step),crs,gdal)
    finally:del a;backing.unlink(missing_ok=True)
    return count


def infer_grid_cell(path,c,feedback=None):
    """Infer a square regular-grid interval from conventional row-major XYZ points."""
    previous=None;dx=dy=math.inf
    for x,y,z in text_rows(path,c,feedback):
        if previous is not None:
            px,py=previous
            change_x=abs(x-px);change_y=abs(y-py)
            if change_x>1e-9: dx=min(dx,change_x)
            if change_y>1e-9: dy=min(dy,change_y)
        previous=(x,y)
    candidates=[value for value in (dx,dy) if math.isfinite(value)]
    if not candidates: raise ValueError('Cannot infer forest CSV/XYZ grid interval from fewer than two coordinate positions')
    step=min(candidates)
    if len(candidates)==2 and not math.isclose(dx,dy,rel_tol=1e-5,abs_tol=1e-8):
        raise ValueError('Forest CSV/XYZ has different X/Y intervals; specify forest_grid_cell after checking the grid')
    return step



def normalize_gsi_groups(paths,work,c,gdal,osr,feedback=None):
    groups=[];base=None
    for path in paths:
        ds=gdal.Open(path);srs=ds.GetSpatialRef();gt=ds.GetGeoTransform()
        same=False
        if base is not None:
            bs,bg=base
            same=bool(srs.IsSame(bs) and np.allclose([gt[1],gt[5]],[bg[1],bg[5]],rtol=1e-7,atol=0))
            if same:
                shifts=[(gt[0]-bg[0])/gt[1],(gt[3]-bg[3])/gt[5]]
                same=all(abs(v-round(v))<1e-5 for v in shifts)
        if not same:
            groups.append([]);base=(srs.Clone(),gt)
        groups[-1].append(path);ds=None
    results=[]
    for i,files in enumerate(groups):
        cancelled(feedback)
        mosaicpath=work/f'gsi_group_{i:06d}.vrt'
        mosaic=gdal.BuildVRT(str(mosaicpath),files,resolution='highest',srcNodata=NODATA,VRTNodata=NODATA,strict=True)
        if mosaic is None:raise RuntimeError('GSI mosaic failed')
        mosaic.FlushCache();mosaic=None
        results.append(normalize(mosaicpath,work/f'aligned_group_{i:06d}.tif',
            {**c,'input_crs':'','max_input_pixels':c['max_pixels']},gdal,osr,feedback,False))
    return results


def prepare_inputs(c,work,gdal,osr,feedback=None):
    if c['input_type']=='raster' and not c['normalize_rasters'] and not c['input_crs']:
        return c, {'mode':'raster_passthrough'}
    work=Path(work);work.mkdir(parents=True,exist_ok=False)
    report={'status':'running','input_type':c['input_type'],'sources':[],
            'vertical_reference':c['vertical_reference'],'vertical_conversion':False}
    reportpath=work/'input_report.json'
    def save():reportpath.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    save();outputs=[];ordinal=0
    try:
        if c['input_type']=='lidar' or (c['input_type']=='text' and c['text_mode']=='points'):
            try:from .point_inputs import points_to_rasters
            except ImportError:from point_inputs import points_to_rasters
            outputs,detail=points_to_rasters(c,work,gdal,osr,feedback)
            report.update(detail)
        elif c['input_type']=='forest':
            expanded,archives=expand_sources(c['inputs'],work)
            records=classify_sources(expanded)
            report['archives']=archives;report['detected']=[];report['ignored']=[]
            for record in records:
                cancelled(feedback);ordinal+=1;path=record['path'];kind=record['kind']
                say(feedback,f'Forest DEM ({kind}): {path}')
                if kind=='lem':
                    if not c['input_crs'].strip():
                        raise ValueError('LEM requires input_crs because the zone number does not identify the datum realization')
                    array,gt,detail=read_lem(path,record['metadata'],c['max_input_pixels'],
                        c['forest_lem_coordinate_scale'],feedback)
                    raw=work/f'forest_lem_{ordinal:06d}.tif'
                    raster_write(raw,array,gt,crs_of(c['input_crs'],osr),gdal);del array
                    source={**c,'source_nodata':None}
                    report['detected'].append({'path':str(path),'kind':'lem','metadata':str(record['metadata']),
                                               'zone':detail.get('zone'),'coordinate_scale':detail['coordinate_scale']})
                elif kind=='grid':
                    try: detected=autodetect_xyz(path,c['forest_axis_order'])
                    except ValueError as exc:
                        report['ignored'].append({'path':str(path),'reason':str(exc)});ordinal-=1;continue
                    if not c['input_crs'].strip(): raise ValueError('Forest CSV/XYZ grid requires input_crs')
                    grid_cell=(infer_grid_cell(path,{**c,**detected,'source_nodata':c['forest_nodata']},feedback)
                               if c['forest_grid_cell']=='auto' else float(c['forest_grid_cell']))
                    source={**c,**detected,'source_nodata':c['forest_nodata'],
                            'text_grid_cell':grid_cell}
                    raw=work/f'forest_grid_{ordinal:06d}.tif'
                    count=text_grid(path,raw,source,gdal,osr,feedback)
                    report['detected'].append({'path':str(path),'kind':'xyz_grid','points':count,
                                               'grid_cell':grid_cell,**detected})
                else:
                    ds=gdal.Open(str(path))
                    if ds is None: raise ValueError('Cannot open forestry TIFF: '+str(path))
                    embedded=bool(ds.GetProjection());band_nodata=ds.GetRasterBand(1).GetNoDataValue();ds=None
                    if not embedded and not has_world_file(path):
                        raise ValueError('Non-GeoTIFF forestry TIFF requires a same-stem TFW/TIFW/WLD: '+str(path))
                    if not embedded and not c['input_crs'].strip():
                        raise ValueError('TIFF/world-file input has no embedded CRS; specify input_crs: '+str(path))
                    codes=forest_nodata_codes(band_nodata,c['forest_nodata'])
                    raw,band_nodata,masked=mask_forest_nodata_codes(
                        path,work/f'forest_masked_{ordinal:06d}.tif',codes,band_nodata,gdal,feedback)
                    if masked:
                        say(feedback,f'NoData設定と異なる欠測値（{", ".join(f"{v:g}" for v in codes)}）を'
                                     f'{masked}画素検出し、NoDataとして扱いました: {path}')
                    source={**c,'source_nodata':band_nodata if band_nodata is not None else c['forest_nodata']}
                    report['detected'].append({'path':str(path),'kind':'geotiff' if embedded else 'tiff_worldfile',
                                               'extra_nodata_cells':masked})
                aligned=work/f'aligned_{ordinal:06d}.tif'
                outputs.append(normalize(raw,aligned,source,gdal,osr,feedback,kind=='raster'))
                report['sources'].append({'path':str(path),'kind':kind});save()
        else:
            gsi_datums=set()
            for path in c['inputs']:
                say(feedback,'Input: '+path)
                if c['input_type']=='gsi':
                    members=[]
                    if Path(path).suffix.lower()=='.zip':
                        with zipfile.ZipFile(path) as archive:
                            for entry in sorted(archive.infolist(),key=lambda e:e.filename):
                                if entry.is_dir() or not entry.filename.lower().endswith('.xml'):continue
                                if entry.file_size>256*1024*1024:raise ValueError('ZIP XML member exceeds 256 MiB')
                                # Members are read directly; never extract user-controlled paths.
                                members.append(entry.filename)
                                with archive.open(entry) as f:
                                    data=f.read(256*1024*1024+1)
                                for a,gt,srs,meta in parse_dem(data,c['max_input_pixels'],c['input_crs'],feedback):
                                    ordinal+=1
                                    raw=work/f'gsi_{ordinal:06d}.tif'
                                    crs=crs_of(srs,osr)
                                    if not crs.IsGeographic():raise ValueError('GSI input_crs must be geographic latitude/longitude')
                                    raster_write(raw,a,gt,crs,gdal)
                                    del a
                                    outputs.append(str(raw))
                                    report['sources'].append({'path':path,'member':entry.filename,**meta})
                                    gsi_datums.add(meta['original_srs'].lower())
                    else:
                        if Path(path).stat().st_size>256*1024*1024:raise ValueError('XML exceeds 256 MiB')
                        for a,gt,srs,meta in parse_dem(Path(path).read_bytes(),c['max_input_pixels'],c['input_crs'],feedback):
                            ordinal+=1
                            raw=work/f'gsi_{ordinal:06d}.tif';crs=crs_of(srs,osr)
                            if not crs.IsGeographic():raise ValueError('GSI input_crs must be geographic')
                            raster_write(raw,a,gt,crs,gdal);del a
                            outputs.append(str(raw))
                            report['sources'].append({'path':path,**meta});gsi_datums.add(meta['original_srs'].lower())
                    if len(outputs)>c['max_input_files']:raise ValueError('Too many DEM members')
                else:
                    ordinal+=1;raw=path
                    if c['input_type']=='text':
                        raw=work/f'text_{ordinal:06d}.tif'
                        count=text_grid(path,raw,c,gdal,osr,feedback)
                        report['sources'].append({'path':path,'points':count})
                    else:report['sources'].append({'path':path})
                    outputs.append(normalize(raw,work/f'aligned_{ordinal:06d}.tif',c,gdal,osr,feedback,c['input_type']=='raster'))
                save()
            # A horizontal CRS override cannot make old/new vertical realizations agree.
            if any('jgd2024' in d for d in gsi_datums) and any('jgd2011' in d or 'jgd2000' in d for d in gsi_datums):
                raise ValueError('GSI old/new datum labels are mixed. Unify elevation realizations before merging.')
            if c['input_type']=='gsi' and outputs:
                outputs=normalize_gsi_groups(outputs,work,c,gdal,osr,feedback)
        if not outputs:raise ValueError('No usable DEMs found (XML may contain metadata only)')
        report.update(status='completed',prepared_rasters=outputs);save()
        return {**c,'inputs':outputs,'source_nodata':None},report
    except BaseException as exc:
        report.update(status='cancelled' if feedback is not None and feedback.isCanceled() else 'failed',error=str(exc));save();raise
