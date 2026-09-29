"""Serial GDAL XYZ PNG writer; no external process or gdal2tiles dependency."""
import json
import math
from pathlib import Path

HALF_WORLD = math.pi * 6378137.0
DEFAULTS = dict(xyz_enabled=True, xyz_min_zoom=12, xyz_max_zoom=18,
                xyz_max_tiles=100000)


def validate_xyz(c):
    if not isinstance(c['xyz_enabled'], bool):
        raise ValueError('xyz_enabled must be boolean')
    for key in ('xyz_min_zoom', 'xyz_max_zoom', 'xyz_max_tiles'):
        if type(c[key]) is not int:
            raise ValueError(key+' must be an integer')
    if not 0 <= c['xyz_min_zoom'] <= c['xyz_max_zoom'] <= 24:
        raise ValueError('XYZ zoom requires 0 <= min <= max <= 24')
    if c['xyz_max_tiles'] <= 0:
        raise ValueError('xyz_max_tiles must be positive')


def tile_bounds(z, x, y):
    span = 2*HALF_WORLD/(1 << z)
    return (-HALF_WORLD+x*span, HALF_WORLD-(y+1)*span,
            -HALF_WORLD+(x+1)*span, HALF_WORLD-y*span)


def tile_range(bounds, z):
    west, south, east, north = bounds
    if not all(math.isfinite(v) for v in bounds) or west >= east or south >= north:
        raise ValueError('Invalid XYZ extent')
    n = 1 << z
    span = 2*HALF_WORLD/n
    # Exclusive east/south edges avoid duplicate tiles on exact boundaries.
    x0 = max(0, math.floor((west+HALF_WORLD)/span))
    x1 = min(n-1, math.ceil((east+HALF_WORLD)/span)-1)
    y0 = max(0, math.floor((HALF_WORLD-north)/span))
    y1 = min(n-1, math.ceil((HALF_WORLD-south)/span)-1)
    return x0, x1, y0, y1


def write_xyz(source, destination, c, gdal, osr, feedback=None):
    validate_xyz({**DEFAULTS, **c})
    c = {**DEFAULTS, **c}
    def cancel():
        if feedback is not None and feedback.isCanceled():
            raise RuntimeError('XYZ cancelled; output is incomplete')
    cancel()
    src = gdal.Open(str(source))
    if src is None or src.RasterCount != 4:
        raise ValueError('XYZ input must be a georeferenced RGBA raster')
    if any(src.GetRasterBand(i).DataType != gdal.GDT_Byte for i in range(1, 5)):
        raise ValueError('XYZ input must be Byte RGBA')
    srs = osr.SpatialReference()
    if not src.GetProjection() or srs.ImportFromWkt(src.GetProjection()) != 0:
        raise ValueError('XYZ input CRS is missing or invalid')
    merc = osr.SpatialReference()
    merc.ImportFromEPSG(3857)
    for crs in (srs, merc):
        crs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    gt = src.GetGeoTransform()
    if gt[2] != 0 or gt[4] != 0 or gt[1] <= 0 or gt[5] >= 0:
        raise ValueError('XYZ input must be north-up')
    bounds = osr.CoordinateTransformation(srs, merc).TransformBounds(
        gt[0], gt[3]+gt[5]*src.RasterYSize,
        gt[0]+gt[1]*src.RasterXSize, gt[3], 41)
    ranges = [(z, tile_range(bounds, z)) for z in range(c['xyz_min_zoom'], c['xyz_max_zoom']+1)]
    total = sum(max(0,b-a+1)*max(0,d-e+1) for _,(a,b,e,d) in ranges)
    if total == 0 or total > c['xyz_max_tiles']:
        raise ValueError(f'XYZ candidate tiles={total:,}; limit={c["xyz_max_tiles"]:,}. Reduce extent/max zoom or increase xyz_max_tiles.')
    out = Path(destination)
    out.mkdir(parents=True, exist_ok=False)
    info = dict(status='running', scheme='xyz', crs='EPSG:3857', tile_size=256,
                minzoom=c['xyz_min_zoom'], maxzoom=c['xyz_max_zoom'],
                bounds_3857=list(bounds), template='{z}/{x}/{y}.png',
                candidate_tiles=total, written_tiles=0, skipped_empty_tiles=0)
    def save():
        (out/'xyz.json').write_text(json.dumps(info, indent=2), encoding='utf-8')
    save()
    try:
        done = 0
        for z,(x0,x1,y0,y1) in ranges:
            for x in range(x0,x1+1):
                for y in range(y0,y1+1):
                    cancel()
                    tile = gdal.Warp('', src, format='MEM', dstSRS='EPSG:3857',
                        outputBounds=tile_bounds(z,x,y), width=256, height=256,
                        outputType=gdal.GDT_Byte, srcAlpha=True, dstAlpha=True,
                        resampleAlg='average', errorThreshold=0, warpMemoryLimit=64,
                        callback=lambda fraction,message,data: 0 if feedback is not None and feedback.isCanceled() else 1)
                    cancel()
                    if tile is None or tile.RasterCount != 4:
                        raise RuntimeError('XYZ warp failed')
                    pixels = tile.ReadAsArray()
                    if pixels is None:
                        raise RuntimeError('XYZ pixel read failed')
                    if pixels[3].any():
                        folder = out/str(z)/str(x)
                        folder.mkdir(parents=True,exist_ok=True)
                        # Plain image dataset avoids sidecar georeferencing files.
                        plain = gdal.GetDriverByName('MEM').Create('',256,256,4,gdal.GDT_Byte)
                        for i, ci in enumerate((gdal.GCI_RedBand,gdal.GCI_GreenBand,gdal.GCI_BlueBand,gdal.GCI_AlphaBand)):
                            band = plain.GetRasterBand(i+1)
                            band.WriteArray(pixels[i])
                            band.SetColorInterpretation(ci)
                        png = gdal.GetDriverByName('PNG').CreateCopy(str(folder/f'{y}.png'),plain)
                        if png is None:
                            raise RuntimeError('PNG write failed')
                        png.FlushCache()
                        png = plain = None
                        info['written_tiles'] += 1
                    else:
                        info['skipped_empty_tiles'] += 1
                    tile = None
                    done += 1
                    if feedback is not None:
                        feedback.setProgress(85+14*done/total)
            save()
            message = f'XYZ z={z}: {done}/{total}, PNG={info["written_tiles"]}'
            if feedback is not None: feedback.pushInfo(message)
            else: print(message,flush=True)
        if info['written_tiles'] == 0:
            raise ValueError('No visible XYZ tiles; choose a higher zoom or check alpha')
        info['status'] = 'completed'
        save()
        return info
    except BaseException as exc:
        info.update(status='cancelled' if feedback is not None and feedback.isCanceled() else 'failed', error=str(exc))
        save()
        raise
    finally:
        src = None


if __name__ == '__main__':
    import argparse
    from osgeo import gdal, osr
    p = argparse.ArgumentParser(description='Existing CS RGBA GeoTIFF -> XYZ PNG')
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True,help='New folder, must not exist')
    p.add_argument('--min-zoom',type=int,default=12)
    p.add_argument('--max-zoom',type=int,default=18)
    p.add_argument('--max-tiles',type=int,default=100000)
    a = p.parse_args()
    gdal.UseExceptions()
    write_xyz(a.input,a.output,dict(xyz_min_zoom=a.min_zoom,xyz_max_zoom=a.max_zoom,xyz_max_tiles=a.max_tiles),gdal,osr)
