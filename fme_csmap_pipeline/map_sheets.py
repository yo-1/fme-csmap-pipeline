"""Japanese national large-scale map sheets in plane rectangular coordinates.

Rule: GSI Public Survey Standard Symbols, Appendix 7, Articles 84-85.
https://www.gsi.go.jp/common/000258741.pdf (PDF pages 26-30, 1-based)
All function arguments use easting/northing (survey Y/X), in metres.
"""
from dataclasses import dataclass
import math
from pathlib import Path

SUBDIVISIONS = {5000: 1, 2500: 2, 1000: 5, 500: 10}


@dataclass(frozen=True)
class Sheet:
    code: str
    west: float
    south: float
    east: float
    north: float


def dimensions(level):
    if level not in SUBDIVISIONS:
        raise ValueError("sheet_level must be 5000, 2500, 1000 or 500")
    div = SUBDIVISIONS[level]
    return 4000 // div, 3000 // div


def sheet_at(easting, northing, zone, level):
    """West/north edges included, east/south edges excluded for point lookup."""
    if not isinstance(zone, int) or isinstance(zone, bool) or not 1 <= zone <= 19:
        raise ValueError("plane_zone must be an integer in 1..19")
    if not (-160000 <= easting < 160000 and -300000 < northing <= 300000):
        raise ValueError("Coordinate outside standard sheet range: E[-160000,160000), N(-300000,300000]")
    width, height = dimensions(level)
    col = math.floor((easting + 160000) / width)
    row = math.floor((300000 - northing) / height)
    div = SUBDIVISIONS[level]
    row5, subrow = divmod(row, div)
    col5, subcol = divmod(col, div)
    majorrow, minorrow = divmod(row5, 10)
    majorcol, minorcol = divmod(col5, 10)
    code = f"{zone:02d}{chr(65+majorrow)}{chr(65+majorcol)}{minorrow}{minorcol}"
    if level == 2500:
        code += str(subrow*2 + subcol + 1)
    elif level == 1000:
        code += f"{subrow}{chr(65+subcol)}"
    elif level == 500:
        code += f"{subrow}{subcol}"
    west = -160000 + col*width
    north = 300000 - row*height
    return Sheet(code, west, north-height, west+width, north)


def intersecting_sheets(bounds, zone, level, max_sheets=100000):
    """Yield full sheets intersecting (west,south,east,north) with positive area."""
    west, south, east, north = bounds
    if not all(math.isfinite(v) for v in bounds) or not (west < east and south < north):
        raise ValueError("Invalid bounds")
    if west < -160000 or east > 160000 or south < -300000 or north > 300000:
        raise ValueError("DEM bounds exceed standard national sheet domain; split by the correct plane zone")
    width, height = dimensions(level)
    c0 = math.floor((west+160000)/width)
    c1 = math.ceil((east+160000)/width)
    r0 = math.floor((300000-north)/height)
    r1 = math.ceil((300000-south)/height)
    if (r1-r0)*(c1-c0) > max_sheets:
        raise ValueError("Too many sheets; check extent, sheet_level and max_sheets")
    for r in range(r0, r1):
        for c in range(c0, c1):
            yield sheet_at(-160000+(c+.5)*width, 300000-(r+.5)*height, zone, level)


def pixel_window(sheet, transform, source_width, source_height):
    """Integer source window and target offsets; no interpolation or rounding away errors."""
    west, dx, rx, north, ry, negdy = transform
    if rx != 0 or ry != 0 or dx <= 0 or negdy >= 0:
        raise ValueError("Expected a north-up raster")
    values = ((sheet.west-west)/dx, (north-sheet.north)/(-negdy),
              (sheet.east-sheet.west)/dx, (sheet.north-sheet.south)/(-negdy))
    if any(abs(v-round(v)) > 1e-6 for v in values):
        raise ValueError("Sheet boundary is not aligned with raster pixel edges")
    x, y, width, height = map(lambda v: int(round(v)), values)
    sx, sy = max(x, 0), max(y, 0)
    ex, ey = min(x+width, source_width), min(y+height, source_height)
    return dict(src_x=sx, src_y=sy, width=max(0, ex-sx), height=max(0, ey-sy),
                dst_x=sx-x, dst_y=sy-y, tile_width=width, tile_height=height)


def cut_sheets(source_path, output_dir, c, gdal, ogr, osr):
    """Copy RGBA pixels to full-extent sheets; transparent padding; output GPKG index."""
    src = gdal.Open(str(source_path))
    if src.RasterCount != 4:
        raise ValueError("Expected four-band RGBA CS image")
    gt = src.GetGeoTransform()
    bounds = (gt[0], gt[3]+gt[5]*src.RasterYSize,
              gt[0]+gt[1]*src.RasterXSize, gt[3])
    out = Path(output_dir)
    out.mkdir(exist_ok=False)
    index = ogr.GetDriverByName("GPKG").CreateDataSource(str(out.parent/"sheet_index.gpkg"))
    srs = src.GetSpatialRef().Clone()
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    layer = index.CreateLayer("sheets", srs, ogr.wkbPolygon)
    for name, field_type in (("sheet_code", ogr.OFTString), ("plane_zone", ogr.OFTInteger),
                             ("sheet_level", ogr.OFTInteger), ("file", ogr.OFTString),
                             ("cell_m", ogr.OFTReal), ("valid_px", ogr.OFTInteger64)):
        layer.CreateField(ogr.FieldDefn(name, field_type))
    written, skipped = 0, 0
    block = c["block_size"]
    for sheet in intersecting_sheets(bounds, c["plane_zone"], c["sheet_level"], c["max_sheets"]):
        p = pixel_window(sheet, gt, src.RasterXSize, src.RasterYSize)
        if p["tile_width"]*p["tile_height"] > c["max_sheet_pixels"]:
            raise ValueError("Sheet exceeds max_sheet_pixels")
        # Inspect alpha blockwise to avoid producing empty files.
        count = 0
        for y in range(0, p["height"], block):
            for x in range(0, p["width"], block):
                a = src.GetRasterBand(4).ReadAsArray(p["src_x"]+x, p["src_y"]+y,
                    min(block, p["width"]-x), min(block, p["height"]-y))
                count += int((a > 0).sum())
        if count == 0:
            skipped += 1
            continue
        dest = out/f"{sheet.code}.tif"
        dst = gdal.GetDriverByName("GTiff").Create(str(dest), p["tile_width"], p["tile_height"],
            4, gdal.GDT_Byte, options=["TILED=YES", f"COMPRESS={c['compression']}",
                                      "BIGTIFF=IF_SAFER", "PHOTOMETRIC=RGB", "ALPHA=YES"])
        dst.SetProjection(src.GetProjection())
        dst.SetGeoTransform((sheet.west, gt[1], 0, sheet.north, 0, gt[5]))
        dst.SetMetadataItem("SHEET_CODE", sheet.code)
        dst.SetMetadataItem("SHEET_LEVEL", str(c["sheet_level"]))
        dst.SetMetadataItem("PLANE_ZONE", str(c["plane_zone"]))
        for b in range(1, 5):
            dst.GetRasterBand(b).SetColorInterpretation(src.GetRasterBand(b).GetColorInterpretation())
            dst.GetRasterBand(b).Fill(0)
        for y in range(0, p["height"], block):
            for x in range(0, p["width"], block):
                bw, bh = min(block, p["width"]-x), min(block, p["height"]-y)
                for b in range(1, 5):
                    a = src.GetRasterBand(b).ReadAsArray(p["src_x"]+x, p["src_y"]+y, bw, bh)
                    dst.GetRasterBand(b).WriteArray(a, p["dst_x"]+x, p["dst_y"]+y)
        dst.FlushCache()
        dst = None
        # World-file origin is the centre of the upper-left pixel, not its corner.
        tfw = (gt[1], 0, 0, gt[5], sheet.west+gt[1]/2, sheet.north+gt[5]/2)
        dest.with_suffix(".tfw").write_text("\n".join(f"{v:.12f}" for v in tfw)+"\n", encoding="ascii")
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for e, n in ((sheet.west,sheet.south), (sheet.east,sheet.south),
                     (sheet.east,sheet.north), (sheet.west,sheet.north), (sheet.west,sheet.south)):
            ring.AddPoint_2D(e, n)
        polygon = ogr.Geometry(ogr.wkbPolygon)
        polygon.AddGeometry(ring)
        feature = ogr.Feature(layer.GetLayerDefn())
        for key, value in dict(sheet_code=sheet.code, plane_zone=c["plane_zone"],
            sheet_level=c["sheet_level"], file=f"sheets/{dest.name}", cell_m=c["cell_size"], valid_px=count).items():
            feature.SetField(key, value)
        feature.SetGeometry(polygon)
        if layer.CreateFeature(feature) != 0:
            raise RuntimeError("Failed to write sheet index")
        feature = None
        written += 1
        print(f"Sheet {sheet.code}: {p['tile_width']} x {p['tile_height']}", flush=True)
    layer = index = src = None
    if written == 0:
        raise ValueError("No nonempty sheets produced")
    return dict(written=written, skipped_empty=skipped)
