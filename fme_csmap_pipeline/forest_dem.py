"""Readers for Japanese forestry airborne-laser DEM deliverables. GPL-3.0-only.

Supported products are LEM+CSV metadata pairs, XYZ/CSV regular grids,
TIFF+world-file rasters, GeoTIFF rasters, and ZIP containers of those files.
"""
import csv
import hashlib
import math
from pathlib import Path, PurePosixPath
import re
import zipfile

import numpy as np

NODATA = -999999.0
ARCHIVE_MEMBER_LIMIT = 512 * 1024 * 1024
ARCHIVE_TOTAL_LIMIT = 8 * 1024 * 1024 * 1024
SIDECARS = {'.tfw', '.tifw', '.wld', '.prj', '.aux.xml'}
PRIMARY = {'.lem', '.csv', '.txt', '.xyz', '.tif', '.tiff'}


def _decode(data):
    for encoding in ('utf-8-sig', 'cp932', 'shift_jis'):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            pass
    raise ValueError('Text encoding is neither UTF-8 nor CP932/Shift_JIS')


def _key(value):
    return re.sub(r'[\s　_()（）・]+', '', value).lower()


def read_lem_metadata(path):
    """Read the companion CSV header defined for LEM mesh elevation files."""
    text, encoding = _decode(Path(path).read_bytes())
    result = {}
    for row in csv.reader(text.splitlines()):
        if len(row) >= 2 and row[0].strip():
            result[_key(row[0])] = row[1].strip()
    aliases = {
        'nx': ('東西方向の点数', '東西方向点数'),
        'ny': ('南北方向の点数', '南北方向点数', '記録レコード数'),
        'dx': ('東西方向のデータ間隔', '東西方向データ間隔'),
        'dy': ('南北方向のデータ間隔', '南北方向データ間隔'),
        'south_n': ('区画左下x座標',),
        'west_e': ('区画左下y座標',),
        'north_n': ('区画右上x座標',),
        'east_e': ('区画右上y座標',),
        'zone': ('平面直角座標系番号',),
        'survey_year': ('測量年',),
        'sheet': ('図名',),
    }
    values = {}
    for name, names in aliases.items():
        value = next((result.get(_key(n)) for n in names if result.get(_key(n)) not in (None, '')), None)
        if value is not None:
            values[name] = value
    for required in ('nx', 'ny', 'dx', 'dy', 'south_n', 'west_e', 'north_n', 'east_e'):
        if required not in values:
            raise ValueError(f'LEM metadata CSV is missing: {required} ({path})')
    for name in ('nx', 'ny', 'zone'):
        if name in values: values[name] = int(float(values[name]))
    for name in ('dx', 'dy', 'south_n', 'west_e', 'north_n', 'east_e'):
        values[name] = float(values[name])
    values['encoding'] = encoding
    return values


def is_lem_metadata(path):
    try:
        metadata = read_lem_metadata(path)
        return metadata['nx'] > 0 and metadata['ny'] > 0
    except (OSError, UnicodeError, ValueError, csv.Error):
        return False


def _coordinate_scale(meta, requested='auto'):
    if requested != 'auto':
        scale = float(requested)
        if scale <= 0 or not math.isfinite(scale): raise ValueError('Invalid LEM coordinate scale')
        return scale
    # Standard files normally store plane rectangular X/Y as centimetre integers.
    # Select the scale whose stated bounds best match point count and spacing.
    target_n = meta['ny'] * meta['dy']; target_e = meta['nx'] * meta['dx']
    candidates = (1.0, 0.01, 0.001)
    def error(scale):
        dn = abs((meta['north_n']-meta['south_n'])*scale-target_n) / max(target_n, 1e-9)
        de = abs((meta['east_e']-meta['west_e'])*scale-target_e) / max(target_e, 1e-9)
        return dn+de
    scale = min(candidates, key=error)
    if error(scale) > 0.05:
        raise ValueError('LEM metadata bounds do not agree with point counts/spacing; set forest_lem_coordinate_scale')
    return scale


def read_lem(path, metadata_path, max_pixels, coordinate_scale='auto', feedback=None):
    """Return (array, geotransform, metadata) for a fixed-width LEM file."""
    meta = read_lem_metadata(metadata_path)
    nx, ny = meta['nx'], meta['ny']
    if nx <= 0 or ny <= 0 or nx*ny > max_pixels:
        raise ValueError('LEM grid exceeds input pixel limit')
    scale = _coordinate_scale(meta, coordinate_scale)
    raw = Path(path).read_bytes()
    text, encoding = _decode(raw)
    lines = [line.rstrip('\r\n') for line in text.splitlines() if line.strip()]
    if len(lines) != ny:
        raise ValueError(f'LEM row count mismatch: expected {ny}, found {len(lines)}')
    array = np.full((ny, nx), NODATA, dtype='float32')
    for row, line in enumerate(lines):
        if feedback is not None and row % 256 == 0 and feedback.isCanceled():
            raise RuntimeError('Cancelled during LEM conversion')
        if len(line) < 10 + nx*5:
            raise ValueError(f'LEM row {row+1} is shorter than the fixed-width definition')
        record = line[6:10].strip()
        if record and record.lstrip('+-').isdigit() and int(record) != row+1:
            raise ValueError(f'LEM record number mismatch at row {row+1}')
        for col in range(nx):
            field = line[10+col*5:15+col*5].strip()
            try: value = int(field)
            except ValueError as exc: raise ValueError(f'LEM invalid elevation at row {row+1}, column {col+1}') from exc
            if value not in (-9999, -1111): array[row, col] = value * 0.1
    west = meta['west_e']*scale; north = meta['north_n']*scale
    gt = (west, meta['dx'], 0.0, north, 0.0, -meta['dy'])
    detail = {**meta, 'coordinate_scale': scale, 'lem_encoding': encoding,
              'metadata_path': str(metadata_path), 'nodata_codes': [-9999, -1111]}
    return array, gt, detail


def _safe_member(name):
    p = PurePosixPath(name.replace('\\', '/'))
    return bool(name and not p.is_absolute() and '..' not in p.parts and not re.match(r'^[A-Za-z]:', name))


def extract_archive(path, destination):
    """Safely extract only supported DEM members and required sidecars."""
    destination = Path(destination); destination.mkdir(parents=True, exist_ok=False)
    extracted=[]; total=0
    with zipfile.ZipFile(path) as archive:
        for entry in sorted(archive.infolist(), key=lambda e:e.filename.lower()):
            if entry.is_dir(): continue
            if not _safe_member(entry.filename): raise ValueError('Unsafe ZIP member path: '+entry.filename)
            suffix = Path(entry.filename).suffix.lower()
            compound = Path(entry.filename).name.lower().endswith('.aux.xml')
            if suffix not in PRIMARY|SIDECARS and not compound: continue
            if entry.file_size > ARCHIVE_MEMBER_LIMIT: raise ValueError('ZIP member exceeds 512 MiB: '+entry.filename)
            total += entry.file_size
            if total > ARCHIVE_TOTAL_LIMIT: raise ValueError('Supported ZIP members exceed 8 GiB uncompressed')
            target = destination.joinpath(*PurePosixPath(entry.filename).parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(entry) as source, target.open('wb') as sink:
                while True:
                    block=source.read(1024*1024)
                    if not block: break
                    sink.write(block)
            extracted.append(target)
    if not extracted: raise ValueError('ZIP contains no supported forestry DEM files')
    return extracted


def expand_sources(paths, work):
    expanded=[]; archives=[]
    for ordinal, value in enumerate(paths, 1):
        path=Path(value)
        if path.suffix.lower()=='.zip':
            token=hashlib.sha256(str(path.resolve()).encode()).hexdigest()[:10]
            members=extract_archive(path, Path(work)/f'archive_{ordinal:06d}_{token}')
            expanded.extend(members); archives.append({'path':str(path),'members':len(members)})
        else:
            expanded.append(path)
            if path.suffix.lower()=='.lem':
                companion=path.with_suffix('.csv')
                if not companion.exists():
                    alternatives=[p for p in path.parent.glob('*') if p.suffix.lower()=='.csv' and p.stem.lower()==path.stem.lower()]
                    companion=alternatives[0] if alternatives else companion
                if companion.exists() and companion not in expanded: expanded.append(companion)
    return expanded, archives


def classify_sources(paths):
    """Identify primary files and avoid treating LEM metadata CSV as an XYZ grid."""
    files=[Path(p) for p in paths]
    by_key={(p.parent, p.stem.lower()):p for p in files if p.suffix.lower()=='.csv'}
    result=[]
    for p in files:
        suffix=p.suffix.lower()
        if suffix in SIDECARS or p.name.lower().endswith('.aux.xml'): continue
        if suffix=='.lem':
            header=by_key.get((p.parent,p.stem.lower()))
            if header is None or not is_lem_metadata(header):
                raise ValueError('LEM requires its companion metadata CSV with the same stem: '+str(p))
            result.append({'kind':'lem','path':p,'metadata':header})
        elif suffix=='.csv' and is_lem_metadata(p):
            if not any(q.suffix.lower()=='.lem' and q.parent==p.parent and q.stem.lower()==p.stem.lower() for q in files):
                raise ValueError('LEM metadata CSV has no companion .lem file: '+str(p))
        elif suffix in ('.csv','.txt','.xyz'):
            result.append({'kind':'grid','path':p})
        elif suffix in ('.tif','.tiff'):
            result.append({'kind':'raster','path':p})
    if not result: raise ValueError('No primary forestry DEM files were detected')
    return result


def autodetect_xyz(path, axis_order='north_east'):
    """Detect delimiter/header/columns for a conventional three-column XYZ grid."""
    text, encoding=_decode(Path(path).read_bytes())
    lines=[line for line in text.splitlines() if line.strip()]
    for skip, line in enumerate(lines[:100]):
        delimiter='comma' if ',' in line else 'tab' if '\t' in line else 'semicolon' if ';' in line else 'space'
        row=line.split() if delimiter=='space' else next(csv.reader([line],delimiter={'comma':',','tab':'\t','semicolon':';'}[delimiter]))
        numeric=[]
        for index, value in enumerate(row):
            try: float(value); numeric.append(index+1)
            except ValueError: pass
        if len(numeric)>=3:
            return dict(text_encoding=encoding,text_delimiter=delimiter,text_skip_rows=skip,
                        text_columns=numeric[:3],text_axis_order=axis_order)
    raise ValueError('Could not detect three numeric XYZ columns: '+str(path))
