#!/usr/bin/env python3
"""DEM mosaic -> independent CS-style RGBA relief -> Japanese map-sheet GeoTIFFs.

Read README_ja.md. This is NOT a pixel-identical reproduction of Nagano CSMap.
Run with a dedicated GDAL Python environment, from FME SystemCaller.
Copyright (C) 2026 Yoichi Wada. GPL-3.0-only.
"""
import argparse
import glob
import json
import math
from pathlib import Path
import sys
import traceback

import numpy as np
from scipy.ndimage import gaussian_filter, minimum_filter
from map_sheets import dimensions, cut_sheets, intersecting_sheets

VERSION = "0.8.1"
from xyz_tiles import DEFAULTS as XYZ_DEFAULTS, validate_xyz, write_xyz

from input_sources import DEFAULTS as INPUT_DEFAULTS, validate_input, discover, prepare_inputs

NODATA = -999999.0
COLOR_DEFAULTS = dict(
    valley_rgb=[74, 140, 198], ridge_rgb=[192, 119, 71], neutral_rgb=[248, 247, 242],
    elevation_low_rgb=[222, 236, 244], elevation_high_rgb=[250, 242, 223],
    curvature_strength=1.0, elevation_mix=0.15, slope_darkness=0.55,
    brightness=1.0, saturation=1.0, contrast=1.0, gamma=1.0)

FME_COLOR_DEFAULTS = dict(
    elevation_low_rgb=[36, 36, 36], elevation_high_rgb=[246, 246, 246],
    valley_rgb=[50, 96, 207], neutral_rgb=[255, 254, 190], ridge_rgb=[198, 72, 59])
FME_SLOPE_RED = ([247, 213, 213], [134, 28, 33])
FME_SLOPE_BW = ([246, 246, 246], [36, 36, 36])
FME_CURV_BLUE = ([42, 95, 131], [208, 223, 230])
FME_STRETCH = ((65, 234), (57, 229), (66, 216))
STRETCH_MODES = ('none', 'nagano_reference', 'custom')

# FMEマニュアル別紙3の明記値（曲率±10 <-> 本プラグイン±0.1 1/mはFMEワークスペースの
# RasterConvolverカーネル除数 cell_size**2*0.01 から確認済み。VALIDATION.txt参照）。
FME_MANUAL_CURVATURE_LIMIT = 0.1
# 暫定値（FMEマニュアルの値ではない）。
# FME標準値(上記)を上書きしてはならず、別プロファイル(config.forestry_tuned.json)としてのみ使用する。
FORESTRY_TUNED_CURVATURE_LIMIT = 0.03


def color_settings(settings=None):
    if settings is None:
        settings = {}
    if not isinstance(settings, dict):
        raise ValueError("color must be an object")
    if set(settings)-set(COLOR_DEFAULTS):
        raise ValueError(f"Unknown color keys: {sorted(set(settings)-set(COLOR_DEFAULTS))}")
    result = {**COLOR_DEFAULTS, **settings}
    for key, value in result.items():
        if key.endswith('_rgb'):
            if not isinstance(value, (list, tuple)) or len(value) != 3 or any(
                not isinstance(v, int) or isinstance(v, bool) or not 0 <= v <= 255 for v in value):
                raise ValueError(f"color.{key} requires three integer RGB values in 0..255")
        else:
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
                raise ValueError(f"color.{key} must be finite")
            low, high = ((0, 1) if key in ('curvature_strength','elevation_mix','slope_darkness')
                         else (0.1, 5) if key == 'gamma' else (0, 3))
            if not low <= value <= high:
                raise ValueError(f"color.{key} must be within {low}..{high}")
    return result


def fme_color_settings(settings=None):
    settings = settings or {}
    if not isinstance(settings, dict):
        raise ValueError('color must be an object')
    unknown = set(settings) - set(FME_COLOR_DEFAULTS)
    if unknown:
        raise ValueError(f'Unknown FME color keys: {sorted(unknown)}')
    result = {**FME_COLOR_DEFAULTS, **settings}
    for key, value in result.items():
        if not isinstance(value, (list, tuple)) or len(value) != 3 or any(
                not isinstance(v, int) or isinstance(v, bool) or not 0 <= v <= 255 for v in value):
            raise ValueError(f'color.{key} requires three integer RGB values in 0..255')
    return result


def _lerp(lo, hi, t):
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    return lo + t[..., None] * (hi - lo)


def validate_stretch(mode, custom_range):
    """Resolve stretch_mode/stretch_range into None or 3 (lo, hi) tuples per band.

    none: no stretch. nagano_reference: FME manual's Nagano-empirical values
    (this is a documented regional result, not a universal spec — see README).
    custom: caller-supplied per-band [low, high] pairs.
    """
    if mode not in STRETCH_MODES:
        raise ValueError(f"stretch_mode must be one of {STRETCH_MODES}")
    if mode == 'none':
        return None
    if mode == 'nagano_reference':
        return FME_STRETCH
    if (not isinstance(custom_range, (list, tuple)) or len(custom_range) != 3
            or any(not isinstance(pair, (list, tuple)) or len(pair) != 2 for pair in custom_range)):
        raise ValueError("stretch_range must be three [low, high] pairs, one per R/G/B band")
    resolved = []
    for lo, hi in custom_range:
        if (not isinstance(lo, (int, float)) or isinstance(lo, bool)
                or not isinstance(hi, (int, float)) or isinstance(hi, bool)
                or not math.isfinite(lo) or not math.isfinite(hi)):
            raise ValueError("stretch_range values must be finite numbers")
        if not 0 <= lo < hi <= 255:
            raise ValueError("stretch_range requires 0 <= low < high <= 255 for each band")
        resolved.append((float(lo), float(hi)))
    return tuple(resolved)


def _fme_relief(raw, slope, curvature, curvature_limit, slope_max, elev_range, color,
                stretch_mode='nagano_reference', stretch_range=None):
    """FME manual palette, five independent layers, weighted blend and optional per-channel stretch."""
    palette = fme_color_settings(color)
    stretch = validate_stretch(stretch_mode, stretch_range)
    emin, emax = elev_range
    height_t = np.clip((raw - emin) / (emax - emin), 0, 1)
    slope_t = np.clip(slope / slope_max, 0, 1)
    curve_t = np.clip((curvature + curvature_limit) / (2 * curvature_limit), 0, 1)
    height = _lerp(palette['elevation_low_rgb'], palette['elevation_high_rgb'], height_t)
    slope_red = _lerp(*FME_SLOPE_RED, slope_t)
    slope_bw = _lerp(*FME_SLOPE_BW, slope_t)
    curve_blue = _lerp(*FME_CURV_BLUE, curve_t)
    left = _lerp(palette['valley_rgb'], palette['neutral_rgb'], np.minimum(curve_t * 2, 1))
    right = _lerp(palette['neutral_rgb'], palette['ridge_rgb'], np.maximum(curve_t * 2 - 1, 0))
    curve_ryb = np.where((curve_t < .5)[..., None], left, right)
    rgb = height * .125 + slope_red * .25 + slope_bw * .25 + curve_blue * .125 + curve_ryb * .25
    if stretch is not None:
        for band, (lo, hi) in enumerate(stretch):
            rgb[..., band] = (rgb[..., band] - lo) / (hi - lo) * 255
    return np.clip(rgb, 0, 255)


SLOPE_ALGORITHMS = ('horn', 'central_difference')
SLOPE_DESCRIPTIONS = {
    'horn': 'Horn (1981) 3x3 weighted differences of unsmoothed elevation (degrees)',
    'central_difference': '2-point central differences of unsmoothed elevation (degrees)',
}


def slope_gradients(raw, cell, slope_algorithm='horn'):
    """Return (dz/dx, dz/dy) for the selected slope algorithm.

    Both need only a 1-cell neighbourhood, which the existing
    ceil(4*sigma_m/cell)+1 halo always covers. Only the magnitude is used
    downstream (np.hypot), so the sign convention does not affect results.
    'central_difference' reproduces v0.7.2 and earlier exactly. Same
    formulas as the QGIS plugin csmap-sheets v0.10.0.
    """
    if slope_algorithm not in SLOPE_ALGORITHMS:
        raise ValueError(f"slope_algorithm must be one of {SLOPE_ALGORITHMS}")
    if slope_algorithm == 'central_difference':
        dx = (np.roll(raw, -1, 1) - np.roll(raw, 1, 1)) / (2 * cell)
        dy = (np.roll(raw, -1, 0) - np.roll(raw, 1, 0)) / (2 * cell)
        return dx, dy
    n = np.roll(raw, 1, 0); s = np.roll(raw, -1, 0)
    w = np.roll(raw, 1, 1); e = np.roll(raw, -1, 1)
    nw = np.roll(n, 1, 1); ne = np.roll(n, -1, 1)
    sw = np.roll(s, 1, 1); se = np.roll(s, -1, 1)
    dx = ((ne + 2 * e + se) - (nw + 2 * w + sw)) / (8 * cell)
    dy = ((sw + 2 * s + se) - (nw + 2 * n + ne)) / (8 * cell)
    return dx, dy


SIGMA_UNITS = ('m', 'px')


def gaussian_sigma(c):
    """Return (sigma in calculation-grid pixels, kernel radius in pixels) for a config.

    v0.9.0: the Gaussian sigma can be given as a ground distance (``sigma_unit='m'``,
    ``sigma_m``) or as a number of calculation-grid pixels (``sigma_unit='px'``,
    ``sigma_px``), the same as csmap-sheets v0.12.0. The calculation grid is the
    reprojected DEM grid of ``cell_size`` (always square) on which relief() smooths.
    The 'm' branch keeps the exact float expressions used up to v0.8.1.
    """
    if c.get('sigma_unit', 'm') == 'px':
        sigma_px = c['sigma_px']
        return sigma_px, math.ceil(4 * sigma_px)
    return c['sigma_m'] / c['cell_size'], math.ceil(4 * c['sigma_m'] / c['cell_size'])


def smoothing_record(c):
    """Smoothing conditions for logs and run.json (setting and effective values kept apart)."""
    sigma_px, radius = gaussian_sigma(c)
    unit = c.get('sigma_unit', 'm')
    return {
        'gaussian_sigma_m': sigma_px * c['cell_size'] if unit == 'px' else c['sigma_m'],
        'sigma_unit': unit,
        'sigma_setting': c['sigma_px'] if unit == 'px' else c['sigma_m'],
        'cell_size_m': c['cell_size'],
        'effective_sigma_m': sigma_px * c['cell_size'] if unit == 'px' else c['sigma_m'],
        'effective_sigma_px': sigma_px,
        'kernel_radius_px': radius if sigma_px > 0 else 0,
        'kernel_size_px': 2 * radius + 1 if sigma_px > 0 else 0,
    }


def smoothing_summary(c):
    r = smoothing_record(c)
    if r['effective_sigma_px'] == 0:
        return f"Smoothing: none (sigma=0, unit {r['sigma_unit']})"
    return (f"Smoothing: sigma {r['sigma_setting']:g} {r['sigma_unit']} on a "
            f"{r['cell_size_m']:g} m grid = {r['effective_sigma_m']:g} m = "
            f"{r['effective_sigma_px']:g} px; kernel {r['kernel_size_px']}x{r['kernel_size_px']} px")


def relief(z, valid, cell, sigma_m, curvature_limit, slope_max, elev_range, color=None,
           color_model='legacy', stretch_mode='nagano_reference', stretch_range=None,
           slope_algorithm='horn', sigma_px=None):
    """Return RGBA, slope degrees and negative-Laplacian proxy (1/m).

    All samples touching a missing value within the full processing support
    are transparent. Callers must supply ceil(4*sigma_m/cell)+1 halo cells.
    ``sigma_px`` (v0.9.0) gives sigma directly in grid pixels and then takes
    precedence over ``sigma_m``; the halo is ceil(4*sigma_px)+1 in that case.
    """
    if color_model not in ('legacy', 'fme'):
        raise ValueError("color_model must be legacy or fme")
    if slope_algorithm not in SLOPE_ALGORITHMS:
        raise ValueError(f"slope_algorithm must be one of {SLOPE_ALGORITHMS}")
    tone = color_settings(color) if color_model == 'legacy' else None
    if sigma_px is None:
        radius = math.ceil(4 * sigma_m / cell)
        sigma_grid = sigma_m / cell
    else:
        radius = math.ceil(4 * sigma_px)
        sigma_grid = sigma_px
    halo = radius + 1
    valid = valid & np.isfinite(z)
    raw = np.where(valid, z, 0).astype(np.float64)
    smooth = (gaussian_filter(raw, sigma=sigma_grid, radius=radius,
                              mode="constant", cval=0)
              if sigma_grid > 0 else raw.copy())
    dx, dy = slope_gradients(raw, cell, slope_algorithm)
    slope = np.degrees(np.arctan(np.hypot(dx, dy)))
    curvature = -(np.roll(smooth, -1, 1) + np.roll(smooth, 1, 1)
                  + np.roll(smooth, -1, 0) + np.roll(smooth, 1, 0)
                  - 4 * smooth) / cell**2
    safe = minimum_filter(valid.astype(np.uint8), size=2 * halo + 1,
                          mode="constant", cval=0).astype(bool)
    if color_model == 'fme':
        rgb = _fme_relief(raw, slope, curvature, curvature_limit, slope_max, elev_range, color,
                         stretch_mode, stretch_range)
    else:
        t = np.clip(curvature / curvature_limit, -1, 1)
        neutral = np.array(tone['neutral_rgb'], dtype=float)
        warm = np.array(tone['ridge_rgb'], dtype=float)
        cool = np.array(tone['valley_rgb'], dtype=float)
        target = np.where((t >= 0)[..., None], warm, cool)
        rgb = neutral + tone['curvature_strength'] * np.abs(t)[..., None] * (target - neutral)
        emin, emax = elev_range
        h = np.clip((raw - emin) / (emax - emin), 0, 1)
        low = np.array(tone['elevation_low_rgb'], dtype=float)
        high = np.array(tone['elevation_high_rgb'], dtype=float)
        tint = low + h[..., None] * (high-low)
        mix = tone['elevation_mix']
        rgb = (1-mix) * rgb + mix * tint
        shade = 1 - tone['slope_darkness'] * np.clip(slope / slope_max, 0, 1)
        rgb *= shade[..., None]
        if tone['saturation'] != 1:
            gray = np.sum(rgb * np.array([.2126, .7152, .0722]), axis=2, keepdims=True)
            rgb = gray + tone['saturation']*(rgb-gray)
        if tone['contrast'] != 1:
            rgb = 127.5 + tone['contrast']*(rgb-127.5)
        rgb = np.clip(rgb*tone['brightness'], 0, 255)
        if tone['gamma'] != 1:
            rgb = 255 * (rgb/255)**(1/tone['gamma'])
    rgb = np.rint(rgb).astype(np.uint8)
    rgb[~safe] = 0
    rgba = np.concatenate([rgb, (safe.astype(np.uint8) * 255)[..., None]], axis=2)
    return rgba, slope, curvature


def rendering_settings(c):
    """Build the complete, serialisable rendering record for run.json.

    settings=c alone omits the FME palette/weights that are module constants
    (not user-configurable in this CLI build), so recompute them here.
    """
    if c.get('color_model') == 'fme':
        palette = fme_color_settings(c.get('color'))
        stretch_mode = c.get('stretch_mode', 'nagano_reference')
        stretch_range = c.get('stretch_range')
        resolved = validate_stretch(stretch_mode, stretch_range)
        return {
            'mode_id': 'fme_manual',
            'normalization': {
                'elevation_m': list(c['elevation_range']),
                'slope_degrees': [0.0, c['slope_max']],
                'curvature_1_per_m': [-c['curvature_limit'], c['curvature_limit']],
            },
            'colors': {
                'elevation_low': palette['elevation_low_rgb'], 'elevation_high': palette['elevation_high_rgb'],
                'slope_a_low': list(FME_SLOPE_RED[0]), 'slope_a_high': list(FME_SLOPE_RED[1]),
                'slope_b_low': list(FME_SLOPE_BW[0]), 'slope_b_high': list(FME_SLOPE_BW[1]),
                'curvature_a_low': list(FME_CURV_BLUE[0]), 'curvature_a_high': list(FME_CURV_BLUE[1]),
                'curvature_b_low': palette['valley_rgb'], 'curvature_b_mid': palette['neutral_rgb'],
                'curvature_b_high': palette['ridge_rgb'],
            },
            'weights': {'elevation': 0.125, 'slope_a': 0.25, 'slope_b': 0.25,
                       'curvature_a': 0.125, 'curvature_b': 0.25},
            'stretch': {
                'mode': stretch_mode,
                'ranges_rgb': [list(pair) for pair in resolved] if resolved is not None else None,
                'nagano_reference_is_optional_empirical_value': True,
            },
            'curvature_sign': 'negative=concave/valley/blue; zero=yellow in curvature_b; positive=convex/ridge/red',
            'smoothing': smoothing_record(c) if 'cell_size' in c else {'gaussian_sigma_m': c['sigma_m']},
            'terrain_calculation': {
                'slope': SLOPE_DESCRIPTIONS[c.get('slope_algorithm', 'horn')],
                'slope_algorithm': c.get('slope_algorithm', 'horn'),
                'curvature': 'negative Laplacian of Gaussian-smoothed elevation (1/m); FME workspace confirmed same sign/kernel shape (VALIDATION.txt)',
            },
        }
    tone = color_settings(c.get('color'))
    return {
        'mode_id': 'independent_v040',
        'normalization': {
            'elevation_m': list(c['elevation_range']),
            'slope_degrees': [0.0, c['slope_max']],
            'curvature_1_per_m': [-c['curvature_limit'], c['curvature_limit']],
        },
        'colors': tone,
        'smoothing': smoothing_record(c) if 'cell_size' in c else {'gaussian_sigma_m': c['sigma_m']},
        'terrain_calculation': {
            'slope': SLOPE_DESCRIPTIONS[c.get('slope_algorithm', 'horn')],
            'slope_algorithm': c.get('slope_algorithm', 'horn'),
        },
    }


def read_config(path, overrides=None):
    c = json.loads(path.read_text(encoding="utf-8-sig"))
    c.update(overrides or {})
    provided = set(c)
    # 裸のデフォルトはindependent_v040相当（従来の独自方式）の値。FMEマニュアル値(0.1)や
    # 暫定値(0.03)をここに混在させない。color_model=fmeでcurvature_limit省略時のみ、
    # 下でFMEマニュアル公式値へ切り替える（バグ修正: 過去にここが0.03のままだったため、
    # 独自方式でcurvature_limitを省略した外部configが誤って暫定値相当になっていた）。
    defaults = dict(sigma_m=3.0, curvature_limit=0.05, slope_max=60.0,
                    elevation_range=[0, 3000], block_size=512, sheet_level=5000,
                    max_sheets=100000, max_sheet_pixels=100000000, compression="DEFLATE",
                    max_pixels=1000000000, source_nodata=None,
                    confirm_elevation_metres=False, color={}, color_model='legacy',
                    stretch_mode='nagano_reference', stretch_range=None,
                    # v0.8.0: 傾斜計算方式（ユーザー決定、2026-10-02）。既定Horn法。
                    # 中央差分法はv0.7.2以前の出力の再現用。
                    slope_algorithm='horn',
                    # v0.9.0: Gaussian σの指定方式（csmap-sheets v0.12.0と同じ仕様、
                    # 2026-10-03）。'm'＝地上距離（sigma_m。従来どおり・既定）、'px'＝計算格子の
                    # 画素数（sigma_px）。方式ごとの値を別々に保持し、読み替えない。
                    sigma_unit='m', sigma_px=3.0)
    defaults.update(XYZ_DEFAULTS)
    defaults.update(INPUT_DEFAULTS)
    unknown = set(c) - set(defaults) - {"inputs", "output_dir", "target_crs", "cell_size", "plane_zone"}
    if unknown:
        raise ValueError(f"Unknown configuration keys: {sorted(unknown)}")
    c = {**defaults, **c}
    validate_input(c)
    validate_xyz(c)
    if c['color_model'] not in ('legacy', 'fme'):
        raise ValueError('color_model must be legacy or fme')
    if c['slope_algorithm'] not in SLOPE_ALGORITHMS:
        raise ValueError(f"slope_algorithm must be one of {SLOPE_ALGORITHMS}")
    if 'slope_algorithm' not in provided:
        # v0.7.2以前の設定は中央差分法で作られていたため、既定のHorn法が黙って
        # 適用されると色の変化に気づきにくい。再現方法とあわせて警告する。
        print("WARNING: slope_algorithm is not set in this configuration; using the default "
              "'horn'. To reproduce output from v0.7.2 or earlier, set "
              "\"slope_algorithm\": \"central_difference\".", flush=True)
    if c['color_model'] == 'fme' and 'curvature_limit' not in provided:
        c['curvature_limit'] = FME_MANUAL_CURVATURE_LIMIT
    validate_stretch(c['stretch_mode'], c['stretch_range'])
    c['color'] = (color_settings(c['color']) if c['color_model'] == 'legacy'
                  else fme_color_settings(c['color']))
    for name in ("inputs", "output_dir", "target_crs", "cell_size"):
        if name not in c:
            raise ValueError(f"Missing configuration: {name}")
    if not c["confirm_elevation_metres"]:
        raise ValueError("Confirm all source elevations use metres and the same vertical datum; then set confirm_elevation_metres=true")
    for key in ("cell_size", "curvature_limit", "slope_max"):
        if not math.isfinite(c[key]) or c[key] <= 0:
            raise ValueError(f"{key} must be finite and positive")
    if c["sigma_unit"] not in SIGMA_UNITS:
        raise ValueError(f"sigma_unit must be one of {SIGMA_UNITS}")
    for key in ("sigma_m", "sigma_px"):
        if isinstance(c[key], bool) or not isinstance(c[key], (int, float)) \
                or not math.isfinite(c[key]) or c[key] < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    if not 0 < c["slope_max"] <= 90:
        raise ValueError("slope_max must be <= 90 degrees")
    er = c["elevation_range"]
    if len(er) != 2 or not all(math.isfinite(x) for x in er) or er[0] >= er[1]:
        raise ValueError("elevation_range must contain two increasing finite values")
    for key in ("block_size", "max_pixels", "max_sheets", "max_sheet_pixels"):
        if not isinstance(c[key], int) or c[key] <= 0:
            raise ValueError(f"{key} must be a positive integer")
    if not 64 <= c["block_size"] <= 2048:
        raise ValueError("block_size:64..2048")
    if c.get("plane_zone") is not None and (type(c["plane_zone"]) is not int or not 1 <= c["plane_zone"] <= 19):
        raise ValueError("plane_zone must be omitted/auto or an integer in 1..19")
    width, height = dimensions(c["sheet_level"])
    if any(abs(v/c["cell_size"]-round(v/c["cell_size"])) > 1e-6 for v in (width, height)):
        raise ValueError("cell_size must divide sheet width and height exactly")
    if round(width/c["cell_size"])*round(height/c["cell_size"]) > c["max_sheet_pixels"]:
        raise ValueError("Sheet dimensions exceed max_sheet_pixels")
    if c["compression"] not in ("DEFLATE", "NONE"):
        raise ValueError("compression must be DEFLATE or NONE")
    if gaussian_sigma(c)[1] > 512:
        raise ValueError("Gaussian radius >512 cells; reduce sigma or use coarser DEM")
    nd = c["source_nodata"]
    if nd is not None and not math.isfinite(nd):
        raise ValueError("source_nodata override must be a finite number")
    files = discover(c["inputs"], path.parent, c)
    p = Path(c["output_dir"]).expanduser()
    out = (path.parent / p if not p.is_absolute() else p).resolve()
    if out.exists():
        raise ValueError(f"Output directory already exists; choose a new directory: {out}")
    c["inputs"] = files
    c["output_dir"] = str(out)
    return c


def check_sources(c, gdal, osr):
    target = osr.SpatialReference()
    target.SetFromUserInput(c["target_crs"])
    if not target.IsProjected() or abs(target.GetLinearUnits() - 1) > 1e-9:
        raise ValueError("target_crs must be a projected CRS with metre units")
    if target.GetAuthorityCode(None) in ("3857", "3395"):
        raise ValueError("Use a suitable local projected CRS, not Web/World Mercator, for terrain derivatives")
    if target.IsCompound() or target.GetAttrValue("PROJECTION") != "Transverse_Mercator":
        raise ValueError("target_crs must be a two-dimensional Japan Plane Rectangular CRS")
    # GSI Notice: https://www.gsi.go.jp/LAW/heimencho.html
    origins = [(33,129.5), (33,131), (36,132+10/60), (33,133.5), (36,134+20/60),
               (36,136), (36,137+10/60), (36,138.5), (36,139+50/60), (40,140+50/60),
               (44,140.25), (44,142.25), (44,144.25), (26,142), (26,127.5),
               (26,124), (26,131), (20,136), (26,154)]
    zone = None
    for number, (lat0, lon0) in enumerate(origins, 1):
        if all(abs(target.GetProjParm(parm)-expected) <= 1e-8 for parm, expected in (
                ("latitude_of_origin",lat0),("central_meridian",lon0),
                ("scale_factor",.9999),("false_easting",0),("false_northing",0))):
            zone=number;break
    if zone is None:raise ValueError("Selected CRS is not one of Japan Plane Rectangular zones I-XIX")
    if c.get("plane_zone") is not None and c["plane_zone"] != zone:
        raise ValueError(f"target_crs is zone {zone}, but legacy plane_zone={c['plane_zone']}")
    c["plane_zone"] = zone
    lat, lon = origins[zone-1]
    for parm, expected in (("latitude_of_origin", lat), ("central_meridian", lon),
                           ("scale_factor", .9999), ("false_easting", 0), ("false_northing", 0)):
        if abs(target.GetProjParm(parm)-expected) > 1e-8:
            raise ValueError(f"target_crs does not match automatically detected plane zone {zone}: {parm}")
    base_srs, base_gt = None, None
    records = []
    for path in c["inputs"]:
        ds = gdal.Open(path)
        if ds is None or ds.RasterCount != 1:
            raise ValueError(f"Expected a single-band DEM: {path}")
        srs = ds.GetSpatialRef()
        gt = ds.GetGeoTransform()
        if srs is None or not ds.GetProjection():
            raise ValueError(f"Missing CRS: {path}")
        if srs.IsCompound():
            raise ValueError(f"Normalize compound/vertical CRS explicitly before input: {path}")
        if abs(gt[2]) > 1e-12 or abs(gt[4]) > 1e-12 or gt[1] <= 0 or gt[5] >= 0:
            raise ValueError(f"Rotated/south-up raster: normalize first: {path}")
        if base_srs is not None:
            if not srs.IsSame(base_srs):
                raise ValueError("Mixed source CRSs: reproject onto one common grid before this pipeline")
            if not np.allclose([gt[1], gt[5]], [base_gt[1], base_gt[5]], rtol=1e-7, atol=0):
                raise ValueError("Mixed source resolutions: normalize onto one common grid first")
            offsets = [(gt[0]-base_gt[0])/gt[1], (gt[3]-base_gt[3])/gt[5]]
            if any(abs(v-round(v)) > 1e-5 for v in offsets):
                raise ValueError(f"Source grids are not aligned: {path}")
        else:
            base_srs, base_gt = srs.Clone(), gt
        band = ds.GetRasterBand(1)
        if band.GetScale() not in (None, 1) or band.GetOffset() not in (None, 0):
            raise ValueError(f"Apply band scale/offset to actual elevation values first: {path}")
        if band.GetNoDataValue() is None and c["source_nodata"] is None:
            raise ValueError(f"No NoData metadata: set source_nodata after checking data: {path}")
        stat = Path(path).stat()
        records.append(dict(path=path, size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                            nodata=repr(band.GetNoDataValue()), width=ds.RasterXSize,
                            height=ds.RasterYSize, crs=ds.GetProjection(), transform=gt))
        ds = None
    return records


def make_relief(dem_path, output_path, c, gdal):
    src = gdal.Open(str(dem_path))
    w, h = src.RasterXSize, src.RasterYSize
    dst = gdal.GetDriverByName("GTiff").Create(str(output_path), w, h, 4,
        gdal.GDT_Byte, options=["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=IF_SAFER",
                               "PHOTOMETRIC=RGB", "ALPHA=YES"])
    dst.SetProjection(src.GetProjection())
    dst.SetGeoTransform(src.GetGeoTransform())
    for i, ci in enumerate((gdal.GCI_RedBand, gdal.GCI_GreenBand,
                             gdal.GCI_BlueBand, gdal.GCI_AlphaBand), 1):
        dst.GetRasterBand(i).SetColorInterpretation(ci)
    dst.SetMetadataItem("METHOD", "Independent CS-style; negative Laplacian, raw "
                        + {"horn": "Horn", "central_difference": "central"}[c.get("slope_algorithm", "horn")]
                        + " slope, elevation tint")
    dst.SetMetadataItem("SETTINGS", json.dumps(c, ensure_ascii=True))
    sigma_px, radius = gaussian_sigma(c)
    halo = radius+1
    block = c["block_size"]
    band = src.GetRasterBand(1)
    valid_count = 0
    for y in range(0, h, block):
        for x in range(0, w, block):
            bw, bh = min(block, w-x), min(block, h-y)
            x0, y0 = max(0, x-halo), max(0, y-halo)
            x1, y1 = min(w, x+bw+halo), min(h, y+bh+halo)
            a = band.ReadAsArray(x0, y0, x1-x0, y1-y0).astype(np.float64)
            valid = np.isfinite(a) & (a != NODATA)
            rgba, _, _ = relief(a, valid, c["cell_size"], c["sigma_m"],
                c["curvature_limit"], c["slope_max"], c["elevation_range"], c.get('color'),
                c.get('color_model', 'legacy'), c.get('stretch_mode', 'nagano_reference'),
                c.get('stretch_range'), c.get('slope_algorithm', 'horn'),
                sigma_px if c.get('sigma_unit', 'm') == 'px' else None)
            tile = rgba[y-y0:y-y0+bh, x-x0:x-x0+bw]
            valid_count += int(np.count_nonzero(tile[:, :, 3]))
            for b in range(4):
                dst.GetRasterBand(b+1).WriteArray(tile[:, :, b], x, y)
        print(f"CS relief: {min(y+block,h)}/{h} rows", flush=True)
    dst.FlushCache()
    dst = src = None
    if valid_count == 0:
        raise ValueError("No valid relief pixels: check DEM coverage, NoData and smoothing radius")
    return valid_count


def render_sheets(projected_path, out, c, gdal, ogr, osr):
    """Render a virtual province mosaic one national map sheet at a time."""
    src=gdal.Open(str(projected_path));gt=src.GetGeoTransform()
    bounds=(gt[0],gt[3]+gt[5]*src.RasterYSize,gt[0]+gt[1]*src.RasterXSize,gt[3])
    candidates=list(intersecting_sheets(bounds,c['plane_zone'],c['sheet_level'],c['max_sheets']))
    sheets_dir=out/'sheets';work=out/'sheet_work';sheets_dir.mkdir();work.mkdir()
    index=ogr.GetDriverByName('GPKG').CreateDataSource(str(out/'sheet_index.gpkg'))
    srs=src.GetSpatialRef().Clone();srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    layer=index.CreateLayer('sheets',srs,ogr.wkbPolygon)
    for name,kind in (('sheet_code',ogr.OFTString),('plane_zone',ogr.OFTInteger),
        ('sheet_level',ogr.OFTInteger),('file',ogr.OFTString),('cell_m',ogr.OFTReal),
        ('valid_px',ogr.OFTInteger64)):layer.CreateField(ogr.FieldDefn(name,kind))
    width_m,height_m=dimensions(c['sheet_level']);width=round(width_m/c['cell_size']);height=round(height_m/c['cell_size'])
    halo=gaussian_sigma(c)[1]+1
    written=[];skipped=0;valid_total=0
    for number,sheet in enumerate(candidates,1):
        dem_path=work/(sheet.code+'_dem.tif');relief_path=work/(sheet.code+'_rgba.tif')
        buffered=(sheet.west-halo*c['cell_size'],sheet.south-halo*c['cell_size'],
                  sheet.east+halo*c['cell_size'],sheet.north+halo*c['cell_size'])
        dem=gdal.Warp(str(dem_path),src,format='GTiff',dstSRS=c['target_crs'],outputBounds=buffered,
            width=width+2*halo,height=height+2*halo,outputType=gdal.GDT_Float32,
            resampleAlg='bilinear',srcNodata=NODATA,dstNodata=NODATA,errorThreshold=0,
            creationOptions=['TILED=YES','COMPRESS=DEFLATE','PREDICTOR=3','BIGTIFF=IF_SAFER'])
        if dem is None:raise RuntimeError('Map-sheet DEM warp failed: '+sheet.code)
        band=dem.GetRasterBand(1);coverage=0;block=c['block_size']
        for y in range(halo,halo+height,block):
            for x in range(halo,halo+width,block):
                a=band.ReadAsArray(x,y,min(block,halo+width-x),min(block,halo+height-y))
                coverage+=int(np.count_nonzero(np.isfinite(a)&(a!=NODATA)))
        dem=None
        if coverage==0:dem_path.unlink(missing_ok=True);skipped+=1;continue
        make_relief(dem_path,relief_path,c,gdal);dest=sheets_dir/(sheet.code+'.tif');rgba=gdal.Open(str(relief_path))
        tile=gdal.Translate(str(dest),rgba,srcWin=[halo,halo,width,height],format='GTiff',
            creationOptions=['TILED=YES',f"COMPRESS={c['compression']}",'BIGTIFF=IF_SAFER','PHOTOMETRIC=RGB','ALPHA=YES'])
        if tile is None:raise RuntimeError('Map-sheet CS write failed: '+sheet.code)
        tile.SetMetadataItem('SHEET_CODE',sheet.code);tile.SetMetadataItem('SHEET_LEVEL',str(c['sheet_level']))
        tile.SetMetadataItem('PLANE_ZONE',str(c['plane_zone']));tile.FlushCache();tile=None;rgba=None
        ds=gdal.Open(str(dest));alpha=ds.GetRasterBand(4);count=0
        for y in range(0,height,block):
            for x in range(0,width,block):count+=int(np.count_nonzero(alpha.ReadAsArray(x,y,min(block,width-x),min(block,height-y))))
        ds=None;dem_path.unlink(missing_ok=True);relief_path.unlink(missing_ok=True)
        if count==0:dest.unlink(missing_ok=True);skipped+=1;continue
        tfw=(c['cell_size'],0,0,-c['cell_size'],sheet.west+c['cell_size']/2,sheet.north-c['cell_size']/2)
        dest.with_suffix('.tfw').write_text('\n'.join(f'{v:.12f}' for v in tfw)+'\n',encoding='ascii')
        ring=ogr.Geometry(ogr.wkbLinearRing)
        for e,n in ((sheet.west,sheet.south),(sheet.east,sheet.south),(sheet.east,sheet.north),(sheet.west,sheet.north),(sheet.west,sheet.south)):ring.AddPoint_2D(e,n)
        polygon=ogr.Geometry(ogr.wkbPolygon);polygon.AddGeometry(ring);feature=ogr.Feature(layer.GetLayerDefn())
        for key,value in dict(sheet_code=sheet.code,plane_zone=c['plane_zone'],sheet_level=c['sheet_level'],
            file='sheets/'+dest.name,cell_m=c['cell_size'],valid_px=count).items():feature.SetField(key,value)
        feature.SetGeometry(polygon)
        if layer.CreateFeature(feature)!=0:raise RuntimeError('Failed to write sheet index')
        feature=None;written.append(str(dest));valid_total+=count
        print(f'Sheet {number}/{len(candidates)}: {sheet.code}',flush=True)
    layer=index=src=None
    try:work.rmdir()
    except OSError:pass
    if not written:raise ValueError('No nonempty sheets produced')
    mosaic=out/'cs_relief.vrt';vrt=gdal.BuildVRT(str(mosaic),written,resolution='highest',strict=True)
    if vrt is None:raise RuntimeError('CS sheet VRT build failed')
    vrt.FlushCache();vrt=None
    return str(mosaic),dict(written=len(written),skipped_empty=skipped),valid_total


def run(c):
    from osgeo import gdal, osr, ogr
    import scipy
    gdal.UseExceptions()
    osr.UseExceptions()
    # API floor; no gdal2tiles dependency in the map-sheet workflow.
    ver = int(gdal.VersionInfo("VERSION_NUM"))
    if ver < 3060000:
        raise RuntimeError("GDAL >=3.6 is required. The supplied environment requires GDAL >=3.10.")
    gdal.SetCacheMax(256*1024*1024)
    # Validate the target projection before doing any input conversion.
    check_sources({**c, "inputs": []}, gdal, osr)
    records = []
    out = Path(c["output_dir"])
    out.mkdir(parents=True, exist_ok=False)
    manifest = dict(status="running", algorithm_version=VERSION, settings=c,
        rendering=rendering_settings(c),
        sources=records, python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
        gdal=gdal.VersionInfo("RELEASE_NAME"), overlap_priority="later valid input wins")
    def save():
        (out/"run.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    save()
    try:
        manifest["stage"] = "input_conversion"
        save()
        c, input_report = prepare_inputs(c, out/"input", gdal, osr)
        manifest["input_report"] = input_report
        manifest["sources"] = check_sources(c, gdal, osr)
        manifest["stage"] = "mosaic"
        save()
        print(f"Mosaic: {len(c['inputs'])} DEM files", flush=True)
        print("Slope algorithm: " + c.get("slope_algorithm", "horn"), flush=True)
        print(smoothing_summary(c), flush=True)
        kwargs = dict(resolution="highest", VRTNodata=NODATA, strict=True)
        if c["source_nodata"] is not None:
            kwargs["srcNodata"] = c["source_nodata"]
        vrt = gdal.BuildVRT(str(out/"source_mosaic.vrt"), c["inputs"], **kwargs)
        if vrt is None:
            raise RuntimeError("BuildVRT failed")
        vrt.FlushCache()
        projected_path=out/'projected_dem.vrt'
        projected = gdal.Warp(str(projected_path), vrt, format="VRT", dstSRS=c["target_crs"],
            xRes=c["cell_size"], yRes=c["cell_size"], targetAlignedPixels=True,
            resampleAlg="bilinear", srcNodata=NODATA, dstNodata=NODATA,
            outputType=gdal.GDT_Float32, errorThreshold=0)
        gt = projected.GetGeoTransform()
        bounds = (gt[0], gt[3]+gt[5]*projected.RasterYSize,
                  gt[0]+gt[1]*projected.RasterXSize, gt[3])
        # Validate domain and sheet count before materializing a large mosaic.
        manifest["candidate_sheets"] = sum(1 for _ in intersecting_sheets(bounds,
            c["plane_zone"], c["sheet_level"], c["max_sheets"]))
        print(f"Projected mosaic: {projected.RasterXSize} x {projected.RasterYSize}", flush=True)
        projected.FlushCache();projected=vrt=None
        manifest['stage']='sheet_streaming';save()
        cs_mosaic,manifest['sheets'],manifest['valid_relief_pixels']=render_sheets(projected_path,out,c,gdal,ogr,osr)
        if c.get("xyz_enabled", True):
            manifest["stage"] = "xyz"
            save()
            manifest["xyz"] = write_xyz(cs_mosaic, out/"xyz", c, gdal, osr)
        manifest["stage"] = "finished"
        manifest["status"] = "completed"
        manifest["sheet_template"] = "sheets/{sheet_code}.tif"
        manifest["cs_mosaic"] = "cs_relief.vrt"
        save()
        print(f"COMPLETED: {out}", flush=True)
    except BaseException as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        save()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument('--input-type', choices=['raster','gsi','text','lidar','forest'])
    source_group=parser.add_mutually_exclusive_group()
    source_group.add_argument('--input-folder', type=Path)
    source_group.add_argument('--input-file', type=Path, action='append')
    parser.add_argument('--input-crs')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--lidar-mode', choices=['classified','ground_only','auto'])
    parser.add_argument('--pdal', type=Path)
    args = parser.parse_args()
    overrides = {}
    for attr,key in [('input_type','input_type'),('input_crs','input_crs'),('lidar_mode','lidar_mode')]:
        value=getattr(args,attr)
        if value is not None:overrides[key]=value
    if args.input_folder is not None:overrides['inputs']=[str(args.input_folder.resolve())]
    if args.input_file is not None:overrides['inputs']=[str(path.resolve()) for path in args.input_file]
    if args.output_dir is not None:overrides['output_dir']=str(args.output_dir.resolve())
    if args.pdal is not None:overrides['pdal_path']=str(args.pdal.resolve())
    run(read_config(args.config.resolve(),overrides))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
