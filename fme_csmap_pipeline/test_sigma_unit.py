"""v0.9.0: Gaussian σの指定方式（地上距離m／計算格子の画素数px）のテスト。

QGISプラグイン版csmap-sheets v0.12.0と同じ仕様。依頼書（2026-10-03、ユーザー承認）の
完了条件の数値と、旧設定（方式の指定なし）がm方式として従来と同じ計算条件になることを確認する。
"""
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
import numpy as np
from csmap_pipeline import (relief, run, read_config, rendering_settings, gaussian_sigma,
                            smoothing_record, smoothing_summary)


def write_config(root, **extra):
    root = Path(root)
    (root / "dummy.tif").touch()
    base = dict(inputs=[str(root / "dummy.tif")], output_dir=str(root / "out"),
                target_crs="EPSG:6677", cell_size=1.0, confirm_elevation_metres=True,
                slope_algorithm="horn")
    base.update(extra)
    path = root / "config.json"
    path.write_text(json.dumps(base), encoding="utf-8")
    return path


def config(**extra):
    with tempfile.TemporaryDirectory() as tmp:
        return read_config(write_config(tmp, **extra))


def terrain(shape=(70, 80), seed=3):
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[:shape[0], :shape[1]]
    return 300 + 8*np.sin(x/7.) + 5*np.cos(y/5.) + rng.normal(0, .3, shape)


class SigmaUnitTests(unittest.TestCase):
    def kernel(self, **kw):
        r = smoothing_record(config(**kw))
        return r['effective_sigma_m'], r['effective_sigma_px'], r['kernel_size_px']

    def test_completion_conditions(self):
        self.assertEqual(self.kernel(cell_size=.5, sigma_m=3.), (3., 6., 49))
        self.assertEqual(self.kernel(cell_size=.5, sigma_unit='px', sigma_px=3.), (1.5, 3., 25))
        self.assertEqual(self.kernel(cell_size=1., sigma_m=3.), (3., 3., 25))
        self.assertEqual(self.kernel(cell_size=1., sigma_unit='px', sigma_px=3.), (3., 3., 25))
        self.assertEqual(self.kernel(cell_size=2., sigma_m=3.), (3., 1.5, 13))
        self.assertEqual(self.kernel(cell_size=2., sigma_unit='px', sigma_px=3.), (6., 3., 25))

    def test_defaults_and_old_configs_are_metres(self):
        c = config()
        self.assertEqual((c['sigma_unit'], c['sigma_m'], c['sigma_px']), ('m', 3.0, 3.0))
        c = config(sigma_m=2.5, cell_size=.5)
        self.assertEqual(gaussian_sigma(c), (5.0, 20))
        self.assertEqual(rendering_settings(c)['smoothing']['gaussian_sigma_m'], 2.5)

    def test_metre_mode_keeps_previous_float_expressions(self):
        for cell in (.25, .5, 1., 2., 5.):
            for sigma in (0., .1, .3, .7, 1., 2.5, 3., 7.3):
                c = dict(sigma_unit='m', sigma_m=sigma, sigma_px=3., cell_size=cell)
                self.assertEqual(gaussian_sigma(c), (sigma/cell, math.ceil(4*sigma/cell)))

    def test_one_metre_grid_same_result_in_both_modes(self):
        a = terrain()
        valid = np.isfinite(a)
        m = relief(a, valid, 1., 3., .05, 60, [0, 3000])
        px = relief(a, valid, 1., 999., .05, 60, [0, 3000], sigma_px=3.)
        for left, right in zip(m, px):
            np.testing.assert_array_equal(left, right)

    def test_px_mode_matches_equivalent_metres_and_differs_from_3m(self):
        a = terrain()
        valid = np.isfinite(a)
        m = relief(a, valid, .5, 3., .05, 60, [0, 3000])
        px6 = relief(a, valid, .5, 0., .05, 60, [0, 3000], sigma_px=6.)
        np.testing.assert_array_equal(m[0], px6[0])
        px3 = relief(a, valid, .5, 3., .05, 60, [0, 3000], sigma_px=3.)
        self.assertFalse(np.allclose(m[2][30:40, 30:40], px3[2][30:40, 30:40]))

    def test_zero_sigma_disables_smoothing_in_both_modes(self):
        a = terrain()
        valid = np.isfinite(a)
        m = relief(a, valid, .5, 0., .05, 60, [0, 3000])
        px = relief(a, valid, .5, 3., .05, 60, [0, 3000], sigma_px=0.)
        np.testing.assert_array_equal(m[0], px[0])
        r = smoothing_record(config(sigma_unit='px', sigma_px=0.))
        self.assertEqual((r['effective_sigma_px'], r['kernel_size_px']), (0., 0))
        self.assertIn('none', smoothing_summary(dict(sigma_unit='px', sigma_px=0., sigma_m=3., cell_size=1.)))

    def test_fractional_px_allowed_and_invalid_values_rejected(self):
        self.assertEqual(gaussian_sigma(config(sigma_unit='px', sigma_px=1.25, cell_size=.5)), (1.25, 5))
        for bad in (-1., float('nan'), float('inf'), '3', True, None):
            with self.assertRaises(ValueError, msg=repr(bad)):
                config(sigma_unit='px', sigma_px=bad)
        for bad in (-1., float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                config(sigma_m=bad)
        with self.assertRaisesRegex(ValueError, 'sigma_unit'):
            config(sigma_unit='cm')
        with self.assertRaisesRegex(ValueError, '512'):
            config(sigma_unit='px', sigma_px=129.)

    def test_inactive_mode_value_is_kept_not_reinterpreted(self):
        c = config(sigma_unit='px', sigma_px=4., sigma_m=3., cell_size=.5)
        self.assertEqual((c['sigma_m'], c['sigma_px']), (3., 4.))
        self.assertEqual(gaussian_sigma(c), (4., 16))
        c['sigma_unit'] = 'm'
        self.assertEqual(gaussian_sigma(c), (6., 24))

    def test_px_mode_block_seams(self):
        a = terrain((90, 60))
        valid = np.isfinite(a)
        sigma_px = 2.5
        halo = math.ceil(4*sigma_px)+1
        full = relief(a, valid, .5, 0., .05, 60, [0, 3000], sigma_px=sigma_px)[0]
        part = relief(a[:50], valid[:50], .5, 0., .05, 60, [0, 3000], sigma_px=sigma_px)[0]
        np.testing.assert_array_equal(full[halo:50-halo], part[halo:50-halo])

    def test_record_keeps_setting_and_effective_values(self):
        for model in ('legacy', 'fme'):
            c = config(sigma_unit='px', sigma_px=3., sigma_m=3., cell_size=2., color_model=model)
            self.assertEqual(rendering_settings(c)['smoothing'], dict(
                gaussian_sigma_m=6., sigma_unit='px', sigma_setting=3., cell_size_m=2.,
                effective_sigma_m=6., effective_sigma_px=3., kernel_radius_px=12, kernel_size_px=25))
        self.assertEqual(smoothing_summary(c),
                         'Smoothing: sigma 3 px on a 2 m grid = 6 m = 3 px; kernel 25x25 px')


@unittest.skipUnless(importlib.util.find_spec("osgeo"), "GDAL is not installed")
class SigmaUnitIntegrationTests(unittest.TestCase):
    def run_case(self, root, name, **kw):
        path = root/f"{name}.json"
        path.write_text(json.dumps(dict(inputs=[str(root/"dem_*.tif")],
            output_dir=str(root/name), target_crs="EPSG:6677", cell_size=.5,
            plane_zone=9, sheet_level=500, confirm_elevation_metres=True,
            slope_algorithm="horn", **kw)), encoding="utf-8")
        run(read_config(path))
        return json.loads((root/name/"run.json").read_text(encoding="utf-8"))

    def test_px_mode_sheets_match_equivalent_metre_mode(self):
        from osgeo import gdal, osr
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(6677)
            for i in range(2):
                ds = gdal.GetDriverByName("GTiff").Create(str(root/f"dem_{i}.tif"), 96, 96, 1, gdal.GDT_Float32)
                ds.SetProjection(srs.ExportToWkt())
                ds.SetGeoTransform((i*96-96, 1, 0, 0, 0, -1))
                yy, xx = np.mgrid[:96, :96]
                z = 100+(xx+i*96)*.2+yy*.1+3*np.sin((xx+i*96)/6.)*np.cos(yy/5.)
                ds.GetRasterBand(1).WriteArray(z.astype("float32"))
                ds.GetRasterBand(1).SetNoDataValue(-9999)
                ds = None
            metre = self.run_case(root, "metre", sigma_m=1.)
            pixel = self.run_case(root, "pixel", sigma_unit="px", sigma_px=2., sigma_m=3.)
            self.assertEqual((metre["status"], pixel["status"]), ("completed", "completed"))
            self.assertEqual(pixel["rendering"]["smoothing"]["kernel_size_px"], 17)
            self.assertEqual(pixel["rendering"]["smoothing"]["gaussian_sigma_m"], 1.)
            names = sorted(f.name for f in (root/"metre/sheets").glob("*.tif"))
            self.assertEqual(names, sorted(f.name for f in (root/"pixel/sheets").glob("*.tif")))
            self.assertEqual(len(names), 2)
            for name in names:
                np.testing.assert_array_equal(gdal.Open(str(root/"metre/sheets"/name)).ReadAsArray(),
                                              gdal.Open(str(root/"pixel/sheets"/name)).ReadAsArray())


if __name__ == '__main__':
    unittest.main()
