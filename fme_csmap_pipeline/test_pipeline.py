"""Numerical tests; GDAL integration test runs only when GDAL is installed."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from csmap_pipeline import (relief, run, read_config, validate_stretch, slope_gradients,
                            rendering_settings,
                            FME_MANUAL_CURVATURE_LIMIT, FORESTRY_TUNED_CURVATURE_LIMIT)


class ReliefTests(unittest.TestCase):
    def render(self, a, cell=1, sigma=2):
        return relief(a, np.isfinite(a), cell, sigma, .05, 60, [0, 3000])

    def test_flat(self):
        a = np.full((81, 81), 100.)
        rgba, slope, curv = self.render(a)
        self.assertEqual(int(rgba[40, 40, 3]), 255)
        self.assertAlmostEqual(slope[40, 40], 0)
        self.assertAlmostEqual(curv[40, 40], 0)

    def test_plane_and_units(self):
        y, x = np.mgrid[:81, :81]
        for cell in (1, 5):
            _, slope, curv = self.render(100 + x*cell*.2 + y*cell*.3, cell=cell)
            self.assertAlmostEqual(slope[40, 40], np.degrees(np.arctan(np.hypot(.2, .3))))
            self.assertAlmostEqual(curv[40, 40], 0)

    def test_bowl_hill_sign_and_colour(self):
        y, x = np.mgrid[-40:41, -40:41]
        bowl = 100 + .01*(x*x+y*y)
        hill = 100 - .01*(x*x+y*y)
        for a, sign in ((bowl, -1), (hill, 1)):
            rgba, _, curv = self.render(a)
            self.assertAlmostEqual(curv[40, 40], sign*.04)
            if sign == 1:
                self.assertGreater(int(rgba[40, 40, 0]), int(rgba[40, 40, 2]))
            else:
                self.assertLess(int(rgba[40, 40, 0]), int(rgba[40, 40, 2]))

    def test_nodata_and_outer_edge(self):
        a = np.full((81, 81), 100.)
        a[40, 40] = np.nan
        rgba, _, _ = self.render(a)
        self.assertEqual(int(rgba[40, 49, 3]), 0)  # radius 8 + derivative 1
        self.assertEqual(int(rgba[40, 50, 3]), 255)
        self.assertEqual(int(rgba[0, 40, 3]), 0)

    def test_all_nodata(self):
        rgba, _, _ = self.render(np.full((81, 81), np.nan))
        self.assertFalse(np.any(rgba))

    def test_block_seams(self):
        y, x = np.mgrid[:163, :175]
        a = 100 + np.sin(x/7)*3 + np.cos(y/9)*5
        a[81, 89] = np.nan
        expected, _, _ = self.render(a)
        actual = np.zeros_like(expected)
        halo, block = 9, 37
        for yy in range(0, a.shape[0], block):
            for xx in range(0, a.shape[1], block):
                x0, y0 = max(0, xx-halo), max(0, yy-halo)
                x1, y1 = min(a.shape[1], xx+block+halo), min(a.shape[0], yy+block+halo)
                rgba, _, _ = self.render(a[y0:y1, x0:x1])
                bw, bh = min(block, a.shape[1]-xx), min(block, a.shape[0]-yy)
                actual[yy:yy+bh, xx:xx+bw] = rgba[yy-y0:yy-y0+bh, xx-x0:xx-x0+bw]
        np.testing.assert_array_equal(actual, expected)


class ConfigDefaultsTests(unittest.TestCase):
    """Regression tests for the curvature_limit profile-mixing bug (VALIDATION.txt 2026-09-28).

    color_model=legacy (independent_v040相当) and color_model=fme must never share
    a curvature_limit default: the FME manual's confirmed value (0.1) and the
    provisional value (0.03) are distinct profiles that must not silently
    leak into each other via an omitted curvature_limit key.
    """

    def _write(self, tmp, extra):
        root = Path(tmp)
        (root / "dummy.tif").touch()
        base = dict(inputs=[str(root / "dummy.tif")], output_dir=str(root / "out"),
                    target_crs="EPSG:6677", cell_size=1.0, confirm_elevation_metres=True)
        base.update(extra)
        config = root / "config.json"
        config.write_text(json.dumps(base), encoding="utf-8")
        return config

    def test_bare_defaults_match_independent_v040(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = read_config(self._write(tmp, {}))
            self.assertEqual(c["color_model"], "legacy")
            self.assertEqual(c["curvature_limit"], 0.05)
            self.assertEqual(c["elevation_range"], [0, 3000])

    def test_fme_manual_bare_defaults_unaffected(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = read_config(self._write(tmp, dict(color_model="fme")))
            self.assertEqual(c["curvature_limit"], FME_MANUAL_CURVATURE_LIMIT)

    def test_fme_explicit_forestry_tuned_value_not_overridden(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = read_config(self._write(tmp, dict(color_model="fme",
                            curvature_limit=FORESTRY_TUNED_CURVATURE_LIMIT)))
            self.assertEqual(c["curvature_limit"], FORESTRY_TUNED_CURVATURE_LIMIT)

    def test_curvature_profiles_are_distinct(self):
        self.assertNotEqual(FME_MANUAL_CURVATURE_LIMIT, FORESTRY_TUNED_CURVATURE_LIMIT)
        self.assertAlmostEqual(FME_MANUAL_CURVATURE_LIMIT, 0.1)
        self.assertAlmostEqual(FORESTRY_TUNED_CURVATURE_LIMIT, 0.03)

    def test_slope_algorithm_defaults_to_horn_and_is_validated(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(read_config(self._write(tmp, {}))['slope_algorithm'], 'horn')
        with tempfile.TemporaryDirectory() as tmp:
            c = read_config(self._write(tmp, dict(slope_algorithm='central_difference')))
            self.assertEqual(c['slope_algorithm'], 'central_difference')
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, 'slope_algorithm must be one of'):
                read_config(self._write(tmp, dict(slope_algorithm='bogus')))


class SlopeAlgorithmTests(unittest.TestCase):
    """v0.8.0: 傾斜計算方式の選択（ユーザー決定、2026-10-02）。csmap-sheets v0.10.0と同じ式。"""

    def test_constant_surface_gives_zero_slope_for_both_algorithms(self):
        z = np.full((9, 9), 123.4)
        for algorithm in ('horn', 'central_difference'):
            dx, dy = slope_gradients(z, 2.0, algorithm)
            np.testing.assert_allclose(np.degrees(np.arctan(np.hypot(dx, dy))), 0.0, atol=1e-12)

    def test_planar_slope_matches_analytic_angle(self):
        yy, xx = np.mgrid[:11, :11].astype(float)
        cell, gx, gy = 0.5, 0.7, -0.4
        plane = 50 + gx * xx * cell + gy * yy * cell
        expected = np.degrees(np.arctan(np.hypot(gx, gy)))
        for algorithm in ('horn', 'central_difference'):
            dx, dy = slope_gradients(plane, cell, algorithm)
            self.assertAlmostEqual(float(np.degrees(np.arctan(np.hypot(dx, dy)))[5, 5]), expected, places=9)

    def test_horn_matches_hand_calculation_on_nonlinear_3x3(self):
        z = np.array([[1, 2, 4], [3, 5, 8], [6, 9, 13]], dtype=float)
        dx, dy = slope_gradients(z, 1.0, 'horn')
        self.assertAlmostEqual(abs(float(dx[1, 1])), 2.5)
        self.assertAlmostEqual(abs(float(dy[1, 1])), 3.5)

    def test_central_difference_reproduces_v0_7_2_formula(self):
        rng = np.random.default_rng(7)
        z = rng.normal(300, 20, (17, 19)); cell = 1.5
        dx, dy = slope_gradients(z, cell, 'central_difference')
        np.testing.assert_array_equal(dx, (np.roll(z, -1, 1) - np.roll(z, 1, 1)) / (2 * cell))
        np.testing.assert_array_equal(dy, (np.roll(z, -1, 0) - np.roll(z, 1, 0)) / (2 * cell))

    def test_relief_default_is_horn_and_invalid_raises(self):
        y, x = np.mgrid[:41, :41]
        a = 500 + 5 * np.sin(x / 6) + 3 * np.cos(y / 9)
        default = relief(a, np.isfinite(a), 1, 2, .05, 60, [0, 3000])
        horn = relief(a, np.isfinite(a), 1, 2, .05, 60, [0, 3000], slope_algorithm='horn')
        np.testing.assert_array_equal(default[1], horn[1])
        with self.assertRaisesRegex(ValueError, 'slope_algorithm must be one of'):
            relief(a, np.isfinite(a), 1, 2, .05, 60, [0, 3000], slope_algorithm='bogus')

    def test_horn_nodata_neighbourhood_and_edges_are_transparent(self):
        y, x = np.mgrid[:21, :21]
        a = (200 + 3 * np.sin(x / 3) + 2 * np.cos(y / 4)).astype(float)
        valid = np.ones(a.shape, bool); valid[10, 10] = False; a[10, 10] = np.nan
        for model in ('legacy', 'fme'):
            alpha = relief(a, valid, 1, 0, .05, 60, [0, 3000], {}, model, slope_algorithm='horn')[0][:, :, 3]
            for r, c in ((9, 9), (9, 11), (11, 9), (11, 11), (10, 9), (9, 10)):
                self.assertEqual(alpha[r, c], 0, (model, r, c))
            self.assertGreater(alpha[12, 12], 0)
            self.assertTrue((alpha[0, :] == 0).all() and (alpha[:, -1] == 0).all())

    def test_block_split_has_no_seam_for_both_algorithms(self):
        y, x = np.mgrid[:60, :50]
        a = 500 + 5 * np.sin(x / 6) + 3 * np.cos(y / 9) + .002 * (x - 25) ** 2
        valid = np.isfinite(a); sigma = 2; halo = int(np.ceil(4 * sigma / 1)) + 1
        for algorithm in ('horn', 'central_difference'):
            full = relief(a, valid, 1, sigma, .05, 60, [0, 3000], slope_algorithm=algorithm)[0]
            part = relief(a[:40], valid[:40], 1, sigma, .05, 60, [0, 3000], slope_algorithm=algorithm)[0]
            np.testing.assert_array_equal(part[halo:40 - halo], full[halo:40 - halo])

    def test_rendering_record_reflects_slope_algorithm_for_both_models(self):
        for model in ('legacy', 'fme'):
            for algorithm in ('horn', 'central_difference'):
                c = dict(color_model=model, elevation_range=[0, 3000], slope_max=60,
                         curvature_limit=.1, sigma_m=3, color={}, slope_algorithm=algorithm)
                record = rendering_settings(c)['terrain_calculation']
                self.assertEqual(record['slope_algorithm'], algorithm)


class StretchModeTests(unittest.TestCase):
    def test_none_returns_no_stretch(self):
        self.assertIsNone(validate_stretch("none", None))

    def test_nagano_reference_returns_manual_values(self):
        self.assertEqual(validate_stretch("nagano_reference", None),
                         ((65, 234), (57, 229), (66, 216)))

    def test_custom_range_is_validated_and_resolved(self):
        resolved = validate_stretch("custom", [[10, 200], [0, 255], [20, 220]])
        self.assertEqual(resolved, ((10.0, 200.0), (0.0, 255.0), (20.0, 220.0)))

    def test_invalid_stretch_mode_rejected(self):
        with self.assertRaises(ValueError):
            validate_stretch("bogus", None)

    def test_custom_range_missing_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_stretch("custom", None)

    def test_custom_range_out_of_bounds_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_stretch("custom", [[10, 300], [0, 255], [20, 220]])

    def test_relief_stretch_none_vs_nagano_reference_differ(self):
        flat = np.full((81, 81), 500.)
        none_stretch = relief(flat, np.isfinite(flat), 1, 2, .1, 60, [200, 2000],
                              {}, 'fme', stretch_mode='none')[0]
        nagano = relief(flat, np.isfinite(flat), 1, 2, .1, 60, [200, 2000],
                        {}, 'fme', stretch_mode='nagano_reference')[0]
        self.assertFalse(np.array_equal(none_stretch[:, :, :3], nagano[:, :, :3]))


@unittest.skipUnless(importlib.util.find_spec("osgeo"), "GDAL is not installed")
class IntegrationTests(unittest.TestCase):
    def test_two_dem_mosaic_to_map_sheets(self):
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
                ds.GetRasterBand(1).WriteArray((100+(xx+i*96)*.2+yy*.1).astype("float32"))
                ds.GetRasterBand(1).SetNoDataValue(-9999)
                ds = None
            config = root/"config.json"
            config.write_text(json.dumps(dict(inputs=[str(root/"dem_*.tif")],
                output_dir=str(root/"result"), target_crs="EPSG:6677", cell_size=1,
                sigma_m=2, plane_zone=9, sheet_level=500,
                confirm_elevation_metres=True)), encoding="utf-8")
            run(read_config(config))
            self.assertTrue((root/'result/projected_dem.vrt').exists())
            rgba = gdal.Open(str(root/"result/cs_relief.vrt"))
            self.assertEqual((rgba.RasterXSize,rgba.RasterYSize),(800,300))
            self.assertEqual(int(rgba.GetRasterBand(4).ReadAsArray()[40,400]),255)
            expected=rgba.ReadAsArray()[:,:96,304:496]
            rgba=None
            files=sorted((root/"result/sheets").glob("*.tif"))
            self.assertEqual([f.stem for f in files],["09KD0909","09KE0000"])
            left=gdal.Open(str(files[0]))
            right=gdal.Open(str(files[1]))
            self.assertEqual((left.RasterXSize,left.RasterYSize),(400,300))
            self.assertEqual(left.GetGeoTransform(),(-400.,1.,0.,0.,0.,-1.))
            self.assertEqual(right.GetGeoTransform(),(0.,1.,0.,0.,0.,-1.))
            reconstructed=np.concatenate((left.ReadAsArray()[:,:96,304:],right.ReadAsArray()[:,:96,:96]),axis=2)
            np.testing.assert_array_equal(reconstructed,expected)
            self.assertFalse(np.any(left.GetRasterBand(4).ReadAsArray()[:, :304]))
            self.assertEqual(files[0].with_suffix('.tfw').read_text().splitlines()[-2:],
                             ['-399.500000000000','-0.500000000000'])
            left=right=None
            self.assertTrue((root/'result/sheet_index.gpkg').exists())
            self.assertEqual(json.loads((root/"result/run.json").read_text(encoding="utf-8"))["status"], "completed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
