"""Numerical tests; GDAL integration test runs only when GDAL is installed."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from csmap_pipeline import (relief, run, read_config, validate_stretch,
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
