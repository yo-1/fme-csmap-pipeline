"""nodata_edge_mode is reserved for csmap-sheets compatibility; only safe_mask is implemented."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from csmap_pipeline import read_config, rendering_settings, run


def write_config(root, extra):
    (root / "dummy.tif").touch()
    base = dict(inputs=[str(root / "dummy.tif")], output_dir=str(root / "out"),
                target_crs="EPSG:6677", cell_size=1.0, confirm_elevation_metres=True,
                slope_algorithm="horn")
    base.update(extra)
    config = root / "config.json"
    config.write_text(json.dumps(base), encoding="utf-8")
    return config


class NodataEdgeModeConfigTests(unittest.TestCase):
    def read(self, extra):
        with tempfile.TemporaryDirectory() as tmp:
            return read_config(write_config(Path(tmp), extra))

    def test_missing_key_is_safe_mask(self):
        self.assertEqual(self.read({})["nodata_edge_mode"], "safe_mask")

    def test_explicit_safe_mask(self):
        self.assertEqual(self.read(dict(nodata_edge_mode="safe_mask"))["nodata_edge_mode"],
                         "safe_mask")

    def test_pss_approximation_is_rejected_as_not_implemented(self):
        with self.assertRaisesRegex(ValueError, "not implemented"):
            self.read(dict(nodata_edge_mode="pss_approximation"))

    def test_unknown_values_are_rejected(self):
        for value in ("bogus", "", "SAFE_MASK", "safe_mask "):
            with self.subTest(value=value), \
                    self.assertRaisesRegex(ValueError, "nodata_edge_mode must be one of"):
                self.read(dict(nodata_edge_mode=value))

    def test_non_string_values_are_rejected(self):
        for value in (None, 0, 1, True, ["safe_mask"], {"mode": "safe_mask"}):
            with self.subTest(value=value), \
                    self.assertRaisesRegex(ValueError, "must be a string"):
                self.read(dict(nodata_edge_mode=value))

    def test_rendering_record_for_both_color_models(self):
        for model in ("legacy", "fme"):
            c = dict(color_model=model, elevation_range=[0, 3000], slope_max=60,
                     curvature_limit=.1, sigma_m=3, color={}, slope_algorithm="horn")
            for settings in (c, {**c, "nodata_edge_mode": "safe_mask"}):
                record = rendering_settings(settings)["terrain_calculation"]
                self.assertEqual(record["nodata_edge_mode"], "safe_mask")
                self.assertIn("NoData", record["nodata_edge_handling"])

    def test_bundled_configs_do_not_set_the_key(self):
        configs = sorted(Path(__file__).parent.glob("config*.json"))
        self.assertEqual(len(configs), 14)
        for path in configs:
            with self.subTest(config=path.name):
                data = json.loads(path.read_text(encoding="utf-8-sig"))
                self.assertNotIn("nodata_edge_mode", data)


@unittest.skipUnless(importlib.util.find_spec("osgeo"), "GDAL is not installed")
class NodataEdgeModeIntegrationTests(unittest.TestCase):
    def test_missing_key_and_safe_mask_give_identical_sheets(self):
        from osgeo import gdal, osr
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(6677)
            ds = gdal.GetDriverByName("GTiff").Create(str(root / "dem.tif"), 96, 96, 1,
                                                      gdal.GDT_Float32)
            ds.SetProjection(srs.ExportToWkt())
            ds.SetGeoTransform((-48, 1, 0, 0, 0, -1))
            yy, xx = np.mgrid[:96, :96]
            z = (100 + xx * .2 + 5 * np.sin(yy / 7.0)).astype("float32")
            z[40:45, 50:55] = -9999
            ds.GetRasterBand(1).WriteArray(z)
            ds.GetRasterBand(1).SetNoDataValue(-9999)
            ds = None
            arrays = {}
            for name, extra in (("missing", {}), ("explicit", dict(nodata_edge_mode="safe_mask"))):
                config = dict(inputs=[str(root / "dem.tif")], output_dir=str(root / name),
                              target_crs="EPSG:6677", cell_size=1, sigma_m=2, plane_zone=9,
                              sheet_level=500, confirm_elevation_metres=True,
                              slope_algorithm="horn", **extra)
                path = root / f"{name}.json"
                path.write_text(json.dumps(config), encoding="utf-8")
                run(read_config(path))
                manifest = json.loads((root / name / "run.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["status"], "completed")
                self.assertEqual(manifest["settings"]["nodata_edge_mode"], "safe_mask")
                self.assertEqual(manifest["rendering"]["terrain_calculation"]["nodata_edge_mode"],
                                 "safe_mask")
                image = gdal.Open(str(root / name / "cs_relief.vrt"))
                arrays[name] = image.ReadAsArray()
                image = None
            np.testing.assert_array_equal(arrays["missing"], arrays["explicit"])
            # Mosaic origin is x=-400, so the DEM starts at column 352. The NoData
            # hole stays transparent while valid terrain away from it is opaque.
            alpha = arrays["missing"][3]
            self.assertFalse(np.any(alpha[40:45, 402:407]))
            self.assertEqual(int(alpha[20, 372]), 255)


if __name__ == "__main__":
    unittest.main(verbosity=2)
