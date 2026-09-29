import unittest
import tempfile
from pathlib import Path
import numpy as np
from xyz_tiles import *

class XYZTests(unittest.TestCase):
    def test_world(self):
        h=HALF_WORLD
        self.assertEqual(tile_range((-h,-h,h,h),0),(0,0,0,0))
        self.assertEqual(tile_range((-h,-h,h,h),2),(0,3,0,3))

    def test_north_origin_and_adjacent_edges(self):
        nw=tile_bounds(1,0,0);se=tile_bounds(1,1,1)
        self.assertEqual(nw,(-HALF_WORLD,0,0,HALF_WORLD))
        self.assertEqual(se,(0,-HALF_WORLD,HALF_WORLD,0))
        self.assertEqual(tile_range(nw,1),(0,0,0,0))

    def test_japan(self):
        h=HALF_WORLD
        east=139/180*h
        north=6378137*np.log(np.tan(np.pi/4+np.deg2rad(35)/2))
        self.assertEqual(tile_range((east,north,east+1,north+1),2),(3,3,1,1))

    def test_settings(self):
        validate_xyz(DEFAULTS)
        for patch in ({'xyz_min_zoom':19},{'xyz_max_zoom':25},{'xyz_max_tiles':0},
                      {'xyz_enabled':'true'},{'xyz_min_zoom':True}):
            with self.assertRaises(ValueError): validate_xyz({**DEFAULTS,**patch})

    def test_rgba_png_integration(self):
        try: from osgeo import gdal,osr
        except ImportError: self.skipTest('GDAL runtime unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); src=root/'rgba.tif'
            ds=gdal.GetDriverByName('GTiff').Create(str(src),256,256,4,gdal.GDT_Byte)
            h=HALF_WORLD; ds.SetGeoTransform((-h,h/256,0,h,0,-h/256))
            crs=osr.SpatialReference();crs.ImportFromEPSG(3857);ds.SetProjection(crs.ExportToWkt())
            for i,value in enumerate((180,90,45,255),1):
                ds.GetRasterBand(i).Fill(value)
            ds.GetRasterBand(4).SetColorInterpretation(gdal.GCI_AlphaBand)
            ds.GetRasterBand(4).WriteArray(np.zeros((128,256),dtype='uint8'),0,128)
            ds=None
            result=write_xyz(src,root/'tiles',dict(xyz_min_zoom=1,xyz_max_zoom=1),gdal,osr)
            self.assertEqual(result['status'],'completed')
            self.assertEqual(result['written_tiles'],1)
            tile=gdal.Open(str(root/'tiles/1/0/0.png'))
            self.assertEqual((tile.RasterXSize,tile.RasterYSize,tile.RasterCount),(256,256,4))
            pixels=tile.ReadAsArray()
            np.testing.assert_array_equal(pixels[:3,10,10],[180,90,45])
            self.assertEqual(pixels[3,10,10],255)
            self.assertEqual(pixels[3,200,10],0)
            tile=None
            with self.assertRaises(ValueError):
                write_xyz(src,root/'limited',dict(xyz_min_zoom=1,xyz_max_zoom=2,xyz_max_tiles=1),gdal,osr)
            self.assertFalse((root/'limited').exists())

if __name__=='__main__':unittest.main()
