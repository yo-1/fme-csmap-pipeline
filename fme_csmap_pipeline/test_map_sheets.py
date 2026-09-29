import unittest
from map_sheets import dimensions, sheet_at, intersecting_sheets, pixel_window


class MapSheetTests(unittest.TestCase):
    def test_dimensions(self):
        self.assertEqual([dimensions(x) for x in (5000,2500,1000,500)],
                         [(4000,3000),(2000,1500),(800,600),(400,300)])

    def test_reference_09ld00(self):
        # GSI Appendix 7 Art.84: L is [-60000,-30000] northing;
        # D is [-40000,0] easting; 00 is its NW 5000-level cell.
        s = sheet_at(-39999,-30001,9,5000)
        self.assertEqual(s.code,"09LD00")
        self.assertEqual((s.west,s.south,s.east,s.north),(-40000,-33000,-36000,-30000))

    def test_2500_quadrants(self):
        points = [(-39000,-30500),(-37000,-30500),(-39000,-32500),(-37000,-32500)]
        self.assertEqual([sheet_at(e,n,9,2500).code for e,n in points],
                         ["09LD001","09LD002","09LD003","09LD004"])

    def test_1000_suffix(self):
        self.assertEqual(sheet_at(-39999,-30001,9,1000).code,"09LD000A")
        self.assertEqual(sheet_at(-36001,-32999,9,1000).code,"09LD004E")

    def test_500_suffix(self):
        self.assertEqual(sheet_at(-39999,-30001,9,500).code,"09LD0000")
        self.assertEqual(sheet_at(-36001,-32999,9,500).code,"09LD0099")

    def test_origin_and_boundary_ownership(self):
        self.assertEqual(sheet_at(0,0,9,5000).code,"09KE00")
        self.assertEqual(sheet_at(-.01,0,9,5000).code,"09KD09")
        self.assertEqual(sheet_at(0,.01,9,5000).code,"09JE90")
        sheets=list(intersecting_sheets((-40000,-33000,-36000,-30000),9,5000))
        self.assertEqual([s.code for s in sheets],["09LD00"])

    def test_adjacency_and_coverage(self):
        for level in (5000,2500,1000,500):
            w,h=dimensions(level)
            sheets=list(intersecting_sheets((-w,-h,w,h),9,level))
            self.assertEqual(len(sheets),4)
            self.assertEqual(sheets[0].east,sheets[1].west)
            self.assertEqual(sheets[0].south,sheets[2].north)
            self.assertEqual(len({s.code for s in sheets}),4)
            self.assertEqual(sum((s.east-s.west)*(s.north-s.south) for s in sheets),4*w*h)

    def test_limits(self):
        self.assertEqual(sheet_at(-160000,300000,1,5000).code,"01AA00")
        self.assertEqual(sheet_at(159999,-299999,19,5000).code,"19TH99")
        for e,n in ((160000,0),(-160001,0),(0,-300000),(0,300001)):
            with self.assertRaises(ValueError):sheet_at(e,n,9,5000)
        with self.assertRaises(ValueError): list(intersecting_sheets((-4000,-3000,4000,3000),9,500,1))

    def test_padding_window_and_axis_order(self):
        s=sheet_at(-1,-1,9,500)
        p=pixel_window(s,(-96,1,0,0,0,-1),192,96)
        self.assertEqual((p['src_x'],p['src_y'],p['width'],p['height']),(0,0,96,96))
        self.assertEqual((p['dst_x'],p['dst_y'],p['tile_width'],p['tile_height']),(304,0,400,300))
        with self.assertRaises(ValueError):pixel_window(s,(-95.5,1,0,0,0,-1),192,96)


if __name__=="__main__":unittest.main(verbosity=2)
