import unittest
import numpy as np
from csmap_pipeline import relief, color_settings


class ColorTests(unittest.TestCase):
    def render(self, tone, a=None):
        if a is None:
            y,x=np.mgrid[:81,:81]
            a=500+5*np.sin(x/6)+3*np.cos(y/9)
        return relief(a,np.isfinite(a),1,2,.05,60,[0,3000],tone)

    def test_only_colours_change(self):
        base,slope,curvature=self.render({})
        edited,s2,c2=self.render(dict(saturation=.3,gamma=1.2,slope_darkness=.2))
        np.testing.assert_array_equal(slope,s2)
        np.testing.assert_array_equal(curvature,c2)
        np.testing.assert_array_equal(base[:,:,3],edited[:,:,3])
        self.assertFalse(np.array_equal(base[:,:,:3],edited[:,:,:3]))

    def test_grayscale_and_transparency(self):
        a=np.full((81,81),100.)
        a[40,40]=np.nan
        img,_,_=self.render(dict(saturation=0,brightness=2,gamma=1.5),a)
        np.testing.assert_array_equal(img[:,:,0],img[:,:,1])
        np.testing.assert_array_equal(img[:,:,1],img[:,:,2])
        self.assertFalse(img[40,40].any())

    def test_custom_palette(self):
        img,_,_=self.render(dict(neutral_rgb=[40,80,120],curvature_strength=0,
                                elevation_mix=0,slope_darkness=0))
        np.testing.assert_array_equal(img[40,40], [40,80,120,255])

    def test_middle_colour_is_an_independent_stop(self):
        tone=dict(valley_rgb=[0,0,255],neutral_rgb=[255,255,0],ridge_rgb=[255,0,0],
                  curvature_strength=1,elevation_mix=0,slope_darkness=0)
        flat=np.full((81,81),500.)
        img,_,curvature=self.render(tone,flat)
        self.assertAlmostEqual(float(curvature[40,40]),0.0)
        np.testing.assert_array_equal(img[40,40], [255,255,0,255])

    def test_fme_middle_colour_changes_flat_curvature_output(self):
        flat=np.full((81,81),500.)
        base=relief(flat,np.isfinite(flat),1,2,.03,60,[200,2000],{},'fme')[0]
        edited=relief(flat,np.isfinite(flat),1,2,.03,60,[200,2000],
                      {'neutral_rgb':[0,255,0]},'fme')[0]
        self.assertFalse(np.array_equal(base[40,40,:3],edited[40,40,:3]))
        np.testing.assert_array_equal(base[:,:,3],edited[:,:,3])

    def test_brightness_direction(self):
        base,_,_=self.render({})
        dark,_,_=self.render(dict(brightness=.7))
        light,_,_=self.render(dict(gamma=1.4))
        self.assertTrue(np.all(dark[40,40,:3]<base[40,40,:3]))
        self.assertTrue(np.all(light[40,40,:3]>base[40,40,:3]))

    def test_validation(self):
        for invalid in ({'valley_rgb':[256,0,0]}, {'ridge_rgb':[0,1]},
                        {'brightness':float('nan')}, {'gamma':0},
                        {'slope_darkness':1.1}, {'saturation':-1}, {'typo':1}):
            with self.assertRaises(ValueError):color_settings(invalid)


if __name__=='__main__':unittest.main(verbosity=2)
