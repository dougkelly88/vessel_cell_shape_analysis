# unit tests for cell shape analysis code
#
# D. J. Kelly, 2018-12-05, douglas.kelly@riken.jp

import os, sys, unittest, math

from ij import IJ
from ij.gui import EllipseRoi

release = True;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "Unwrapping and shape analysis";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import ellipse_fitting

class TestEllipseFit(unittest.TestCase):

	def do_ellipse_testing(self, imw, imh, centre, angle, axl):
		"""Utility function to provide test ellipse data"""
		imp = IJ.createImage("test", imw, imh, 1, 8);
		imp.show()
		roi = ellipse_fitting.generate_ellipse_roi(centre, angle, axl);
		imp.setRoi(roi);
		IJ.run(imp, "Set...", "value=128");
		IJ.setThreshold(imp, 1, 255);
		IJ.run(imp, "Convert to Mask", "method=Default background=Dark list");
		IJ.run(imp, "Outline", "stack");
		IJ.run(imp, "Create Selection", "");
		roi = imp.getRoi();
		pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()];
		imp.changes = False;
		imp.close();
		return pts;

	def test_fit_ellipse_center_of_image(self):
		"""Verify that outputs match test input parameters for non-occluded ellipse, with allowance for pixelation"""
		c_in = (250, 250);
		ang_in = 0;
		axl_in = (150,50);
		pts = self.do_ellipse_testing(500, 500, c_in, ang_in, axl_in);
		centre, angle, axl = ellipse_fitting.fit_ellipse(pts);
		print("centre = " + str(centre));
		print("angle (deg) = " + str(angle * 180 / math.pi));
		print("axes = " + str(axl));
		#self.assertEqual(centre, c);
		#self.assertEqual(angle, ang_in);
		#self.assertEqual(axl, axl_in);

	def test_fit_ellipse_edge_of_image(self):
		"""Verify that outputs match test input parameters for partially occluded ellipse, with allowance for pixelation"""
		c_in = (250, 475);
		ang_in = math.pi / 2;
		axl_in = (150, 50);
		pts = self.do_ellipse_testing(500, 500, c_in, ang_in, axl_in);
		centre, angle, axl = ellipse_fitting.fit_ellipse(pts);
		print("centre = " + str(centre));
		print("angle (deg) = " + str(angle * 180 / math.pi));
		print("axes = " + str(axl));
		#self.assertEqual(centre, c);
		#self.assertEqual(angle, ang_in);
		#self.assertEqual(axl, axl_in);

	def test_fit_ellipse_oblique(self):
		"""Verify that outputs match test input parameters for vertically oriented ellipse, with allowance for pixelation"""
		c_in = (250, 250);
		ang_in = math.pi / 2;
		axl_in = (150,50);
		pts = self.do_ellipse_testing(500, 500, c_in, ang_in, axl_in);
		centre, angle, axl = ellipse_fitting.fit_ellipse(pts);
		print("centre = " + str(centre));
		print("angle (deg) = " + str(angle * 180 / math.pi));
		print("axes = " + str(axl));
		#self.assertEqual(centre, c);
		#self.assertEqual(angle, ang_in);
		#self.assertEqual(axl, axl_in);

if __name__ in ['__builtin__','__main__']:
	suite = unittest.TestLoader().loadTestsFromTestCase(TestEllipseFit)
	unittest.TextTestRunner(verbosity=2).run(suite)
