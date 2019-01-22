import os, sys
from ij import IJ, Prefs

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "marcksl1 shape prescreener";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, script_path);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import tubefitter as tbf
import projection as prj

Prefs.blackBackground = True;

im_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.tif";
metadata_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.json";
output_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\straightening output";

imp = tbf.do_tubefitting(im_path, metadata_path, output_path);
print("done tubefitting");
print(imp);
Prefs.blackBackground = True;
out_imp, _, ring_rois, centres = prj.do_angular_projection(imp, generate_roi_stack=True);
print("done projection")
IJ.setAutoThreshold(out_imp, "Intermodes dark");
threshold_val = out_imp.getProcessor().getMinThreshold();
unwrapped_projection_imp, unwrap_axis = prj.twist_and_unwrap(out_imp);
print("done unwrapping")
r_imp, mask_imp = prj.generate_r_image(out_imp, ring_rois, centres, unwrap_axis, threshold_val);
print("done generating r image")
r_imp.show();
out_imp.close();
unwrapped_projection_imp.show();
area, aspect_ratio = prj.calculate_area_and_aspect_ratio(r_imp, mask_imp, 0.108333);
print("done calculating outputs")
print("A = " + str(area));
print("Aspect ratio = " + str(aspect_ratio));