import os, sys
from datetime import datetime
from ij import IJ, Prefs
from ij.io import FileSaver

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "marcksl1 shape prescreener";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

import tubefitter as tbf
import projection as prj

Prefs.blackBackground = True;

# bulbous data
im_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Ml1bEGFP\\Cropped\\2019-03-01 10-18-42 output\\Cropped dextran-rhodamine ml1b-egfp 2019-02-27 cell shape analysis ee2a xISV 1.tif";
metadata_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Ml1bEGFP\\Cropped\\2019-03-01 10-18-42 output\\Cropped dextran-rhodamine ml1b-egfp 2019-02-27 cell shape analysis ee2a xISV 1.json";
output_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Ml1bEGFP\\Straightening output";

# original data
im_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.tif";
metadata_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.json";
output_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\straightening output";


print("starting tube fitting for vessel axis extraction...");
t1 = datetime.now();
tt1 = t1;
imp = tbf.do_tubefitting(im_path, metadata_path=metadata_path, output_path=output_path, save_output=True);
t2 = datetime.now();
print("done tubefitting in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
print(imp);
Prefs.blackBackground = True;
print("starting angular projection of vessel surface...");
t1 = datetime.now();
cal = imp.getCalibration();
print("pre-projection pixel wxhxd = {}x{}x{}".format(cal.pixelWidth, cal.pixelHeight, cal.pixelDepth));
out_imp, _, ring_rois, centres = prj.do_angular_projection(imp, generate_roi_stack=True, min_r_pix=5, max_r_pix=30);
t2 = datetime.now();
print("done projection in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
print("starting unwrapping of cell from around vessel wall...");
t1 = datetime.now();
IJ.setAutoThreshold(out_imp, "Intermodes dark");
threshold_val = out_imp.getProcessor().getMinThreshold();
unwrapped_projection_imp, unwrap_axis = prj.twist_and_unwrap(out_imp);
FileSaver(unwrapped_projection_imp).saveAsTiff(os.path.join(output_path, "projection.tif"));
t2 = datetime.now();
print("done unwrapping in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
print("start generation of image colormapped by radius...");
t1 = datetime.now();
r_imp, mask_imp = prj.generate_r_image(out_imp, ring_rois, centres, unwrap_axis, threshold_val);
t2 = datetime.now();
print("done generating r image in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
FileSaver(r_imp).saveAsTiff(os.path.join(output_path, "radius image.tif"));
r_imp.show();
out_imp.close();
unwrapped_projection_imp.show();
print("start calculation of cell inner surface area and aspect ratio...");
t1 = datetime.now();
area, aspect_ratio = prj.calculate_area_and_aspect_ratio(r_imp, mask_imp, 0.108333);
t2 = datetime.now();
print("done calculating outputs in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
print("A = " + str(area));
print("Aspect ratio = " + str(aspect_ratio));
t2 = datetime.now();
print("Finished in {} s. ".format(round(0.001*(t2.getTime() - tt1.getTime()))));