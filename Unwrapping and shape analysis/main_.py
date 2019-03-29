import os, sys
from datetime import datetime
from ij import IJ, Prefs
from ij.io import FileSaver

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

import tubefitter as tbf
import projection as prj
import file_io as io

Prefs.blackBackground = True;

#root_input_path = "D:\\source\\vascular_morphogenesis_ij\\marcksl1 cell shape analysis\\v1.x.x\\test_data";
#root_output_path = "D:\\source\\vascular_morphogenesis_ij\\marcksl1 cell shape analysis\\v1.x.x\\test_data_out";

#root_input_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP\\Cropped\\2019-03-01 10-59-05 output";
root_input_path = "D:\\data\2019-02-22 Lumen stained cell shape analysis\\Ml1bEGFP\\Cropped\\2019-03-01 10-18-42 output"
root_output_path = "D:\\analysis\\Data analysed for cell shape analysis\\";


root_input_path = io.projection_input_folder_chooser(root_input_path);
root_output_path = io.projection_output_folder_chooser(root_output_path);
fs = [f for f in os.listdir(root_input_path) if (os.path.splitext(f)[1]=='.tif')];

for f in fs:
	im_path = os.path.join(root_input_path, f);
	metadata_path = os.path.join(root_input_path, os.path.splitext(f)[0] + '.json');
	output_path = root_output_path; # just to make reversion easier...
	print("im_path = {}".format(im_path));
	print("metadata_path = {}".format(metadata_path));
	print("output_path = {}".format(output_path));
#	continue;
	
	# basic test data
	#im_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.tif";
	#metadata_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.json";
	#output_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\straightening output";
	#load_centers_path = os.path.join(output_path, os.path.splitext(os.path.basename(im_path))[0], "vessel_axis.csv");
	
	# working data
	#im_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP\\Cropped\\2019-03-01 10-59-05 output\\Cropped dextran-rhodamine lynEGFP 2019-02-22 cell shape analysis control e1d xISV 1.tif"
	#metadata_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP\\Cropped\\2019-03-01 10-59-05 output\\Cropped dextran-rhodamine lynEGFP 2019-02-22 cell shape analysis control e1d xISV 1.json"
	#output_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP\\shape analysis out";
	#load_centers_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP\\shape analysis out\\vessel_axis.csv"
	load_centers_path = None;
	
	output_path = os.path.join(output_path, os.path.splitext(os.path.basename(im_path))[0]);
	if not os.path.isdir(output_path):
		os.mkdir(output_path);
	load_centers_path = os.path.join(output_path, "vessel_axis.csv");
	if not os.path.isfile(load_centers_path):
		load_centers_path = None;
	
	print("starting tube fitting for vessel axis extraction...");
	t1 = datetime.now();
	tt1 = t1;
	imp, info, vessel_centres, max_axis, min_axis = tbf.do_tubefitting(im_path, metadata_path, output_path, load_centers_path=load_centers_path);
	t2 = datetime.now();
	print("done tubefitting in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
	print(imp);
	Prefs.blackBackground = True;
	print("starting angular projection of vessel surface...");
	t1 = datetime.now();
	# TODO: set limits based on actual (predicted) vessel radii in microns
	#if "isv" in info.get_vessel_type().lower():
	#	max_r_pix=int(1.2 * max_axis/2);
	#	min_r_pix=int(0.8 * min_axis/2);
	#else:
	#	max_r_pix = int(1.2 * max_axis/2);
	#	min_r_pix=int(0.8 * min_axis/2);
	max_r_pix=int(1.2 * max_axis/2);
	min_r_pix=int(0.8 * min_axis/2);
	print("max_r_pix = {}, min_r_pix = {}".format(max_r_pix, min_r_pix));
	smooth_radius_um = 1.0;
	out_imp, merged_radius_imp, ring_rois, centres = prj.do_angular_projection(imp, output_path, vessel_centres, generate_roi_stack=True, max_r_pix=max_r_pix, min_r_pix=min_r_pix, smooth_radius_pix=int(smooth_radius_um/info.get_xy_pixel_size_um()));
	print("Memory usage = {:.2E}".format(IJ.currentMemory()));
	imp.changes = False;
	imp.close();
	print("Memory usage = {:.2E}".format(IJ.currentMemory()));
	FileSaver(merged_radius_imp).saveAsTiffStack(os.path.join(output_path, "flythrough with radius identified.tif"));
	merged_radius_imp.changes = False;
	merged_radius_imp.close();
	t2 = datetime.now();
	print("done projection in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
	print("starting unwrapping of cell from around vessel wall...");
	t1 = datetime.now();
	unwrapped_projection_imp, unwrap_axis = prj.twist_and_unwrap(out_imp);
	IJ.setAutoThreshold(unwrapped_projection_imp, "Intermodes dark");
	threshold_val = unwrapped_projection_imp.getProcessor().getMinThreshold();
	print("threshold_val = {}".format(threshold_val));
	FileSaver(unwrapped_projection_imp).saveAsTiff(os.path.join(output_path, "unwrapped cell.tif"));
	t2 = datetime.now();
	print("done unwrapping in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
	print("start generation of image colormapped by radius...");
	t1 = datetime.now();
	r_imp, mask_imp, r_imp_unmasked = prj.generate_r_image(out_imp, unwrapped_projection_imp, ring_rois, centres, unwrap_axis, threshold_val, smooth_radius_pix=int(smooth_radius_um/info.get_xy_pixel_size_um()));
	FileSaver(r_imp).saveAsTiff(os.path.join(output_path, "radius image.tif"));
	FileSaver(r_imp_unmasked).saveAsTiff(os.path.join(output_path, "radius image before masking.tif"));
	FileSaver(mask_imp).saveAsTiff(os.path.join(output_path, "binary cell.tif"));
	t2 = datetime.now();
	print("done generating r image in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
	r_imp.show();
	out_imp.changes = False;
	out_imp.close();
	unwrapped_projection_imp.show();
	print("start calculation of cell inner surface area and aspect ratio...");
	t1 = datetime.now();
	area, aspect_ratio, min_r, max_r = prj.calculate_area_and_aspect_ratio(r_imp, mask_imp, info.get_xy_pixel_size_um());
	print("Area = {}, \n Aspect ratio = {}, \n Min radius = {}, \n, Max radius = {}\n".format(area, aspect_ratio, min_r, max_r));
	t2 = datetime.now();
	print("done calculating outputs in {} s. ".format(round(0.001*(t2.getTime() - t1.getTime()))));
	print("A = " + str(area));
	print("Aspect ratio = " + str(aspect_ratio));
	io.save_projection_csv(output_path, info, [area, aspect_ratio, min_r, max_r]);
	r_imp.changes = False;
	r_imp.close();
	mask_imp.changes = False;
	mask_imp.close();
	unwrapped_projection_imp.changes = False;
	unwrapped_projection_imp.close();
	t2 = datetime.now();
	print("Finished {} in {} s. ".format(f, round(0.001*(t2.getTime() - tt1.getTime()))));