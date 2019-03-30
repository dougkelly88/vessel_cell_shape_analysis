import os, sys
from ij import IJ
from ij.gui import WaitForUserDialog
from ij.plugin import Duplicator, ZProjector

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

from PrescreenInfo import PrescreenInfo
import file_io as io;
	
def z_crop(imp):
	"""trim a z stack based on interactively-defined start and end points"""
	IJ.setTool("zoom");
	IJ.run("Brightness/Contrast...");
	imp.setZ(1);
	WaitForUserDialog("Choose first z plane and click OK...").show();
	start_z = imp.getZ();
	WaitForUserDialog("Now choose last z plane and click OK...").show();
	end_z = imp.getZ();
	frames = imp.getNFrames();
	channels = imp.getNChannels();
	imp.killRoi();
	dupimp = Duplicator().run(imp, 1, channels, start_z, end_z, 1, frames);
	imp.hide();
	dupimp.show()
	return dupimp, (start_z, end_z);

def perform_cropping(imp, info, output_folder, default_path):
	imp.show();
	print("C = {}, Z = {}, T = {}".format(imp.getNChannels(), imp.getNSlices(), imp.getNFrames()));
	imp.setC(info.get_mosaic_labeled_ch());
	IJ.run("Enhance Contrast", "saturated=0.35");
	zcrop_imp, z_lims = z_crop(imp);
	info.set_z_crop_frames(z_lims);
	info.set_mosaic_labeled_ch(imp.getC());
	IJ.setTool("rect");
	zcrop_imp.hide();
	imp2 = ZProjector.run(zcrop_imp, "max");
	imp2.show();
	imp2.setC(info.get_mosaic_labeled_ch());
	crop_roi = None;
	while crop_roi is None:
		WaitForUserDialog("Select XY region to crop...").show();
		crop_roi = imp2.getRoi();
	if crop_roi.getType():
		info.set_xy_crop_rect(str([(x,y) for x, y in zip(crop_roi.getPolygon().xpoints, crop_roi.getPolygon().ypoints)]));
	info.set_xy_crop_rect(crop_roi.getBounds().toString());
	imp2.close();
	zcrop_imp.show();
	zcrop_imp.setRoi(crop_roi);
	IJ.run("Crop", "");
	if crop_roi.getType():
		IJ.run(zcrop_imp, "Make Inverse", "");
		inv_roi = zcrop_imp.getRoi();
		IJ.run(zcrop_imp, "Set...", "value=0 stack");
		IJ.run(zcrop_imp, "Make Inverse", "");
	return imp, zcrop_imp, info, info.get_input_file_path();
		
def main():
	info = PrescreenInfo();
	default_path = "D:\\data\\Marcksl1 cell shape analysis\\zstacks"
	output_folder = io.output_folder_chooser(default_path);
	used_files = [];
	answer = "";
	imp = None;
	while answer != "stop":
		answer = "stop";
		if imp is None:
			imp, info = io.get_image(info, default_path, used_files);
		imp, zcrop_imp, info, default_path = perform_cropping(imp, info, output_folder, default_path);
		answer = io.save_etc(zcrop_imp, info, output_folder);
		used_files.append(info.get_input_file_path());
		print("answer = " + str(answer));
		if answer == "continue_newimage":
			imp.changes = False;
			imp.close();
			imp = None;
		if zcrop_imp is not None:
			zcrop_imp.changes = False;
			zcrop_imp.close();
	if imp is not None:
		imp.changes = False;
		imp.close();		
	return;
	
# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main();