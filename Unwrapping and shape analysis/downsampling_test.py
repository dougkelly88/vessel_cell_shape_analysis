import os, re, sys
from datetime import datetime
from ij import ImagePlus, IJ, ImageStack
from ij.gui import WaitForUserDialog
from ij.process import FloatProcessor, StackProcessor
from ij.plugin import ChannelSplitter, RGBStackMerge, Duplicator

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

from UpdateRoiImageListener import UpdateRoiImageListener
from PrescreenInfo import PrescreenInfo
import file_io as io
import ellipse_fitting
import ui


def rot_around_x(input_stack):
	"""do rotation around x axis"""
	output_slices = input_stack.getHeight();
	output_width = input_stack.getWidth();
	output_height = input_stack.getSize();
	output_stack = ImageStack(output_width, output_height);
	for yidx in range(input_stack.getHeight()):
		IJ.showProgress(float(yidx)/output_slices);
		output_stack.addSlice(FloatProcessor(output_width, output_height, input_stack.getVoxels(0, yidx, 0, output_width, 1, output_height, [])));
	IJ.showProgress(1.0);
	return output_stack;

def rot3d(imp, axis='x'):
	"""pare back Slicer to implement whole-image, +90* rotations that return ImagePlus"""
	if imp.getType()==ImagePlus.COLOR_256 or imp.getType()==ImagePlus.COLOR_RGB:
		raise NotImplementedError("Handling of colour images isn't implemented yet");
	IJ.showStatus("Rotating around {}-axis...".format(axis))
	title = imp.getTitle();
	original_cal = imp.getCalibration();
	new_cal = original_cal;
	
	if imp.getNChannels > 1:
		split_ch = ChannelSplitter().split(imp);
	else:
		split_ch = [imp];

	out_imps = [];
	if axis=='x':
		for ch_imp in split_ch:
			input_stack = ch_imp.getStack();
			output_stack = rot_around_x(input_stack);
			out_imps.append(ImagePlus(title, output_stack));
			new_cal.pixelHeight = original_cal.pixelDepth;
			new_cal.pixelDepth = original_cal.pixelHeight;
	elif axis=='y' or axis=='-y':
		for ch_imp in split_ch:
			if axis[0]=='-':
				input_stack = StackProcessor(ch_imp.getStack()).rotateLeft();
			else:
				input_stack = StackProcessor(ch_imp.getStack()).rotateRight();
			output_stack = rot_around_x(input_stack);
			if axis[0]=='-':
				final_stack = StackProcessor(output_stack).rotateLeft();
			else:
				final_stack = StackProcessor(output_stack).rotateRight();
			out_imps.append(ImagePlus(title, final_stack));
			new_cal.pixelWidth = original_cal.pixelDepth;
			new_cal.pixelDepth = original_cal.pixelWidth;
	elif axis=='z' or axis=='-z':
		for ch_imp in split_ch:
			if axis[0]=='-':
				output_stack = StackProcessor(ch_imp.getStack()).rotateLeft();
			else:
				output_stack = StackProcessor(ch_imp.getStack()).rotateRight();
			out_imps.append(ImagePlus(title, output_stack));
			new_cal.pixelWidth = original_cal.pixelHeight;
			new_cal.pixelHeight = original_cal.pixelWidth;
	else:
		raise NotImplementedError("Please check which axis you've chosen - if not (x, +/-y, +/-z) then it's not implemented...");
	imp.changes = False;
	imp.close();
	if len(out_imps) > 1:
		out_imp = RGBStackMerge().mergeChannels(out_imps, False);
	else:
		out_imp = out_imps[0];
	out_imp.setCalibration(new_cal);
	return out_imp;

def downsample_for_isotropy(imp, extra_downsample_factor=2.0, info=None):
	"""downsample x, y pixel directions to get cubic voxels"""
	title = imp.getTitle();
	cal = imp.getCalibration();
	if info is None:
		pix_w = cal.pixelWidth;
		pix_h = cal.pixelHeight;
		pix_d = cal.pixelDepth;
	else:
		pix_w = info.get_xy_pixel_size_um();
		pix_h = pix_w;
		pix_d = info.get_z_plane_spacing_um();
		print("pixel whd = ({}, {}, {}".format(pix_w, pix_h, pix_d));
	im_w = imp.getWidth();
	im_h = imp.getHeight();
	im_d = imp.getNSlices();
	im_nch = imp.getNChannels();
	if im_nch > 1:
		split_ch = ChannelSplitter().split(imp);
	else:
		split_ch = [imp];
	print("downsampling {} and making isotropic...".format(title));
	IJ.showStatus("Downsampling and making ~isotropic...");
	xy_scale = pix_h / (pix_d * extra_downsample_factor);
	xy_scaled_h = int(xy_scale * im_h);
	xy_scaled_w = int(xy_scale * im_w);
	z_scale = 1/ extra_downsample_factor;
	z_scaled_h = int(z_scale * im_d);
	out_imps = [];
	for ch_imp in split_ch:
		print(ch_imp.getTitle());
		sp = StackProcessor(ch_imp.getStack());
		print((xy_scaled_w, xy_scaled_h));
		stack = sp.resize(xy_scaled_w, xy_scaled_h, True);
		xz_stack = rot_around_x(stack);
		xz_sp = StackProcessor(xz_stack);
		xz_stack = xz_sp.resize(xy_scaled_w, z_scaled_h, True);
		out_stack = rot_around_x(xz_stack);
		out_imps.append(ImagePlus("Isotropic downsampled {}".format(title), out_stack));
	cal.setUnit('um');
	cal.pixelWidth = im_w/xy_scaled_w * pix_w;
	cal.pixelHeight = im_h/xy_scaled_h * pix_h;
	cal.pixelDepth = im_d/z_scaled_h * pix_d;
	imp.changes = False;
	imp.close();
	for ch_imp in split_ch:
		ch_imp.close();
	if len(out_imps) > 1:
		out_imp = RGBStackMerge().mergeChannels(out_imps, False);
	else:
		out_imp = out_imps[0];
	out_imp.setCalibration(cal);
	print("...done downsampling {} and making isotropic. ".format(title));
	IJ.showStatus("...done downsampling and making ~isotropic. ");
	return out_imp;

test_path = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Ml1bEGFP\\Cropped\\2019-03-01 10-18-42 output";
info = PrescreenInfo().load_info_from_json(os.path.join(test_path, "Cropped dextran-rhodamine ml1b-egfp 2019-02-27 cell shape analysis ee2b xISV 1.json"));
imp = IJ.openImage(os.path.join(test_path, "Cropped dextran-rhodamine ml1b-egfp 2019-02-27 cell shape analysis ee2b xISV 1.tif"));
out_imp = downsample_for_isotropy(imp, extra_downsample_factor=2.0);
out_imp.show()