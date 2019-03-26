import math, os, re
from ij import IJ, ImagePlus, Prefs, ImageStack
from ij.gui import OvalRoi, PolygonRoi, WaitForUserDialog
from ij.gui import WaitForUserDialog
from ij.process import FloatProcessor, StackProcessor, ByteProcessor
from ij.plugin import ChannelSplitter, RGBStackMerge, Duplicator, HyperStackConverter

def convert_multichannel_stack_to_Xbit(imp, bitdepth=8):
	"""convert multichannel z stacks to X bits"""
	if bitdepth not in [8, 16, 32]:
		raise NotImplementedError;
	if imp.getType()==ImagePlus.GRAY32:
		if bitdepth==32:
			return imp;
	elif imp.getType()==ImagePlus.GRAY16:
		if bitdepth==16:
			return imp;
	elif imp.getType()==ImagePlus.GRAY8:
		if bitdepth==8:
			return imp;
	print("Memory usage = {:.2E}".format(IJ.currentMemory()));
	current_slice = imp.getCurrentSlice();
	current_channel = imp.getChannel();
	n_c = imp.getNChannels();
	n_z = imp.getNSlices();
	progress_inc = 1 if n_z/20 < 1 else n_z/20;
	stack = imp.getStack();
	out_stack = ImageStack(imp.getWidth(), imp.getHeight());
	for idx in range(n_z * n_c):
		label = stack.getSliceLabel(1);
		ip = stack.getProcessor(1);
		stack.deleteSlice(1);
		if bitdepth==16:
			out_stack.addSlice(label, ip.convertToShort(False));
		elif bitdepth==8:
			out_stack.addSlice(label, ip.convertToByte(False));
		elif bitdepth==32:
			out_stack.addSlice(label, ip.convertToFloat(False));
		else:
			raise NotImplementedError;
		IJ.showProgress(float(idx)/(n_z * n_c));
		if idx%progress_inc==0:
			IJ.showStatus("Converting stack to 16-bits: {}/{}".format(idx, n_z));
	title = imp.getTitle();
	luts = imp.getLuts();
	#print("luts = {}".format(luts));
	imp.changes = False;
	imp.close();
	imp = ImagePlus(title +  " 16-bit converted", out_stack);
#	for ch in range(n_c):
#		imp.setC(ch+1);
#		imp.setLut(luts[ch]);
#	imp.setDisplayMode(IJ.COMPOSITE);
	if n_c>1:
		imp = HyperStackConverter.toHyperStack(imp, n_c, n_z, 1, "Composite");
	imp.show();
	IJ.showProgress(0);
	print("NEW Memory usage = {:.2E}".format(IJ.currentMemory()));
	return imp;

def rot_around_x(input_stack):
	"""do rotation around x axis"""
	output_slices = input_stack.getHeight();
	output_width = input_stack.getWidth();
	output_height = input_stack.getSize();
	output_stack = ImageStack(output_width, output_height);
	for yidx in range(input_stack.getHeight()):
		output_stack.addSlice(FloatProcessor(output_width, output_height, input_stack.getVoxels(0, yidx, 0, output_width, 1, output_height, [])));
	return output_stack;

def rot3d(imp, axis='x'):
	"""pare back Slicer to implement whole-image, +90* rotations that return ImagePlus"""
	if imp.getType()==ImagePlus.COLOR_256 or imp.getType()==ImagePlus.COLOR_RGB:
		raise NotImplementedError("Handling of colour images isn't implemented yet");

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

Prefs.blackBackground = True;

lateral_side_um = 0.01;
longitudinal_side_um = 1.0;

vessel_r_um = 5;
vessel_l_um = 100;
vessel_r_pix = int(vessel_r_um/lateral_side_um);
vessel_l_pix = int(vessel_l_um/longitudinal_side_um);

imp = IJ.createHyperStack("TestData", 4 * vessel_r_pix, 4 * vessel_r_pix, 2, vessel_l_pix, 1, 16);
imp.show();

n_ch = imp.getNChannels();
n_z = imp.getNSlices();
offset_x = (imp.getWidth() - 2*vessel_r_pix)/2;
offset_y = (imp.getHeight() - 2*vessel_r_pix)/2;
centre_x = imp.getWidth()/2;
centre_y = imp.getHeight()/2;

lumen_roi = OvalRoi(offset_x, offset_y, 2*vessel_r_pix, 2*vessel_r_pix);
cell_roi = PolygonRoi([centre_x + (1+vessel_r_pix) * math.cos(math.radians(theta)) for theta in range(180)], 
					[centre_y + (1+vessel_r_pix) * math.sin(math.radians(theta)) for theta in range(180)], 
					PolygonRoi.POLYLINE);

bp = ByteProcessor(imp.getWidth(), imp.getHeight());
cell_slice_imp = ImagePlus("cell slice", bp);
cell_slice_imp.setRoi(cell_roi);
IJ.run(cell_slice_imp, "Line to Area", "");
IJ.run(cell_slice_imp, "Set...", "value={} slice".format(imp.getProcessor().maxValue()));
cell_slice_imp.show();
cell_slice_imp.killRoi();
IJ.run(cell_slice_imp, "Dilate", "");
IJ.run(cell_slice_imp, "Create Selection", "");
#IJ.run(cell_slice_imp, "Make Inverse", "");
cell_roi = cell_slice_imp.getRoi();
cell_slice_imp.changes = False;
cell_slice_imp.close();

for zidx in range(n_z):
	print(zidx);
	imp.setC(1);
	imp.setZ(zidx+1);
	imp.setRoi(lumen_roi);
	IJ.run(imp, "Set...", "value={} slice".format(float(imp.getProcessor().maxValue())/2));
	imp.setC(2);
	imp.setRoi(cell_roi);
	IJ.run(imp, "Rotate...", "rotate angle={}".format(zidx));
	IJ.run(imp, "Set...", "value={} slice".format(float(imp.getProcessor().maxValue())/2));

imp.killRoi();
imp = rot3d(imp, axis='x');
imp.show();
imp = convert_multichannel_stack_to_Xbit(imp, bitdepth=16);
imp.show();
IJ.run(imp, "Size...", "width=200 height=1000 depth=100 average interpolation=Bilinear");