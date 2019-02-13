from ij import ImagePlus, IJ, ImageStack
from ij.gui import WaitForUserDialog
from ij.process import FloatProcessor, StackProcessor
from ij.plugin import ChannelSplitter, RGBStackMerge

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
	elif axis=='y':
		for ch_imp in split_ch:
			input_stack = StackProcessor(ch_imp.getStack()).rotateRight();
			output_stack = rot_around_x(input_stack);
			final_stack = StackProcessor(output_stack).rotateRight();
			out_imps.append(ImagePlus(title, final_stack));
			new_cal.pixelWidth = original_cal.pixelDepth;
			new_cal.pixelDepth = original_cal.pixelWidth;
	elif axis=='z':
		for ch_imp in split_ch:
			output_stack = StackProcessor(ch_imp.getStack()).rotateRight();
			out_imps.append(ImagePlus(title, output_stack));
			new_cal.pixelWidth = original_cal.pixelHeight;
			new_cal.pixelHeight = original_cal.pixelWidth;
	else:
		raise NotImplementedError("Please check which axis you've chosen - if not (x, y, z) than it's not implemented...");
	imp.changes = False;
	imp.close();
	if len(out_imps) > 1:
		out_imp = RGBStackMerge().mergeChannels(out_imps, False);
	else:
		out_imp = out_imps[0];
	out_imp.setCalibration(new_cal);
	return out_imp;

#path = "C:\\Users\\dougk\\Desktop\\test image 2.tif";
#imp = IJ.openImage(path);
#imp.show();
#
#axis = 'x'
#out_imp = rot3d(imp, axis=axis);
#out_imp.show();
#WaitForUserDialog("rotated around " + axis).show();
#out_imp2.close();
#imp.close();