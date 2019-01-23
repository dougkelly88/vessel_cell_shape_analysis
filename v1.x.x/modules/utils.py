from ij import ImagePlus, IJ, ImageStack
from ij.gui import WaitForUserDialog
from ij.process import FloatProcessor, StackProcessor

def rot_around_x(input_stack):
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
	if axis=='x':
		input_stack = imp.getStack();
		output_stack = rot_around_x(input_stack);
		out_imp = ImagePlus(title, output_stack)
	elif axis=='y':
		input_stack = StackProcessor(imp.getStack()).rotateRight();
		output_stack = rot_around_x(input_stack);
		final_stack = StackProcessor(output_stack).rotateRight();
		out_imp = ImagePlus(title, final_stack);
	elif axis=='z':
		output_stack = StackProcessor(imp.getStack()).rotateRight();
		out_imp = ImagePlus(title, output_stack)
	else:
		raise NotImplementedError("Please check which axis you've chosen - if not (x, y, z) than it's not implemented...");
	imp.close();
	return out_imp;

path = "C:\\Users\\dougk\\Desktop\\test image.tif";
imp = IJ.openImage(path);
imp.show();

axis = 'y'
out_imp = rot3d(imp, axis=axis);
out_imp.show();
WaitForUserDialog("rotated around " + axis).show();
#out_imp2.close();
#imp.close();
