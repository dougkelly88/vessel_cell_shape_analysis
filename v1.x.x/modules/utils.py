import os, re
from datetime import datetime
from ij import ImagePlus, IJ, ImageStack
from ij.gui import WaitForUserDialog
from ij.process import FloatProcessor, StackProcessor
from ij.plugin import ChannelSplitter, RGBStackMerge, Duplicator

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

def factor_tuples(n, lowest_factor=1, highest_small_factor=None, force_large_factor_even=False):
	"""return a list of tuples that are factors of n
	
	Keyword arguments:
	lowest_factor:				minimum value of factors to consider, used to ignore trivial cases (default 1)
	highest_small_factor:		set upper bound on the value of the smaller factor in factor pairs (default None)
	force_large_factor_even:	return only tuples in which the large factor is even, in reverse order - useful for 
								restricting how many images to load in case of two-channel image stacks
	"""
	f = 2 if force_large_factor_even else 1;
	if highest_small_factor is not None:
		return [(i, n//i) for i in range(highest_small_factor, lowest_factor-1, -1) if n % i == 0 and n//i % f == 0];
	else:
		return [(i, n//i) for i in range(lowest_factor, int(n**0.5)+1) if n % i == 0 and n//i % f == 0];

def rename_files(folder, previous_text, new_text):
	"""save function for quickly renaming files if iQ3 mangles filenames"""
	#folder = ("D:\\data\\2019-02-14 Lumen stained cell shape analysis\\" + 
	#		"Control\\e2\\Lyn Control1020190214_175250_20190214_181736");
	#previous_text = "Lyn Control10";
	#new_text = "2dpf TL inj lyn-EGFP + dextran-rhodamine e2d";
	files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and previous_text in f];
	for f in files:
		os.rename(os.path.join(folder, f), os.path.join(folder, f.replace(previous_text, new_text)));
	return;

def cross_planes_approx_median_filter(stack_imp, filter_radius_um=5.0):
	"""relatively computationally cheap, slightly crude approximation of median filter"""
	title = imp.getTitle();
	xy_imp = Duplicator().run(stack_imp);
	xy_imp.setTitle("xy {}".format(title));
	xz_imp = Duplicator().run(stack_imp);
	xz_imp.setTitle("xz {}".format(title));
	xz_imp = rot3d(xz_imp, 'x');
	zy_imp = Duplicator().run(stack_imp);
	zy_imp.setTitle("zy {}".format(title));
	zy_imp = rot3d(zy_imp, 'y');
	stack_imp.changes = False;
	stack_imp.close();
	xy_imp = stack_median_filter(xy_imp, radius_um=filter_radius_um);
	xz_imp = stack_median_filter(xz_imp, radius_um=filter_radius_um);
	zy_imp = stack_median_filter(zy_imp, radius_um=filter_radius_um);
	xz_imp = rot3d(xz_imp, 'x');
	zy_imp = rot3d(zy_imp, '-y');
	ic = ImageCalculator();
	dummy = ic.run("Add create 32-bit stack", xy_imp, xz_imp);
	xz_imp.close();
	xy_imp.close();
	output_imp = ic.run("Add create 32-bit stack", dummy, zy_imp);
	zy_imp.close();
	dummy.close();
	output_imp.show();
	IJ.run(output_imp, "Divide...", "value=3 stack");
	output_imp.setTitle("Cross-plane median filtered {} (r={} um)".format(title, filter_radius_um).replace(".tif", ""));
	return output_imp;

def downsample_for_isotropy(imp, extra_downsample_factor=2.0):
	"""downsample x, y pixel directions to get cubic voxels"""
	title = imp.getTitle();
	cal = imp.getCalibration();
	pix_w = cal.pixelWidth;
	pix_h = cal.pixelHeight;
	pix_d = cal.pixelDepth;
	im_w = imp.getWidth();
	im_h = imp.getHeight();
	im_d = imp.getNSlices();
	print("downsampling {} and making isotropic...".format(title));
	IJ.showStatus("Downsampling and making ~isotropic...");
	xy_scale = pix_h / (pix_d * extra_downsample_factor);
	xy_scaled_h = int(xy_scale * im_h);
	xy_scaled_w = int(xy_scale * im_w);
	z_scale = 1/ extra_downsample_factor;
	z_scaled_h = int(z_scale * im_d);
	sp = StackProcessor(imp.getStack());
	stack = sp.resize(xy_scaled_w, xy_scaled_h, True);
	xz_stack = rot_around_x(stack);
	xz_sp = StackProcessor(xz_stack);
	xz_stack = xz_sp.resize(xy_scaled_w, z_scaled_h, True);
	out_stack = rot_around_x(xz_stack);
	cal.setUnit('um');
	cal.pixelWidth = im_w * pix_w/xy_scaled_w;
	cal.pixelHeight = im_h * pix_h/xy_scaled_w;
	cal.pixelDepth = im_d * pix_d/z_scaled_h;
	imp.changes = False;
	imp.close();
	out_imp = ImagePlus("Isotropic downsampled {}".format(title), out_stack);
	out_imp.setCalibration(cal);
	print("...done downsampling {} and making isotropic. ".format(title));
	IJ.showStatus("...done downsampling and making ~isotropic. ");
	return out_imp;

def stack_median_filter(imp, radius_um=5.0):
	"""perform 2D median filtering slicewise on a stack"""
	t1 = datetime.now();
	stack = imp.getStack();
	title = imp.getTitle();
	n_z = imp.getNSlices();
	progress_inc = 1 if n_z/20 < 1 else n_z/20;
	filt = RankFilters();
	for idx in range(n_z):
		#print("Median filtering stack, {} % complete...".format(100 * round(float(zidx)/n_z, 3)));
		ip = stack.getProcessor(idx+1);
		filt.rank(ip, radius_um, RankFilters.MEDIAN);
		IJ.showProgress(float(idx)/n_z);
		if idx%progress_inc==0:
			IJ.showStatus("Median filtering stack in 2D: {}/{}".format(idx, n_z));
	imp.updateAndDraw();
	imp.setTitle("Median filtered (r={}) {}".format(radius_um, imp.getTitle()).replace(".tif", ""));
	t2 = datetime.now();
	IJ.showProgress(1.0);
	print("Median filtering took {} s".format((t2-t1).total_seconds()));
	return imp;

def threshold_and_binarise(imp):
	"""Return thresholded stack"""
	print("performing segmentation on channel: " + imp.getTitle());
	stats = StackStatistics(imp);
	IJ.run(imp, "Subtract...", "value={} stack".format(stats.min));
	IJ.run(imp, "Divide...", "value={} stack".format((stats.max - stats.min)/255.0));
	WaitForUserDialog("pre 8-bit").show();
	imp = robust_convertStackToGray8(imp);
	WaitForUserDialog("post 8-bit").show();

	# Apply automatic THRESHOLD to differentiate cells from background
	histo = StackStatistics(imp).histogram;
	thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Default, histo);
	print(thresh_lev);
	min_volume = int(float(imp.getHeight() * imp.getWidth() * imp.getNSlices())/100);
	IJ.run(imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
												" min_size=" + str(min_volume) + " max_size=-1");
	fit_basis_imp = WM.getImage("Seg");
	bin_imp = WM.getImage("Bin");
	bin_imp.changes=False;
	bin_imp.close();
	IJ.setThreshold(fit_basis_imp, 1, 65535);
	IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");
	#IJ.run(fit_basis_imp, "Fill Holes", "stack");
#	IJ.run("3D Fill Holes", "");
	return fit_basis_imp;

def robust_convertStackToGrayX(imp, x=8):
	"""simplified from https://github.com/imagej/imagej1/blob/master/ij/process/StackConverter.java, 
	avoiding complications of scaling based on LUT taken from the current frame. Assumes conversion
	from grayscale image"""
	if x!=8 and x!=16 and x!=32:
		raise NotImplementedError("can't convert to the specified bit depth");
	if imp.getBitDepth()==x:
		return imp;
	current_slice = imp.getCurrentSlice();
	stack = imp.getStack();
	out_stack = ImageStack(imp.getWidth(), imp.getHeight());
	n_z = imp.getNSlices();
	progress_inc = 1 if n_z/20 < 1 else n_z/20;
	for idx in range(n_z):
		# always use 1 as old stack index since we remove a slice at each iteration
		label = stack.getSliceLabel(1);
		ip = stack.getProcessor(1);
		stack.deleteSlice(1);
		if x==8:
			out_stack.addSlice(label, ip.convertToByte(False));
		elif x==16:
			out_stack.addSlice(label, ip.convertToShort(False));
		elif x==32:
			out_stack.addSlice(label, ip.convertToFloat());
		IJ.showProgress(float(idx)/n_z);
		if idx%progress_inc==0:
			IJ.showStatus("Converting stack to {}-bits: {}/{}".format(x, idx, n_z));
	imp.setStack(out_stack);
	IJ.showProgress(1.0);
	imp.setSlice(current_slice);
	return imp;

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
