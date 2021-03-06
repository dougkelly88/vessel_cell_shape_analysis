# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math, sys, os

from ij import IJ, ImageStack, ImagePlus, Prefs
from ij.gui import EllipseRoi, WaitForUserDialog, Roi, PolygonRoi
from ij.io import FileSaver
from ij.plugin import ChannelSplitter, ImageCalculator, Duplicator, SubstackMaker, Straightener, RGBStackMerge
from ij.process import StackProcessor, AutoThresholder, StackStatistics, FloatProcessor
from ij import WindowManager as WM
from loci.plugins import BF as bf


#im_test_path = "D:\\data\\Marcksl1 cell shape analysis\\e27 ISV1.tif";
#metadata_test_path = "D:\\data\Marcksl1 cell shape analysis\\Cropped\\2018-12-05 16-06-29 output\\Cropped UAS-marcksl1b-delED e27 xISV 1.json";
#im_test_path = "D:\\data\\Structural imaging\\Imaging protocol tests\\2018-09-14 vessel structure imaging\\vessel_structure_protocol3_20180914_174448\\Cropped vessels\\2018-12-12 16-38-31 output\\Cropped dextran-rhodamine kdrl-egfp vessel_structure etest_embryo aISV 1.tif";
#metadata_test_path = "D:\\data\\Structural imaging\\Imaging protocol tests\\2018-09-14 vessel structure imaging\\vessel_structure_protocol3_20180914_174448\\Cropped vessels\\2018-12-12 16-38-31 output\\Cropped dextran-rhodamine kdrl-egfp vessel_structure etest_embryo aISV 1.json";
im_test_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.tif";
metadata_test_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.json";
output_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\straightening output";

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
import utils

def generate_smoothed_vessel_axis(centres, pixel_size_um=0.1625, smooth_parameter_um=5.0):
	"""From a list of fitted vessel centres, generate the vessel path"""
	out_centres = [];
	tangent_vectors = [];
	smooth_planes = smooth_parameter_um/pixel_size_um;
	xs = [x for (x,y) in centres];
	ys = [y for (x,y) in centres];
	for idx, centre in enumerate(centres):
		high_idx = idx + int(round(smooth_planes/2));
		high_idx = high_idx if high_idx < len(centres)-1 else len(centres)-2;
		low_idx = idx - int(round(smooth_planes/2));
		low_idx = low_idx if low_idx > 0 else 0;
		out_centres.append((sum(xs[low_idx:high_idx+1])/(high_idx+1-low_idx), sum(ys[low_idx:high_idx+1])/(high_idx+1-low_idx)));
		if idx > 0:
			tangent_vectors.append((out_centres[-1][0] - out_centres[-2][0], out_centres[-1][1] - out_centres[-2][1], 1));
	tangent_vectors.insert(0, tangent_vectors[0]);
	return out_centres, tangent_vectors;

def threshold_and_binarise(imp, z_xy_ratio, approx_median=True, prune_slicewise=True):
	"""Return thresholded stack"""
	print("performing segmentation on channel: " + imp.getTitle());
	if approx_median:
		imp = utils.cross_planes_approx_median_filter(imp);
		imp = utils.robust_convertStackToGrayXbit(imp, x=8, normalise=True);
	else:
		filter_radius = 3.0;
		IJ.run(imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(math.ceil(filter_radius / z_xy_ratio)) + " z=" + str(filter_radius));
		IJ.run(imp, "8-bit", "");
	imp.show();

	# Apply automatic THRESHOLD to differentiate cells from background
	# get threshold value from stack histogram using otsu
	histo = StackStatistics(imp).histogram;
	thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.IJ_IsoData, histo);
	#thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histo);
	max_voxel_volume = int(float(imp.getHeight() * imp.getWidth() * imp.getNSlices())/100);
	IJ.run(imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
												" min_size=" + str(max_voxel_volume) + " max_size=-1");
	fit_basis_imp = WM.getImage("Seg");
	bin_imp = WM.getImage("Bin");
	bin_imp.changes=False;
	bin_imp.close();
	IJ.setThreshold(fit_basis_imp, 1, 65535);
	IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");
	#IJ.run(fit_basis_imp, "Fill Holes", "stack");
	IJ.run("3D Fill Holes", "");
	if prune_slicewise:
		fit_basis_imp = utils.keep_largest_blob(fit_basis_imp); # note this will be unstable if "branches" of approx equal x-section exist
	return fit_basis_imp;

def combine_two_channel_segments(imp1, imp2, binary_imp1, binary_imp2):
	"""Use both channels to segment the vessel..."""
	# two choices - 
	# EITHER: combine median-filtered channels BEFORE segmentation: this means that correlation between 
	# planes is maintained thanks to 3d segmentation
	# OR: make choice about which channel to use (or whether to combine - how?) AFTER segmentation - potentially easier, so start with this
	stack = ImageStack(imp1.getWidth(), imp1.getHeight());
	for zidx in range(1, imp1.getNSlices()+1):
		WaitForUserDialog("about to combine channels for slice " + str(zidx)).show(); 
		print("zidx = " + str(zidx));
		ch1_snr = calculate_snr(imp1, binary_imp1, zidx);
		print("Ch1 snr = " + str(ch1_snr));
		ch2_snr = calculate_snr(imp2, binary_imp2, zidx);
		print("Ch2 snr = " + str(ch2_snr));
		if (0.75 * ch1_snr) > ch2_snr:
			binary_imp1.setZ(zidx);
			ip = binary_imp1.getProcessor();
		else:
			print("USING MCH CHANNEL!");
			binary_imp2.setZ(zidx);
			ip = binary_imp2.getProcessor();
		stack.addSlice(ip);
	return ImagePlus("Combined thresholded channels", stack);

def calculate_snr(signal_imp, thresholded_imp, zidx):
	"""get the snr, calculated as the average signal in the thresholded region vs the average signal outwith the thresholded region, for a given slice"""
	signal_imp.setZ(zidx);
	thresholded_imp.setZ(zidx);
	IJ.run(thresholded_imp, "Create Selection", "");
	roi = thresholded_imp.getRoi();
	signal_imp.setRoi(roi);
	signal = signal_imp.getStatistics().mean;
	IJ.run(signal_imp, "Make Inverse", "");
	noise = signal_imp.getStatistics().mean;
	return float(signal)/noise;

def normalise_to_fill_range(imp, max_val=255):
	"""Return an imp with minimum=0 and maximum=max_val"""
	stats = StackStatistics(imp);
	IJ.run(imp, "Subtract...", "value=" + str(stats.min) + " stack");
	IJ.run(imp, "Multiply...", "value=" + str(float(max_val)/(stats.max - stats.min)) + " stack");
	return imp;

def convex_hull_pts(pts):
	"""Return points describing effective convex hull from non-contiguous outline points"""
	clean_pts = [];
	temp_pts = [];
	ys = [y for (x,y) in pts]
	for yy in set(ys):
	    xs = sorted([x for (x,y) in pts if y==yy]);
	    temp_pts.append((xs[0], yy));
	    if len(xs)>1:
	        temp_pts.append((xs[-1], yy));
	xs = [x for (x,y) in temp_pts];
	for xx in set(xs):
	    ys = sorted([y for (x,y) in temp_pts if x==xx]);
	    clean_pts.append((xx, ys[0]));
	    if len(xs)>1:
	        clean_pts.append((xx, ys[-1]));
	return clean_pts;

def split_and_rotate(imp, info):
	"""return image to segment on, image to project out, and images to display"""
	# for now, assume that these are ISVs and that embryo is mounted in familiar fashion. First of these can be developed out...
	print("Image dimensions at start of split and rotate = ({}x{}x{}) pix".format(imp.getWidth(), 
																			   imp.getHeight(), 
																			   imp.getNSlices()));
	if imp.isVisible():
		IJ.run("Enhance Contrast", "saturated=0.35");
	seg_ch_idx, proj_ch_idx = ui.choose_segmentation_and_projection_channels(info);
	channels  = ChannelSplitter().split(imp);
	seg_imp = Duplicator().run(channels[seg_ch_idx]); # use Duplicator to decouple - can do smarter to save memory?
	proj_imp = Duplicator().run(channels[proj_ch_idx]);
	rot_seg_imp = utils.rot3d(seg_imp, axis='x');
	rot_seg_imp.setTitle("rot_seg_imp");
	rot_proj_imp = utils.rot3d(proj_imp, axis='x');
	rot_proj_imp.setTitle("rot_proj_imp");
	egfp_mch_imps = [];
	egfp_idx = 0 if "gfp" in info.ch1_label.lower() else 1;
	mch_idx = int(not(egfp_idx)); # assume two channel...
	for ch_idx in [egfp_idx, mch_idx]:
		if ch_idx==seg_ch_idx:
			egfp_mch_imps.append(Duplicator().run(rot_seg_imp));
		elif ch_idx==proj_ch_idx:
			egfp_mch_imps.append(Duplicator().run(rot_proj_imp));
		else:
			egfp_mch_imps.append(utils.rot3d(Duplicator().run(channels[ch_idx]), axis='x'));
	imp.changes=False;
	imp.close();
	seg_imp.changes = False;
	proj_imp.changes = False;
	seg_imp.close();
	proj_imp.close();
	print("Image dimensions at start of split and rotate = ({}x{}x{}) pix".format(rot_proj_imp.getWidth(), 
																			   rot_proj_imp.getHeight(), 
																			   rot_proj_imp.getNSlices()));
	return rot_seg_imp, rot_proj_imp, egfp_mch_imps

def lin_interp_1d(old_x, old_y, new_x):
	"""stretch a curve described by points [(old_x, old_y)] to cover domain given by [new_x]"""
	xpp = [o * float(len(old_x)-1)/(len(new_x)-1) for o in new_x]
	new_y = [old_y[0]];
	for o in xpp:
	    if (o!=xpp[0]) & (o!=xpp[-1]):
	        xa = int(math.floor(o));
	        xb = int(math.ceil(o));
	        ya = old_y[xa];
	        yb = old_y[xb];
	        new_y.append(ya + (yb - ya)*(o - xa)/(xb - xa));
	new_y.append(old_y[-1]);
	return new_y;
	
def straighten_vessel(imp, smooth_centres, it=1, save_output=False):
	"""use IJ straigtening tool to deal with convoluted vessels"""
	print("straighten vessel input image dims = " + str(imp.getWidth()) + "x" + str(imp.getHeight()));
	rot_imp = utils.rot3d(imp, axis='x');
	if it==1:
		roi = PolygonRoi([x for x, y, z in smooth_centres], [z for x, y, z in smooth_centres], Roi.FREELINE);
		print("len interp polygon = " + str(roi.getInterpolatedPolygon().npoints));
	elif it==2:
		new_zs = [z for z in range(rot_imp.getWidth())];
		new_ys = lin_interp_1d([z for x, y, z in smooth_centres], [y for x, y, z in smooth_centres], new_zs);
		roi = PolygonRoi(new_zs, new_ys, Roi.FREELINE);
	
	split_ch = ChannelSplitter().split(rot_imp);
	mch_imp = split_ch[0];
	egfp_imp = split_ch[1];
	roi_imp = split_ch[2];

	roi_imp.setRoi(roi);

	for zidx in range(egfp_imp.getNSlices()):
		for chidx in range(3):
			split_ch[chidx].setZ(zidx+1);
			split_ch[chidx].setRoi(roi);
			ip = Straightener().straightenLine(split_ch[chidx], 150);
			if chidx==1:
				if zidx==0:
					egfp_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
				egfp_straight_stack.addSlice(ip);
			elif chidx==0:
				if zidx==0:
					mch_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
				mch_straight_stack.addSlice(ip);
			else:
				if zidx==0:
					roi_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
				roi_straight_stack.addSlice(ip);

	egfp_out_imp = ImagePlus("Straightened EGFP", egfp_straight_stack);
	mch_out_imp = ImagePlus("Straightened mCh", mch_straight_stack);
	roi_out_imp = ImagePlus("Straightened ROI", roi_straight_stack);
	if it==2: 
		egfp_out_imp = utils.rot3d(egfp_out_imp, axis='y');
		mch_out_imp = utils.rot3d(mch_out_imp, axis='y');
		roi_out_imp = utils.rot3d(roi_out_imp, axis='y');
	egfp_out_imp.show();
	mch_out_imp.show();
	roi_out_imp.show();
	IJ.run("Merge Channels...", "c1=[" + mch_out_imp.getTitle() + 
										"] c2=[" + egfp_out_imp.getTitle() + 
										"] c7=[" + roi_out_imp.getTitle() + "] create keep");
#	WaitForUserDialog("pause").show();
#	if it==1:
	egfp_out_imp.close()
	mch_out_imp.close()
	roi_out_imp.close()
	new_composite = IJ.getImage();
	if save_output:
		FileSaver(new_composite).saveAsTiffStack(os.path.join(output_path, "after rotation " + str(it) + ".tif"));
	return new_composite;

def do_tubefitting(im_path=im_test_path, metadata_path=metadata_test_path, output_path=output_path, save_output=False):
	# todo: fix things so that all operations use a consistent definition of background rather than changing Prefs on the fly...
	Prefs.blackBackground = False;
	info = PrescreenInfo();
	info.load_info_from_json(metadata_path);
	z_xy_ratio = abs(info.get_z_plane_spacing_um()) / info.get_xy_pixel_size_um();
	#z_xy_ratio = 1.0;
	bfimp = bf.openImagePlus(im_path);
	imp = bfimp[0];
	imp.show();
	IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	imp = utils.downsample_for_isotropy(imp, extra_downsample_factor=1.0, info=info);
	rot_seg_imp, rot_proj_imp, egfp_mch_imps = split_and_rotate(imp, info);
	depth = rot_seg_imp.getNSlices() if rot_seg_imp.getNSlices() > rot_seg_imp.getNFrames() else rot_seg_imp.getNFrames();
	width = rot_seg_imp.getWidth();
	height = int(round(rot_seg_imp.getHeight() * z_xy_ratio));

	# Apply 3d MEDIAN FILTER to denoise and emphasise vessel-associated voxels
	fit_basis_imp = threshold_and_binarise(rot_seg_imp, z_xy_ratio);
	fit_basis_imp.setTitle("fit_basis_imp");
	fit_basis_imp.show();

	# plane-wise, use binary-outline
	# say the non-zero points then make up basis for fitting to be performed per http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
	rois = [];
	centres = [];
	major_axes = [];
	roi_imp = IJ.createImage("rois", width, height, depth, 32);
	pts_stack = ImageStack(width, height+1);
	IJ.run(imp, "Line Width...", "line=3");
	for zidx in range(fit_basis_imp.getNSlices()):
		fit_basis_imp.setZ(zidx+1);
		IJ.run(fit_basis_imp, "Outline", "slice");
		IJ.run(fit_basis_imp, "Create Selection", "");
		roi = fit_basis_imp.getRoi();
		fit_basis_imp.killRoi();
		pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()];
		clean_pts = convex_hull_pts(pts);
		clean_pts = [(x, z_xy_ratio * y) for (x, y) in clean_pts];
		# make a stack of clean points...
		ip = FloatProcessor(width, height+1)
		pix = ip.getPixels();
		for pt in clean_pts:
			pix[int(pt[1]) * width + int(pt[0])] = 128;
		pts_stack.addSlice(ip);
		centre, angle, axl = ellipse_fitting.fit_ellipse(clean_pts);
		major_axes.append(max(axl));
		centres.append(centre);
		rot_seg_imp.setZ(zidx+1);
		ellipse_roi = ellipse_fitting.generate_ellipse_roi(centre, angle, axl);
		rois.append(ellipse_roi);
	IJ.run(imp, "Line Width...", "line=1");
	cal = imp.getCalibration();
	smooth_centres, tangent_vecs =  generate_smoothed_vessel_axis(centres, pixel_size_um=cal.pixelDepth);
	for zidx in range(fit_basis_imp.getNSlices()):
		centre = smooth_centres[zidx];
		major_axis = major_axes[zidx];
		ellipse_roi = EllipseRoi(centre[0]-2, centre[1], centre[0]+2, centre[1], 1.0);
		roi_imp.setZ(zidx+1);
		roi_imp.setRoi(ellipse_roi);
		IJ.run(roi_imp, "Set...", "value=" + str(roi_imp.getProcessor().maxValue()) + " slice");

	pts_stack_imp = ImagePlus("Cleaned points", pts_stack);
	pts_stack_imp.setTitle("pts_stack_imp");
	pts_stack_imp.show();

	rot_seg_imp.changes = False;
	rot_seg_imp.close();
	egfp_imp = egfp_mch_imps[0];
	mch_imp = egfp_mch_imps[1];
	imps_to_combine = [egfp_mch_imps[1], egfp_mch_imps[0], roi_imp];
	egfp_imp.show()
	mch_imp.show();
	roi_imp.show();
	print("box height um = " + str(roi_imp.getNSlices() * info.get_xy_pixel_size_um()));
	IJ.run(egfp_imp, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
	IJ.run(mch_imp, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
	#IJ.run("Merge Channels...", "c1=[" + mch_imp.getTitle() + 
	#								"] c2=[" + egfp_imp.getTitle() + 
	#								"] c7=[" + roi_imp.getTitle() + "] create keep");
	composite_imp = RGBStackMerge().mergeChannels(imps_to_combine, False);
	print(composite_imp);
	composite_imp.show();
	print("end of vessel centerline id step, image dims = ({}x{}x{})".format(composite_imp.getWidth(), 
																		  composite_imp.getHeight(), 
																		  composite_imp.getNSlices()));
	WaitForUserDialog("pause").show();
	# do qc here?

	#WM.getImage("Composite").addImageListener(UpdateRoiImageListener(rois));
	IJ.run(roi_imp, "8-bit", "");

	if save_output:
		FileSaver(composite_imp).saveAsTiffStack(os.path.join(output_path, "segmentation result.tif"));
		print(roi_imp);
		FileSaver(roi_imp).saveAsTiff(os.path.join(output_path, "vessel axis.tif"));

	egfp_imp.changes=False;
	mch_imp.changes=False;
	roi_imp.changes=False;
	fit_basis_imp.changes=False;
	pts_stack_imp.changes = False;
	egfp_imp.close();
	mch_imp.close();
	#roi_imp.close();
	fit_basis_imp.close();
	pts_stack_imp.close();
	
	zcoords = [i for i in range(composite_imp.getNSlices())];
	xyz_smooth_centres = [(x, y, z) for ((x, y), z) in zip(smooth_centres, zcoords)];
	
	composite_imp2 = straighten_vessel(composite_imp, xyz_smooth_centres, save_output=True);
	composite_imp3 = straighten_vessel(composite_imp2, xyz_smooth_centres, it=2, save_output=True);
	return composite_imp3;

#hsimp.addImageListener(UpdateRoiImageListener(rois));

def main():
	do_tubefitting();

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()