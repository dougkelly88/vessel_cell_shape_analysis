# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math, sys, os

from ij import IJ, ImageStack, ImagePlus
from ij.gui import EllipseRoi, WaitForUserDialog, Roi
from ij.plugin import ChannelSplitter, Slicer, HyperStackConverter, ImageCalculator, Duplicator, SubstackMaker
from ij.process import StackProcessor, AutoThresholder, StackStatistics, FloatProcessor
from ij import WindowManager as WM
from loci.plugins import BF as bf

#im_test_path = "D:\\data\\Marcksl1 cell shape analysis\\e27 ISV1.tif";
#metadata_test_path = "D:\\data\Marcksl1 cell shape analysis\\Cropped\\2018-12-05 16-06-29 output\\Cropped UAS-marcksl1b-delED e27 xISV 1.json";
#im_test_path = "D:\\data\\Structural imaging\\Imaging protocol tests\\2018-09-14 vessel structure imaging\\vessel_structure_protocol3_20180914_174448\\Cropped vessels\\2018-12-12 16-38-31 output\\Cropped dextran-rhodamine kdrl-egfp vessel_structure etest_embryo aISV 1.tif";
#metadata_test_path = "D:\\data\\Structural imaging\\Imaging protocol tests\\2018-09-14 vessel structure imaging\\vessel_structure_protocol3_20180914_174448\\Cropped vessels\\2018-12-12 16-38-31 output\\Cropped dextran-rhodamine kdrl-egfp vessel_structure etest_embryo aISV 1.json";
im_test_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.tif";
metadata_test_path = "D:\\data\Marcksl1 cell shape analysis\\2018-12-31 Lumen stained samples\\AB inj marcksl1b-EGFP, rhodamine-dextran uangiography\\Cropped\\2019-01-10 16-03-47 output\\Cropped Dextran-Rhodamine Marcksl1b-EGFP Lumen staining test eE1 xISV 1.json";


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

def generate_smoothed_vessel_axis(centres, pixel_size_um=0.1625):
	"""From a list of fitted vessel centres, generate the vessel path"""
	smooth_parameter_um = 5.0;
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

def threshold_and_binarise(imp, z_xy_ratio):
	"""Return thresholded stack"""
	print("performing segmentation on channel: " + imp.getTitle());
	filter_radius = 3.0;
	IJ.run(imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(math.ceil(filter_radius / z_xy_ratio)) + " z=" + str(filter_radius));
	IJ.run(imp, "8-bit", "");
	imp.show();

	# Apply automatic THRESHOLD to differentiate cells from background
	# get threshold value from stack histogram using otsu
	histo = StackStatistics(imp).histogram;
	#thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.IJ_IsoData, histo);
	thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histo);
	max_voxel_volume = int(float(imp.getHeight() * imp.getWidth() * imp.getNSlices())/100);
	IJ.run(imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
												" min_size=" + str(max_voxel_volume) + " max_size=-1");
	fit_basis_imp = WM.getImage("Seg");
	IJ.setThreshold(fit_basis_imp, 1, 65535);
	IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");
	#IJ.run(fit_basis_imp, "Fill Holes", "stack");
	IJ.run("3D Fill Holes", "");
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
	#print(pts);
	clean_pts = [];
	temp_pts = [];
	ys = [y for (x,y) in pts]
	for yy in set(ys):
	    xs = sorted([x for (x,y) in pts if y==yy]);
	    temp_pts.append((xs[0], yy));
	    if len(xs)>1:
	        temp_pts.append((xs[-1], yy));
	#print(temp_pts)
	xs = [x for (x,y) in temp_pts];
	for xx in set(xs):
	    ys = sorted([y for (x,y) in temp_pts if x==xx]);
	    clean_pts.append((xx, ys[0]));
	    if len(xs)>1:
	        clean_pts.append((xx, ys[-1]));
	#print(clean_pts);
	return clean_pts;

def split_and_rotate(imp, info):
	"""return image to segment on, image to project out, and images to display"""
	# for now, assume that these are ISVs and that embryo is mounted in familiar fashion. First of these can be developed out...
	seg_ch_idx, proj_ch_idx = ui.choose_segmentation_and_projection_channels(info);
	channels  = ChannelSplitter().split(imp);
	seg_imp = Duplicator().run(channels[seg_ch_idx]); # use Duplicator to decouple - can do smarter to save memory?
	proj_imp = Duplicator().run(channels[proj_ch_idx]);
	rot_seg_imp = Slicer().reslice(seg_imp);
	rot_seg_imp.setTitle("rot_seg_imp");
	rot_proj_imp = Slicer().reslice(proj_imp);
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
			egfp_mch_imps.append(Slicer().reslice(Duplicator().run(channels[ch_idx])));
	imp.changes=False;
	imp.close();
	seg_imp.changes = False;
	proj_imp.changes = False;
	seg_imp.close();
	proj_imp.close();
	return rot_seg_imp, rot_proj_imp, egfp_mch_imps

def main():
	info = PrescreenInfo();
	info.load_info_from_json(metadata_test_path);
	z_xy_ratio = abs(info.get_z_plane_spacing_um()) / info.get_xy_pixel_size_um();
	print(z_xy_ratio);
	bfimp = bf.openImagePlus(im_test_path);
	imp = bfimp[0];
	imp.show();
	cal = imp.getCalibration();
	IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

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
	roi_stack = IJ.createImage("rois", width, height, depth, 16);
	#roi_stack = IJ.createImage("rois", rot_imp.getWidth(), rot_imp.getHeight(), rot_imp.getNSlices(), 16);
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

	smooth_centres, tangent_vecs =  generate_smoothed_vessel_axis(centres, pixel_size_um=info.get_xy_pixel_size_um());
	for zidx in range(fit_basis_imp.getNSlices()):
		centre = smooth_centres[zidx];
		major_axis = major_axes[zidx];
		ellipse_roi = EllipseRoi(centre[0]-2, centre[1], centre[0]+2, centre[1], 1.0);
		roi_stack.setZ(zidx+1);
		roi_stack.setRoi(ellipse_roi);
		IJ.run(roi_stack, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");

	pts_stack_imp = ImagePlus("Cleaned points", pts_stack);
	pts_stack_imp.setTitle("pts_stack_imp");
	pts_stack_imp.show();

	rot_seg_imp.changes = False;
	rot_seg_imp.close();
	egfp_imp = egfp_mch_imps[0];
	mch_imp = egfp_mch_imps[1];
	egfp_imp.show()
	mch_imp.show();
	roi_stack.show();
	print("box height um = " + str(roi_stack.getNSlices() * info.get_xy_pixel_size_um()));
	IJ.run(egfp_imp, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
	IJ.run(mch_imp, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
	IJ.run("Merge Channels...", "c1=[" + mch_imp.getTitle() + 
									"] c2=[" + egfp_imp.getTitle() + 
									"] c7=[" + roi_stack.getTitle() + "] create keep");

	WM.getImage("Composite").addImageListener(UpdateRoiImageListener(rois));
	IJ.run(roi_stack, "8-bit", "");

	egfp_imp.changes=False;
	mch_imp.changes=False;
	roi_stack.changes=False;
	fit_basis_imp.changes=False;


#hsimp.addImageListener(UpdateRoiImageListener(rois));


# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()