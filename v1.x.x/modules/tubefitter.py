# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math, sys, os, csv

from ij import IJ, ImageStack, ImagePlus, Prefs, ImageListener
from ij.gui import EllipseRoi, WaitForUserDialog, Roi, PolygonRoi, NonBlockingGenericDialog, OvalRoi, PointRoi
from ij.io import FileSaver
from ij.plugin import ChannelSplitter, ImageCalculator, Duplicator, SubstackMaker, Straightener, RGBStackMerge
from ij.process import StackProcessor, AutoThresholder, StackStatistics, FloatProcessor
from ij import WindowManager as WM
from loci.plugins import BF as bf
from java.awt import Point, Color


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

from PrescreenInfo import PrescreenInfo
import file_io as io
import ellipse_fitting
import ui
import utils

def MyWaitForUser(title, message):
	"""non-modal dialog with option to cancel the analysis"""
	dialog = NonBlockingGenericDialog(title);
	dialog.setCancelLabel("Cancel analysis");
	if type(message) is list:
		for line in message:
			dialog.addMessage(line);
	else:
		dialog.addMessage(message);
	dialog.showDialog();
	if dialog.wasCanceled():
		raise KeyboardInterrupt("Run canceled");
	return;

class ManualSegmentationImageListener(ImageListener):
	"""class to support updating ROI from list upon change of frame"""
	def __init__(self, roi_list):
		self.last_slice = 1;
		self.roi_list = roi_list;
		print("ManualSegmentationImageListener started");

	def imageUpdated(self, imp):
		roi = imp.getRoi();
		if imp.getNSlices() > imp.getNFrames():
			slc = imp.getZ();
		else:
			slc = imp.getT();
		if roi is not None:
			if not roi.isArea():
				pt = roi.getContainedPoints()[0];
			else:
				rot_center = roi.getRotationCenter(); 
				x = int(rot_center.xpoints[0]);
				y = int(rot_center.ypoints[0]);
				print("x = {}".format(x));
				print("y = {}".format(y));
				pt = Point(x, y);
			print("adding {} to list at slice {}".format(pt, self.last_slice));
			if self.last_slice in [z for (x, y, z) in self.roi_list]:
				self.roi_list[[z for (x, y, z) in self.roi_list].index(self.last_slice)] = (pt.x, pt.y, slc);
				print(self.roi_list);
			else:
				self.roi_list.append((pt.x, pt.y, self.last_slice));
			imp.killRoi();
		
#		if slc in [z for (x, y, z) in self.roi_list]:
#			roi_list_entry = self.roi_list[[z for (x, y, z) in self.roi_list].index(slc)];
#			print("found entry in list ({}), displaying...".format(roi_list_entry));
#			roi = PointRoi(roi_list_entry[0], roi_list_entry[1]);
#			imp.setRoi(roi);
		self.last_slice = slc;
		
#		print(self.roi_list);
		
	def imageOpened(self, imp):
		print("ManualSegmentationImageListener: image opened");
			
	def imageClosed(self, imp):
		print("ManualSegmentationImageListener: image closed");
		imp.removeImageListener(self);

	def getRoiList(self):
		return self.roi_list;

def catmull_rom_spline(P0, P1, P2, P3, n_points=100):
	"""
	P0, P1, P2, and P3 should be (x,y,z)-tuples that are knots defining the Catmull-Rom spline.
	nPoints is the number of points to include in this curve segment.
	Modified from https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline#Code_example_in_Python to be 
	independent of Numpy
	"""
	# Calculate t0 to t4
	alpha = 0.5 # 0: uniform sampling; 0.5: centripetal sampling; 1: chordal sampling
	def tj(ti, Pi, Pj):
		xi, yi, zi = Pi
		xj, yj, zj = Pj
		return ( ( (xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)**0.5 )**alpha + ti

	t0 = 0
	t1 = tj(t0, P0, P1)
	t2 = tj(t1, P1, P2)
	t3 = tj(t2, P2, P3)

	# Only calculate points between P1 and P2
	t = [t1+(t2-t1)*float(idx)/(n_points-1) for idx in range(n_points)];

	# Implement Catmull-Rom in pure python :/
	A1 = [[(t1-tt)/(t1-t0)*P0[idx] + (tt-t0)/(t1-t0)*P1[idx] for idx in range(3)] for tt in t];
	A2 = [[(t2-tt)/(t2-t1)*P1[idx] + (tt-t1)/(t2-t1)*P2[idx] for idx in range(3)] for tt in t];
	A3 = [[(t3-tt)/(t3-t2)*P2[idx] + (tt-t2)/(t3-t2)*P3[idx] for idx in range(3)] for tt in t];
	B1 = [[(t2-tt)/(t2-t0)*A1[tidx][idx] + (tt-t0)/(t2-t0)*A2[tidx][idx] for idx in range(3)] for tidx, tt in enumerate(t)];
	B2 = [[(t3-tt)/(t3-t1)*A2[tidx][idx] + (tt-t1)/(t3-t1)*A3[tidx][idx] for idx in range(3)] for tidx, tt in enumerate(t)];
	C  = [[(t2-tt)/(t2-t1)*B1[tidx][idx] + (tt-t1)/(t2-t1)*B2[tidx][idx] for idx in range(3)] for tidx, tt in enumerate(t)];
	return C

def linear_interp(p0, p1, n_points=100):
	"""
	Return linear interpolation for the segment between knots P0 and P1
	"""
	t = [float(idx)/(n_points - 1) for idx in range(n_points)];
	c = [[p0[idx] + tt*(p1[idx] - p0[idx]) for idx in range(3)] for tt in t];
	return c;

def catmull_rom_chain(p, n_points=100):
	"""
	Calculate Catmull Rom for a list of (x, y, z)-tuples and return the combined curve, adding linear segments at start
	and end. 
	"""
	sz = len(p)

	# The curve C will contain an array of (x,y) points.
	c = []
	c.extend(linear_interp(p[0], p[1], n_points=n_points));
	for i in range(sz-3):
		cc = catmull_rom_spline(p[i], p[i+1], p[i+2], p[i+3], n_points=n_points)
		c.extend(cc)
	lininterp = linear_interp(p[-2], p[-1], n_points=n_points);
	c.extend(lininterp);
	return c;

def new_resample_z(spline_points, target_zs):
	spline_zs = [z for _,_,z in spline_points];
	out_points = [];
	closest_idxs = [];
	for idx, z_target in enumerate(target_zs):
	    lf = [abs(spline_z - z_target) for spline_z in spline_zs];
	    # this might fail when spline has regions of decreasing z!
	    closest_idx = lf.index(min(lf));
	    print("Target z = {}, achieved point = {}".format(z_target, spline_points[closest_idx]));
	    closest_idxs.append(closest_idx);
	    out_points.append((spline_points[closest_idx][0], spline_points[closest_idx][1], z_target));
	# DO SMOOTHING OF OUTPUT POINTS?
	return out_points, closest_idxs;

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
	bin_imp = WM.getImage("Bin");
	bin_imp.changes=False;
	bin_imp.close();
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
	
def straighten_vessel(imp, smooth_centres, it=1):
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
	FileSaver(new_composite).saveAsTiffStack(os.path.join(output_path, "after rotation " + str(it) + ".tif"));
	return new_composite;

def manual_axis_definition(imp, output_path=None):
	"""perform manual steps to determine vessel centreline"""
	listener = ManualSegmentationImageListener([]);
	imp.addImageListener(listener);
	IJ.setTool("elliptical");
	step_size = 100;
	next_z = 1;
	while next_z < imp.getNSlices():
		imp.setZ(next_z)
		MyWaitForUser("Manual segmentation...", "Add approx centre around every {}th frame by inscribing an ellipse, ensuring that the start (lowest z) of the region for unwrapping is included".format(step_size));
		next_z = imp.getZ() + step_size
		listener.imageUpdated(imp);
		#print(listener.getRoiList());
	imp.setZ(imp.getNSlices());
	MyWaitForUser("Manual segmentation...", "Ensure that end (highest z) of the region for unwrapping has a manually defined point...");
	listener.imageUpdated(imp);
	manual_centers = listener.getRoiList();
	imp.removeImageListener(listener);

	#manual_centers = [(270, 135, 1), (258, 102, 101), (286, 63, 201), (299, 60, 301), (311, 77, 401), (335, 131, 501), (364, 156, 601), (382, 132, 701), (402, 120, 801), (382, 142, 901), (322, 200, 1001)];
	manual_centers = sorted(manual_centers, key=lambda xyz: xyz[2]);
	spline_interp_centers = catmull_rom_chain(manual_centers, n_points=5*step_size);
	target_zs = [z for z in range(manual_centers[0][2], manual_centers[-1][2]+1)];
	output_positions, _ = new_resample_z(spline_interp_centers, target_zs);
	print(output_positions);
	print("len(output_positions) = {}".format(len(output_positions)));
	vessel_ax_imp = IJ.createImage("Spline fitted central axis", "16-bit", imp.getWidth(), imp.getHeight(), imp.getNSlices());
	for zidx, vessel_center in enumerate(output_positions):
		vessel_ax_imp.setZ(zidx+1);
		roi = OvalRoi(vessel_center[0], vessel_center[1], 5, 5);
		vessel_ax_imp.setRoi(roi);
		IJ.run(vessel_ax_imp, "Set...", "value=" + str(vessel_ax_imp.getProcessor().maxValue()) + " slice");
	vessel_ax_imp.show();
	channels  = ChannelSplitter().split(imp);
	out_imp = RGBStackMerge().mergeChannels([channels[0], channels[1], vessel_ax_imp], True);
	out_imp.setTitle("Manually defined vessel axis");
	vessel_ax_imp.changes = False;
	vessel_ax_imp.close();
	out_imp.show();
	try:
		FileSaver(out_imp).saveAsTiffStack(os.path.join(output_path, "manually determined vessel axis.tif"));
	except IOError as e:
		print(e.message);
	save_centers(output_positions, output_path);
	return output_positions, out_imp;

def save_centers(centers, output_folder):
	"""save vessel axis to csv for later loading"""
	out_path = os.path.join(output_folder, "vessel_axis.csv");
	f = open(out_path, 'wb');
	try:
		writer = csv.writer(f);
		for c in centers:
			writer.writerow([c[0], c[1], c[2]]);
	except IOError as e:
		print("problem saving, {}".format(e));
	finally:
		f.close();

def load_centers(centers_path):
	"""load vessel axis from previously determined data"""
	f = open(centers_path, 'rb');
	centers = [];
	try:
		csvreader = csv.reader(f);
		for row in csvreader:
			centers.append((float(row[0]), float(row[1]), float(row[2])));
	except IOError as e:
		print("problem saving, {}".format(e));
	finally:
		f.close();
	return centers;

def do_tubefitting(im_path=im_test_path, metadata_path=metadata_test_path, output_path=output_path, load_centers_path=None):
	# todo: fix things so that all operations use a consistent definition of background rather than changing Prefs on the fly...
	Prefs.blackBackground = False;
	info = PrescreenInfo();
	info.load_info_from_json(metadata_path);
	z_xy_ratio = abs(info.get_z_plane_spacing_um()) / info.get_xy_pixel_size_um();
	bfimp = bf.openImagePlus(im_path);
	imp = bfimp[0];
	imp.show();
	cal = imp.getCalibration();
	IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

	rot_seg_imp, rot_proj_imp, egfp_mch_imps = split_and_rotate(imp, info);
	if "isv" in info.get_vessel_type().lower():
		depth = rot_seg_imp.getNSlices() if rot_seg_imp.getNSlices() > rot_seg_imp.getNFrames() else rot_seg_imp.getNFrames();
		width = rot_seg_imp.getWidth();
		height = int(round(rot_seg_imp.getHeight() * z_xy_ratio));
	else:
		depth = rot_seg_imp.getNSlices() if rot_seg_imp.getNSlices() > rot_seg_imp.getNFrames() else rot_seg_imp.getNFrames();
		width = int(round(rot_seg_imp.getWidth() * z_xy_ratio));
		height = rot_seg_imp.getHeight();
	print(z_xy_ratio);

	merged_rotated_imp = RGBStackMerge().mergeChannels([rot_seg_imp, rot_proj_imp], False);
	#merged_rotated_imp.show();
	IJ.run(merged_rotated_imp, 
			"Size...", 
			"width={} height={} depth={} average interpolation=Bilinear".format(width, height, depth));
	merged_rotated_imp = utils.convert_multichannel_stack_to_Xbit(merged_rotated_imp, bitdepth=16);
	merged_rotated_imp.setTitle("converted merged_rotated_imp");
	merged_rotated_imp.show();
	#MyWaitForUser("pause", "pause");
	FileSaver(merged_rotated_imp).saveAsTiffStack(os.path.join(output_path, "rotated to align z to DV axis.tif"));
	merge_rot_imp_with_axis = None;
	if load_centers_path is None:
		xyz_smooth_centres, merge_rot_imp_with_axis = manual_axis_definition(merged_rotated_imp, output_path=output_path);
		utils.convert_multichannel_stack_to_Xbit(merge_rot_imp_with_axis, bitdepth=16);
		merged_rotated_imp.close();
	else:
		xyz_smooth_centres = load_centers(load_centers_path);
	
	if merge_rot_imp_with_axis is not None:
		composite_imp2 = straighten_vessel(merge_rot_imp_with_axis, xyz_smooth_centres);
		merge_rot_imp_with_axis.changes = False;
		merge_rot_imp_with_axis.close();
	else:
		composite_imp2 = straighten_vessel(merged_rotated_imp, xyz_smooth_centres);
		merged_rotated_imp.changes = False;
		merged_rotated_imp.close();
	
	#composite_imp2.show();
	#MyWaitForUser("pause", "after straightening in first axis");
	#composite_imp2.hide();
	composite_imp3 = straighten_vessel(composite_imp2, xyz_smooth_centres, it=2);
	composite_imp2.close();
	#composite_imp3.show();
	#MyWaitForUser("pause", "after straightening in second axis");
	#composite_imp3.hide();
	#utils.convert_multichannel_stack_to_16bit(composite_imp3);
	return composite_imp3, info;

#hsimp.addImageListener(UpdateRoiImageListener(rois));

def main():
	do_tubefitting();

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()