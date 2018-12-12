# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math, sys, os

from ij import IJ, ImageStack, ImagePlus
from ij.gui import EllipseRoi, WaitForUserDialog
from ij.plugin import ChannelSplitter, Slicer, HyperStackConverter, ImageCalculator
from ij.process import StackProcessor, AutoThresholder, StackStatistics, FloatProcessor
from ij import WindowManager as WM
from loci.plugins import BF as bf

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

def generate_smoothed_vessel_axis(centres, pixel_size_um=0.1625):
	"""From a list of fitted vessel centres, generate the vessel path"""
	smooth_parameter_um = 1.0;
	out_centres = [];
	smooth_planes = smooth_parameter_um/pixel_size_um;
	# for each element, get elements +/- smooth_parameter_um/2 along the length of the line?
	# or, since this will add additional weight to points that are massively out of line, just do +/- smooth_parameter_um/um_per_z_plane
	#for idx, centre in enumerate(centres):
	#	sub_centres = [c for c in centres if (abs(centres.index(c) - idx) <= round(smooth_planes/2))];
	#	out_centres.append((sum([x for (x,y) in sub_centres])/len(sub_centres), sum([y for (x,y) in sub_centres])/len(sub_centres)));
	# or, avoid looping twice:
	xs = [x for (x,y) in centres];
	ys = [y for (x,y) in centres];
	for idx, centre in enumerate(centres):
		high_idx = idx + int(round(smooth_planes/2));
		high_idx = high_idx if high_idx < len(centres)-1 else len(centres)-2;
		low_idx = idx - int(round(smooth_planes/2));
		low_idx = low_idx if low_idx > 0 else 0;
		out_centres.append((sum(xs[low_idx:high_idx+1])/(high_idx+1-low_idx), sum(ys[low_idx:high_idx+1])/(high_idx+1-low_idx)));
	return out_centres;

def threshold_and_binarise(imp, z_xy_ratio):
	"""Return thresholded stack"""
	filter_radius = 3.0;
	IJ.run(imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(math.ceil(filter_radius / z_xy_ratio)) + " z=" + str(filter_radius));
	IJ.run(imp, "8-bit", "");
	rot_imp.show();

	# Apply automatic THRESHOLD to differentiate cells from background
	# get threshold value from stack histogram using otsu
	histo = StackStatistics(imp).histogram;
	#thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.IJ_IsoData, histo);
	thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histo);
	max_voxel_volume = int(float(imp.getHeight() * imp.getWidth() * imp.getNSlices())/100);
	IJ.run(rot_imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
												" min_size=" + str(max_voxel_volume) + " max_size=-1");
	fit_basis_imp = WM.getImage("Seg");
	IJ.setThreshold(fit_basis_imp, 1, 65535);
	IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");
	IJ.run(fit_basis_imp, "Fill Holes", "stack");
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

im_test_path = "D:\\data\\Marcksl1 cell shape analysis\\e27 ISV1.tif";
metadata_test_path = "D:\\data\Marcksl1 cell shape analysis\\Cropped\\2018-12-05 16-06-29 output\\Cropped UAS-marcksl1b-delED e27 xISV 1.json";
info = PrescreenInfo();
info.load_info_from_json(metadata_test_path);
z_xy_ratio = abs(info.get_z_plane_spacing_um()) / info.get_xy_pixel_size_um();
print(z_xy_ratio);
bfimp = bf.openImagePlus(im_test_path);
imp = bfimp[0];
imp.show();
IJ.run(imp, "Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

# split channels and get EC-label channel (GFP, ch1 - confirm this on a case-wise basis from acq metadata...)
channels  = ChannelSplitter().split(imp);
egfp_imp = channels[0];
mch_imp = channels[1];

# (crudely) ROTATE so that long axis is z axis...
cal = imp.getCalibration();
x_extent = cal.getX(egfp_imp.getWidth());
y_extent = cal.getY(egfp_imp.getHeight());
z_extent = cal.getZ(egfp_imp.getNSlices());
#orientation_str = "Top" if (y_extent > x_extent) else "Left";
#IJ.run(ml1_imp, "Reslice [/]...", "output=" + str(cal.getX(1.0)) + 
#								" start=" + orientation_str + " rotate avoid");
rot_imp = Slicer().reslice(egfp_imp);
#rot_imp.show();
disp_imp_ch1 = rot_imp.clone();
rot_imp_2 = Slicer().reslice(mch_imp);
disp_imp_ch2 = rot_imp_2.clone();
depth = rot_imp.getNSlices() if rot_imp.getNSlices() > rot_imp.getNFrames() else rot_imp.getNFrames();
width = rot_imp.getWidth();
height = int(round(rot_imp.getHeight() * z_xy_ratio));

# Apply 3d MEDIAN FILTER to denoise and emphasise vessel-associated voxels
fit_basis_imp = threshold_and_binarise(rot_imp, z_xy_ratio);

# plane-wise, use binary-outline
# say the non-zero points then make up basis for fitting to be performed per http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
rois = [];
centres = [];
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
		#print(pt);
		#print(int(round(fit_basis_imp.getHeight() * z_xy_ratio)));
		pix[int(pt[1]) * width + int(pt[0])] = 128;
	pts_stack.addSlice(ip);
	centre, angle, axl = ellipse_fitting.fit_ellipse(clean_pts);
	centres.append(centre);
	#print("Slice " + str(zidx+1) + 
	#	", centre = " + str(centre) + 
	#	", angle (deg) = " + str(angle * 180 / math.pi) + 
	#	", axes = " + str(axl));
	rot_imp.setZ(zidx+1);
	ellipse_roi = ellipse_fitting.generate_ellipse_roi(centre, angle, axl);
	rois.append(ellipse_roi);
	#roi_stack.setZ(zidx+1);
	#roi_stack.setRoi(ellipse_roi);
	#IJ.run(roi_stack, "Draw", "slice");
	#IJ.run(roi_stack, "Set...", "value=255 slice");
IJ.run(imp, "Line Width...", "line=1");

smooth_centres =  generate_smoothed_vessel_axis(centres, pixel_size_um=info.get_xy_pixel_size_um());
for zidx in range(fit_basis_imp.getNSlices()):
	centre = smooth_centres[zidx];
	ellipse_roi = EllipseRoi(centre[0]-2, centre[1], centre[0]+2, centre[1], 1.0);
	roi_stack.setZ(zidx+1);
	roi_stack.setRoi(ellipse_roi);
	IJ.run(roi_stack, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");

pts_stack_imp = ImagePlus("Cleaned points", pts_stack);
pts_stack_imp.show();
#WaitForUserDialog("Clean pts").show();
rot_imp.changes = False;
rot_imp.close();
disp_imp_ch2.show()
disp_imp_ch1.show();
roi_stack.show();
IJ.run(disp_imp_ch1, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
IJ.run(disp_imp_ch2, "Size...", "width=" + str(width) + " height=" + str(height) + " depth=" + str(depth) + " average interpolation=Bilinear");
print("disp_imp_ch2_size:"  + str((disp_imp_ch1.getHeight(), disp_imp_ch1.getWidth, disp_imp_ch1.getNSlices())))
IJ.run("Merge Channels...", "c1=[" + disp_imp_ch2.getTitle() + 
								"] c2=[" + disp_imp_ch1.getTitle() + 
								"] c7=[" + roi_stack.getTitle() + "] create keep");

WM.getImage("Composite").addImageListener(UpdateRoiImageListener(rois));
IJ.run(roi_stack, "8-bit", "");

disp_imp_ch1.changes=False;
disp_imp_ch2.changes=False;
imp.changes=False;
roi_stack.changes=False;
fit_basis_imp.changes=False;


#hsimp.addImageListener(UpdateRoiImageListener(rois));
