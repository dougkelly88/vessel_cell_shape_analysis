# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math, sys, os

from ij import IJ, ImageStack, ImagePlus
from ij.gui import EllipseRoi, WaitForUserDialog
from ij.plugin import ChannelSplitter, Slicer, HyperStackConverter
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

disp_imp_ch2 = Slicer().reslice(mch_imp);

# Apply 3d MEDIAN FILTER to denoise and emphasise vessel-associated voxels
filter_radius = 3.0;
IJ.run(rot_imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(math.ceil(filter_radius / z_xy_ratio)) + " z=" + str(filter_radius));
#IJ.run(rot_imp, "Median 3D...", "x=5 y=5 z=" + str(filter_radius));

rot_imp.show();

# Apply automatic THRESHOLD to differentiate cells from background
# get threshold value from stack histogram using otsu
IJ.run(rot_imp, "8-bit", "");
histo = StackStatistics(rot_imp).histogram;
#thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.IJ_IsoData, histo);
thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histo);
n_vox = rot_imp.getHeight() * rot_imp.getWidth() * rot_imp.getNSlices();
IJ.run(rot_imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
											" min_size=" + str(int(float(n_vox)/100)) + " max_size=-1");
fit_basis_imp = WM.getImage("Seg");
IJ.setThreshold(fit_basis_imp, 1, 65535);
IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");
IJ.run(fit_basis_imp, "Fill Holes", "stack");

# plane-wise, use binary-outline
# say the non-zero points then make up basis for fitting to be performed per http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
rois = [];
roi_stack = IJ.createImage("rois", rot_imp.getWidth(), int(round(rot_imp.getHeight() * z_xy_ratio)), rot_imp.getNSlices(), 16);
roi_stack = IJ.createImage("rois", rot_imp.getWidth(), rot_imp.getHeight(), rot_imp.getNSlices(), 16);
#pts_stack = ImageStack(fit_basis_imp.getWidth(), fit_basis_imp.getHeight());
for zidx in range(fit_basis_imp.getNSlices()):
	fit_basis_imp.setZ(zidx+1);
	IJ.run(fit_basis_imp, "Outline", "slice");
	IJ.run(fit_basis_imp, "Create Selection", "");
	roi = fit_basis_imp.getRoi();
	fit_basis_imp.killRoi();
	pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()];
	clean_pts = convex_hull_pts(pts);
	#clean_pts = [(x, z_xy_ratio * y) for (x, y) in clean_pts];
	# make a stack of clean points...
	#ip = FloatProcessor(fit_basis_imp.getWidth(), fit_basis_imp.getHeight())
	#pix = ip.getPixels();
	#for pt in clean_pts:
	#	pix[int(pt[1]) * fit_basis_imp.getWidth() + int(pt[0])] = 128;
	#pts_stack.addSlice(ip);
	centre, angle, axl = ellipse_fitting.fit_ellipse(clean_pts);
	#print("Slice " + str(zidx+1) + 
	#	", centre = " + str(centre) + 
	#	", angle (deg) = " + str(angle * 180 / math.pi) + 
	#	", axes = " + str(axl));
	rot_imp.setZ(zidx+1);
	ellipse_roi = ellipse_fitting.generate_ellipse_roi(centre, angle, axl);
	rois.append(ellipse_roi);
	roi_stack.setZ(zidx+1);
	roi_stack.setRoi(ellipse_roi);
	IJ.run(roi_stack, "Draw", "slice");

#pts_stack_imp = ImagePlus("Cleaned points", pts_stack);
#pts_stack_imp.show();
#WaitForUserDialog("Clean pts").show();
rot_imp.changes = False;
rot_imp.close();
disp_imp_ch2.show()
disp_imp_ch1.show();
roi_stack.show();
#IJ.run(disp_imp_ch1, "Size...", "width=" + str(rot_imp.getWidth()) + 
#								" height=" + str(int(round(rot_imp.getHeight() * z_xy_ratio))) + 
#								" depth=" + str(rot_imp.getNSlices()) + 
#								" average interpolation=Bilinear");
#IJ.run(disp_imp_ch2, "Size...", "width=" + str(rot_imp.getWidth()) + 
#								" height=" + str(int(round(rot_imp.getHeight() * z_xy_ratio))) + 
#								" depth=" + str(rot_imp.getNSlices()) + 
#								" average interpolation=Bilinear");
IJ.run("Merge Channels...", "c1=[" + disp_imp_ch2.getTitle() + 
								"] c2=[" + disp_imp_ch1.getTitle() + 
								"] c7=[" + roi_stack.getTitle() + "] create keep");
								
disp_imp_ch1.addImageListener(UpdateRoiImageListener(rois));


#hsimp.addImageListener(UpdateRoiImageListener(rois));
