# @ImagePlus imp
# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
import math

from ij import IJ, ImageStack, ImagePlus
from ij.gui import EllipseRoi, WaitForUserDialog
from ij.plugin import ChannelSplitter, Slicer, HyperStackConverter
from ij.process import StackProcessor, AutoThresholder, StackStatistics
from ij import WindowManager as WM

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

from UpdateImageRoiListener import UpdateRoiImageListener
import file_io as io
import ellipse_fitting


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
rot_imp.show();
disp_imp_ch1 = rot_imp.clone();
disp_imp_ch2 = Slicer().reslice(mch_imp);

# Apply 3d MEDIAN FILTER to denoise and emphasise vessel-associated voxels
filter_radius = 3.0;
IJ.run(rot_imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(filter_radius) + " z=" + str(filter_radius));
#IJ.run(rot_imp, "Median 3D...", "x=5 y=5 z=" + str(filter_radius));

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

# plane-wise, use binary-outline
# say the non-zero points then make up basis for fitting to be performed per http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
rois = [];
roi_stack = IJ.createImage("rois", fit_basis_imp.getWidth(), fit_basis_imp.getHeight(), fit_basis_imp.getNSlices(), 16);
for zidx in range(fit_basis_imp.getNSlices()):
	fit_basis_imp.setZ(zidx+1);
	
	IJ.run(fit_basis_imp, "Outline", "stack");
	IJ.run(fit_basis_imp, "Create Selection", "");
	roi = fit_basis_imp.getRoi();
	pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()]
	centre, angle, axl = ellipse_fitting.fit_ellipse(pts);
	print("Slice " + str(zidx+1) + 
		", centre = " + str(centre) + 
		", angle (deg) = " + str(angle * 180 / math.pi) + 
		", axes = " + str(axl));
	rot_imp.setZ(zidx+1);
	ellipse_roi = ellipse_fitting.generate_ellipse_roi(centre, angle, axl);
	rois.append(ellipse_roi);
	roi_stack.setZ(zidx+1);
	roi_stack.setRoi(ellipse_roi);
	IJ.run(roi_stack, "Draw", "slice");
	
rot_imp.changes = False;
rot_imp.close();
disp_imp_ch2.show()
disp_imp_ch1.show();
roi_stack.show();
IJ.run("Merge Channels...", "c1=[" + disp_imp_ch2.getTitle() + 
								"] c2=[" + disp_imp_ch1.getTitle() + 
								"] c7=[" + roi_stack.getTitle() + "] create keep");
								
#disp_imp_ch1.addImageListener(UpdateRoiImageListener(rois));


#hsimp.addImageListener(UpdateRoiImageListener(rois));
