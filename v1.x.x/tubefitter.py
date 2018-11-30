# @ImagePlus imp
# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
from ij import IJ, ImageStack, ImagePlus
from ij.plugin import ChannelSplitter, Slicer
from ij.process import StackProcessor, AutoThresholder, StackStatistics
from ij import WindowManager as WM

# split channels and get EC-label channel (GFP, ch1 - confirm this on a case-wise basis from acq metadata...)
channels  = ChannelSplitter().split(imp);
egfp_imp = channels[0];

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

# Apply 3d MEDIAN FILTER to denoise and emphasise vessel-associated voxels
filter_radius = 5.0;
#IJ.run(rot_imp, "Median 3D...", "x=" + str(filter_radius) + " y=" + str(filter_radius) + " z=" + str(filter_radius));
IJ.run(rot_imp, "Median 3D...", "x=1 y=1 z=" + str(filter_radius));

# Apply automatic THRESHOLD to differentiate cells from background
# get threshold value from stack histogram using otsu
IJ.run(rot_imp, "8-bit", "");
histo = StackStatistics(rot_imp).histogram;
thresh_lev = AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histo);
n_vox = rot_imp.getHeight() * rot_imp.getWidth() * rot_imp.getNSlices();
IJ.run(rot_imp, "3D Simple Segmentation", "low_threshold=" + str(thresh_lev + 1) + 
											" min_size=" + str(int(float(n_vox)/100)) + " max_size=-1");
fit_basis_imp = WM.getImage("Seg");
IJ.setThreshold(fit_basis_imp, 1, 65535);
IJ.run(fit_basis_imp, "Convert to Mask", "method=Default background=Dark list");

# plane-wise, use binary-outline
# say the non-zero points then make up basis for fitting to be performed per http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
fit_basis_imp.setZ(262);
IJ.run(fit_basis_imp, "Outline", "stack");
IJ.run(fit_basis_imp, "Create Selection", "");
fit_basis_imp.show();
roi = fit_basis_imp.getRoi();
pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()]
print(pts)
#pts = [(x, y) for x, y in zip(roi.getContainedPoints().xpoints, roi.getContainedPoints().ypoints)];
#print(pts);
