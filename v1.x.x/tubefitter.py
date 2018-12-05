# @ImagePlus imp
# use D:\data\Marcksl1 cell shape analysis\e27 ISV1.tif for testing,,,
from Jama import Matrix, EigenvalueDecomposition
import math

from ij import IJ, ImageStack, ImagePlus
from ij.gui import EllipseRoi, WaitForUserDialog
from ij.plugin import ChannelSplitter, Slicer
from ij.process import StackProcessor, AutoThresholder, StackStatistics
from ij import WindowManager as WM
from ij import ImageListener

class UpdateRoiImageListener(ImageListener):
	"""class to support updating ROI from list upon change of frame"""
	def __init__(self, roi_list):
		self.last_slice = 1;
		self.roi_list = roi_list;
		print("UpdateRoiImageListener started");

	def imageUpdated(self, imp):
		slc = imp.getZ();
		print(slc);
		self.last_slice = slc;
		if slc < len(self.roi_list):
			imp.setRoi(self.roi_list[slc - 1]);

	def imageOpened(self, imp):
		print("UpdateRoiImageListener: image opened");
			
	def imageClosed(self, imp):
		print("UpdateRoiImageListener: image closed");
		imp.removeImageListener(self);

	def getRoiList(self):
		return self.roi_list;

def ellipse_center(a):
	"""return centre of least squares best fit ellipse"""
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	num = b*b-a*c
	x0=(c*d-b*f)/num
	y0=(a*f-b*d)/num
	return (x0,y0)

def ellipse_angle_of_rotation(a):
	"""return rotation angle of least squares best fit ellipse"""
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	return 0.5*math.atan(2*b/(a-c))

def ellipse_axis_length(a):
	"""return axis length of least squares best fit ellipse"""
	b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
	up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
	down1=(b*b-a*c)*( (c-a)*math.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	down2=(b*b-a*c)*( (a-c)*math.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
	res1=math.sqrt(abs(up/down1));
	res2=math.sqrt(abs(up/down2));
	return (res1, res2)

def fit_ellipse(coords):
	"""find a best fit ellipse (least squares sense) from a list of co-ordinates"""
	# based on http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
	xs = [float(x) for (x, y) in pts];
	xmean = sum(xs)/len(xs);
	xs = [x - xmean for x in xs];
	ys = [float(y) for (x, y) in pts];
	ymean = sum(ys)/len(ys);
	ys = [y - ymean for y in ys];

	N = 6;
	d1 = [x * x for x in xs];
	d2 = [x * y for (x, y) in zip(xs, ys)];
	d3 = [y * y for y in ys];
	d4 = [x for x in xs];
	d5 = [y for y in ys];
	d6 = [1 for x in xs];
	d = sum([d1, d2, d3, d4, d5, d6], []);
	D = Matrix(d, len(d1));
	
	S = D.transpose().times(D);
	C = Matrix(6, 6);
	C.set(0, 2, 2);
	C.set(2, 0, 2);
	C.set(1, 1, -1);
	eig = EigenvalueDecomposition(S.inverse().times(C.transpose()));
	E = eig.getRealEigenvalues();
	#print("eigenvalues = " + str([('%.2E' % e) for e in E]));
	V = eig.getV();
	#print("eigenvectors = ");
	#for idx in range(0,N):
	#	print(str([('%.2E' % v) for v in vs[N * idx : (idx + 1) * N]]));
	absE = [abs(e) for e in E];
	n = absE.index(max(absE));
	a = [V.get(idx, n) for idx in range(0, N)];
	
	(xc, yc) = ellipse_center(a);
	xc += xmean;
	yc += ymean;
	centre = (xc, yc);
	angle = ellipse_angle_of_rotation(a);
	axes_l = ellipse_axis_length(a);
	return centre, angle, axes_l;

def generate_ellipse_roi(centre, angle, axes_l):
	"""translate the outputs of ellipse fitting to an ellipse ROI"""
	xc = centre[0]; yc = centre[1];
	majAxL = axes_l[0]; minAxL = axes_l[1];
	x1 = xc - majAxL * math.cos(angle); x2 = xc + majAxL * math.cos(angle);
	y1 = yc - majAxL * math.sin(angle); y2 = yc + majAxL * math.sin(angle);
	aspectRatio = float(minAxL)/float(majAxL);
	roi = EllipseRoi(x1, y1, x2, y2, aspectRatio);
	return roi;

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
disp_title = rot_imp.getTitle();

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
#WaitForUserDialog("Move reslice image to foreground...").show();
for zidx in range(fit_basis_imp.getNSlices()):
#zidx = 225;
	fit_basis_imp.setZ(zidx);
	
	IJ.run(fit_basis_imp, "Outline", "stack");
	IJ.run(fit_basis_imp, "Create Selection", "");
	#fit_basis_imp.show();
	roi = fit_basis_imp.getRoi();
	pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()]
	centre, angle, axl = fit_ellipse(pts);
	print("Slice " + str(zidx) + 
		", centre = " + str(centre) + 
		", angle (deg) = " + str(angle * 180 / math.pi) + 
		", axes = " + str(axl));
	rot_imp.setZ(zidx);
	ellipse_roi = generate_ellipse_roi(centre, angle, axl);
	rot_imp.setRoi(ellipse_roi);
	rois.append(ellipse_roi);
disp_imp = WM.getImage(disp_title);
disp_window = WM.getWindow(disp_title);
WM.setCurrentWindow(disp_window);
disp_imp.addImageListener(UpdateRoiImageListener(rois));
#pts = [(x, y) for x, y in zip(roi.getContainedPoints().xpoints, roi.getContainedPoints().ypoints)];
#print(pts);
