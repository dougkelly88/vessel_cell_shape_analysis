from ij import IJ
from ij.gui import EllipseRoi
import math
from Jama import Matrix, EigenvalueDecomposition


c = (50, 50);
angle = 0;
axl = (30, 10);

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
	res1=math.sqrt(up/down1)
	res2=math.sqrt(up/down2)
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
	print(aspectRatio);
	roi = EllipseRoi(x1, y1, x2, y2, aspectRatio);
	return roi;

imp = IJ.createImage("test", 100, 100, 1, 8);
imp.show()
roi = generate_ellipse_roi(c, angle, axl);
imp.setRoi(roi);
print(roi);
IJ.run(imp, "Set...", "value=128");

# try fitting...
IJ.setThreshold(imp, 1, 255);
IJ.run(imp, "Convert to Mask", "method=Default background=Dark list");
IJ.run(imp, "Dilate", "");
IJ.run(imp, "Outline", "stack");
IJ.run(imp, "Create Selection", "");
roi = imp.getRoi();
pts = [(pt.x, pt.y) for pt in roi.getContainedPoints()]
centre, angle, axl = fit_ellipse(pts);
print("centre = " + str(centre));
print("angle = " + str(angle * 180 / math.pi));
print("axes = " + str(axl));



