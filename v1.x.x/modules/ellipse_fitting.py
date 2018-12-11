# module containing functions to handle fitting ellipses (least squares sense)
#
# D. J. Kelly, 2018-12-05, douglas.kelly@riken.jp

import math
from ij import IJ
from ij.gui import EllipseRoi
from Jama import Matrix, EigenvalueDecomposition

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
	res1=math.sqrt(abs(up/down1))
	res2=math.sqrt(abs(up/down2))
	return res1, res2;

def fit_ellipse(coords):
	"""find a best fit ellipse (least squares sense) from a list of co-ordinates"""
	# based on http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
	xs = [float(x) for (x, y) in coords];
	xmean = sum(xs)/len(xs);
	xs = [x - xmean for x in xs];
	ys = [float(y) for (x, y) in coords];
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
	majAxL = axes_l[1]; minAxL = axes_l[0];
	if axes_l[0] > axes_l[1]:
		aspectRatio = float(axes_l[1])/float(axes_l[0]);
		x1 = xc - axes_l[0] * math.cos(angle); x2 = xc + axes_l[0] * math.cos(angle);
		y1 = yc - axes_l[0] * math.sin(angle); y2 = yc + axes_l[0] * math.sin(angle);
	else:
		aspectRatio = float(axes_l[0])/float(axes_l[1]);
		angle = angle + math.pi/2;
		print("angle = " + str(180*angle/math.pi));
		x1 = xc - axes_l[1] * math.cos(angle); x2 = xc + axes_l[1] * math.cos(angle);
		y1 = yc - axes_l[1] * math.sin(angle); y2 = yc + axes_l[1] * math.sin(angle);
	print("Aspect ratio = " + str(aspectRatio));
	roi = EllipseRoi(x1, y1, x2, y2, aspectRatio);
	return roi;