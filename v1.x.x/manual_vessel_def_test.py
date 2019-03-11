# @ImagePlus imp

from ij import ImageListener, IJ
from ij.gui import WaitForUserDialog, PointRoi
from java.awt import Point

class UpdateRoiImageListener(ImageListener):
	"""class to support updating ROI from list upon change of frame"""
	def __init__(self, roi_list):
		self.last_slice = 1;
		self.roi_list = roi_list;
		print("UpdateRoiImageListener started");

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
		
		print(self.roi_list);
		
	def imageOpened(self, imp):
		print("UpdateRoiImageListener: image opened");
			
	def imageClosed(self, imp):
		print("UpdateRoiImageListener: image closed");
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

    # Reshape so that we can multiply by the points P0 to P3
    # and get a point for each value of t.
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

def catmull_rom_chain(p):
    """
    Calculate Catmull Rom for a list of (x, y, z)-tuples and return the combined curve, adding linear segments at start
    and end. 
    """
    sz = len(p)

    # The curve C will contain an array of (x,y) points.
    c = []
    c.extend(linear_interp(p[0], p[1]));
    for i in range(sz-3):
        cc = catmull_rom_spline(p[i], p[i+1], p[i+2], p[i+3])
        c.extend(cc)
    c.extend(linear_interp(p[-2], p[-1]));
    return c;

def resample_z(catmull_rom_points, z_points):
    """
    Resample the interpolated list of points at given z-points
    """
    _,_,z = zip(*catmull_rom_points);
    out_idx = [[abs(zz-ztarget) for zz in z].index(min([abs(zz-ztarget) for zz in z])) for ztarget in z_points];
    return [catmull_rom_points[idx] for idx in out_idx];


listener = UpdateRoiImageListener([]);
imp.addImageListener(listener);
IJ.setTool("elliptical");
step_size = 100;
while next_z < imp.getNSlices():
	imp.setZ(next_z)
	WaitForUserDialog("Add approx centre around every {}th frame".format(step_size)).show();
	next_z = imp.getZ() + imp.getNSlices()//step_size
manual_centers = listener.getRoiList();
spline_interp_centers = catmull_rom_chain(manual_centers);
output_positions = resample_z(spline_interp_centers, [z+1 for z in range(imp.getNSlices())]);
print(output_positions);
print("len(output_positions) = {}".format(output_positions));