# @ImagePlus imp

from ij import ImageListener, IJ
from ij.gui import WaitForUserDialog, PointRoi, NonBlockingGenericDialog, OvalRoi
from ij.plugin import ChannelSplitter, RGBStackMerge
from java.awt import Point, Color

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
		
#		print(self.roi_list);
		
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
#	print("do linear interp between p[-2] = {} and p[-1] = {}".format(p[-2], p[-1]));
#	MyWaitForUser("pause", "pause");
	lininterp = linear_interp(p[-2], p[-1], n_points=n_points);
#	print("lininterp = {}".format(lininterp));
	c.extend(lininterp);
#	MyWaitForUser("pause", "pause");
	return c;

def resample_z(catmull_rom_points, z_points):
	"""
	Resample the interpolated list of points at given z-points
	EITHER remove interpolated points before first identified index to ensure that path doesn't jump around, 
	OR just do linear interpolation and smoothing...
	"""
	_,_,z = zip(*catmull_rom_points);
	unique_z = list(set(z));
	print("len(unique_z) = {}".format(len(unique_z)));
	unique_z_catmull_rom = [catmull_rom_point for catmull_rom_point in catmull_rom_points if catmull_rom_point[2] in unique_z];
	out_points = [];
	for ztarget in z_points:
		print("ztarget = {}".format(ztarget));
		loss_function = [abs(zz-ztarget) for zz in unique_z];
		print("loss_function = {}".format(loss_function));
		out_idx1 = loss_function.index(min(loss_function));
		print("out_idx1 = {}".format(out_idx1));
		print("unique_z_catmull_rom[out_idx1] = {}".format(unique_z_catmull_rom[out_idx1]));
		if out_idx1 == 0:
			out_coords = (unique_z_catmull_rom[out_idx1][0], unique_z_catmull_rom[out_idx1][1], ztarget);
			out_points.append(out_coords);
			continue;
		try:
			out_idx2 = loss_function.index(min([loss_function[out_idx1 - 1], loss_function[out_idx1 + 1]]));
			print("out_idx2 = {}".format(out_idx2));
			print("unique_z_catmull_rom[out_idx2] = {}".format(unique_z_catmull_rom[out_idx2]));
		except IndexError:
			print("INDEX ERROR! out_idx1 = {}, ztarget = {}".format(out_idx1, ztarget));
			out_coords = (unique_z_catmull_rom[out_idx1][0], unique_z_catmull_rom[out_idx1][1], ztarget);
			out_points.append(out_coords);
			continue;

		print("denominator = {}".format((unique_z_catmull_rom[out_idx2][2] - unique_z_catmull_rom[out_idx1][2])));
		deltaxy = [(unique_z_catmull_rom[out_idx2][idx] - unique_z_catmull_rom[out_idx1][idx]) * 
			(ztarget/(unique_z_catmull_rom[out_idx2][2] - unique_z_catmull_rom[out_idx1][2])) for idx in range(2)];
		out_coords = (unique_z_catmull_rom[out_idx1][0] + deltaxy[0], 
	    				unique_z_catmull_rom[out_idx1][1] + deltaxy[1], 
	    				ztarget);
		out_points.append(out_coords);
		unique_z_catmull_rom = unique_z_catmull_rom[out_idx1:];
		unique_z = unique_z[out_idx1:];
		print("trimmed len(unique_z) = {}".format(len(unique_z)));
		MyWaitForUser("pause", "pause");
	return out_points;

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

listener = UpdateRoiImageListener([]);
imp.addImageListener(listener);
IJ.setTool("elliptical");
step_size = 100;
next_z = 1;
while next_z < imp.getNSlices():
	imp.setZ(next_z)
	WaitForUserDialog("Add approx centre around every {}th frame, ensuring that the start (lowest z) of the region for unwrapping is included".format(step_size)).show();
	next_z = imp.getZ() + step_size
	listener.imageUpdated(imp);
#	imp.setZ(imp.getZ() + 1); # force update of roi list, possibly better to call imageUpdated directly?
	print(listener.getRoiList());
imp.setZ(imp.getNSlices());
WaitForUserDialog("Ensure that end (highest z) of the region for unwrapping has a manually defined point...").show();
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
vessel_ax_imp = IJ.createImage("Spline fitted central axis", "32-bit", imp.getWidth(), imp.getHeight(), imp.getNSlices());
for zidx, vessel_center in enumerate(output_positions):
	vessel_ax_imp.setZ(zidx+1);
	roi = OvalRoi(vessel_center[0], vessel_center[1], 5, 5);
	vessel_ax_imp.setRoi(roi);
	IJ.run(vessel_ax_imp, "Set...", "value=" + str(vessel_ax_imp.getProcessor().maxValue()) + " slice");
vessel_ax_imp.show();
channels  = ChannelSplitter().split(imp);
out_imp = RGBStackMerge().mergeChannels([channels[0], channels[1], vessel_ax_imp], True);
out_imp.setTitle("Manually defined vessel axis");
out_imp.show();
