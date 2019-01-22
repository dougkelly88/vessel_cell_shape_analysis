# @ImagePlus imp
import math, sys, os
from ij import ImageStack, ImagePlus, IJ, Prefs
from ij.gui import PolygonRoi, Roi, WaitForUserDialog, Line, ProfilePlot
from ij.plugin import Duplicator, MontageMaker, ChannelSplitter, ImageCalculator
#from ij.plugin.filter import MaximumFinder
from ij.process import FloatProcessor, AutoThresholder

def twist_and_unwrap(imp):
	"""from the output of angular projection, define an axis along which to unzip the vessel and return this unzipped image"""
	imp.show();
	IJ.setTool("freeline");
	WaitForUserDialog("Input required...", "Please draw the line down which unwrapping should occur...").show();
	roi = imp.getRoi();
	if roi is None:
		raise ValueError;
	unwrap_poly = roi.getPolygon();
	unwrap_poly_xs = [x for x in unwrap_poly.xpoints];
	unwrap_poly_ys = [y for y in unwrap_poly.ypoints];
	# extend to the top and bottom of the image:
	unwrap_poly_xs.insert(0, unwrap_poly_xs[0]);
	unwrap_poly_xs.append(unwrap_poly_xs[-1])
	if unwrap_poly_ys[0] < unwrap_poly_ys[-1]:
		unwrap_poly_ys.insert(0, 1);
		unwrap_poly_ys.append(imp.getHeight());
	else:
		unwrap_poly_ys.insert(0, imp.getHeight());
		unwrap_poly_ys.append(1);
	unwrap_axis = [(x, y) for (x, y) in zip(unwrap_poly_xs, unwrap_poly_ys)];
	
	unwrapped_projection_imp = do_unwrap(imp, unwrap_axis);
	return unwrapped_projection_imp, unwrap_axis;

def do_unwrap(imp, unwrap_axis, colorbar_width=200):
	"""given an unwrap axis extending from top to bottom of an image, generate an unwrapped image"""
	# extend to right hand corners...
	unwrap_poly_xs = [x for x, y in unwrap_axis];
	unwrap_poly_ys = [y for x, y in unwrap_axis];
	unwrap_poly_xs.append(imp.getWidth());
	unwrap_poly_xs.append(imp.getWidth());
	unwrap_poly_ys.append(unwrap_poly_ys[-1]);
	unwrap_poly_ys.append(unwrap_poly_ys[0]);
	unwrap_roi = PolygonRoi(unwrap_poly_xs, unwrap_poly_ys, Roi.POLYGON);

	input_title = imp.getTitle();
	unwrap_imp = IJ.createImage("unwrap", imp.getWidth(), imp.getHeight(), 2, imp.getBitDepth());
	left_crop = unwrap_roi.getBoundingRect().x;
	dummy = Duplicator().run(imp);
	dummy.setRoi(unwrap_roi);
	IJ.run(dummy, "Make Inverse", "");
	IJ.run(dummy, "Set...", "value=0");
	ip = dummy.getProcessor();
	unwrap_imp.setProcessor(ip);
	dummy.close();
	unwrap_imp.setZ(2);
	dummy = Duplicator().run(imp);
	dummy.setRoi(unwrap_roi);
	IJ.run(dummy, "Set...", "value=0");
	ip = dummy.getProcessor();
	unwrap_imp.setProcessor(ip);
	dummy.close();
	tile_imp = MontageMaker().makeMontage2(unwrap_imp, 2, 1, 1, 1, 2, 1, 0, False);
	# do thresholding/binarisation/profile plot to work out where to crop tight? otherwise...
	crop_roi = Roi(left_crop, 1, tile_imp.getWidth() - left_crop + colorbar_width, unwrap_imp.getHeight());
	unwrap_imp.close()
	tile_imp.setRoi(crop_roi);
	final_imp = tile_imp.crop();
	tile_imp.close()
	final_imp.setTitle(input_title + " twisted and unwrapped")
#	final_imp.show();
	return final_imp;

	
def do_angular_projection(imp, max_r_pix=60, min_r_pix=10, generate_roi_stack=False):
	"""perform ray-based projection of vessel wall, c.f. ICY TubeSkinner (Lancino 2018)"""
	split_chs = ChannelSplitter().split(imp);
	mch_imp = split_chs[0];
	IJ.setAutoThreshold(mch_imp, "IsoData dark stack");
	egfp_imp = split_chs[1];
	proj_imp = Duplicator().run(egfp_imp);
	cl_imp = split_chs[2];
	if generate_roi_stack:
		egfp_imp_disp = Duplicator().run(egfp_imp);
		roi_stack = IJ.createImage("rois", egfp_imp.getWidth(), egfp_imp.getHeight(), egfp_imp.getNSlices(), 16);

	centres = [];
	projected_im_pix = [];
	ring_rois = [];
	for zidx in range(cl_imp.getNSlices()):
		if ((zidx+1) % 100)==0:
			print("Progress = " + str(round(100*(float(zidx+1)/cl_imp.getNSlices()))));
		projected_im_row = [];
		proj_imp.setZ(zidx+1);
		mch_imp.setZ(zidx+1);
		bp = mch_imp.createThresholdMask();
		bp.dilate();
		bp.erode();
		bp.erode();
		bp.erode();
		mask_imp = ImagePlus("mask", bp);
		IJ.run(mask_imp, "Create Selection", "");
		roi = mask_imp.getRoi();
		proj_imp.setRoi(roi);
		IJ.run(proj_imp, "Set...", "value=0 slice");
		IJ.run(proj_imp, "Make Inverse", "");
		roi = proj_imp.getRoi();
		centre = (roi.getStatistics().xCentroid, roi.getStatistics().yCentroid);
		centres.append(centre);
		ring_roi_xs = [];
		ring_roi_ys = [];
		for theta in range(360):
			pt1 = (centre[0] + min_r_pix * math.cos(math.radians(theta)), 
					centre[1] + min_r_pix * math.sin(math.radians(theta)));
			pt2 = (centre[0] + max_r_pix * math.cos(math.radians(theta)), 
					centre[1] + max_r_pix * math.sin(math.radians(theta)));
			roi = Line(pt1[0], pt1[1], pt2[0], pt2[1]);
			proj_imp.setRoi(roi);
			profile = roi.getPixels();
			projected_im_row.append(max(profile));
			try:
				ring_roi_xs.append(roi.getContainedPoints()[profile.index(max(profile))].x);
			except IndexError:
				ring_roi_xs.append(pt2[0]);
			try:
				ring_roi_ys.append(roi.getContainedPoints()[profile.index(max(profile))].y);
			except IndexError:
				ring_roi_ys.append(pt2[1]);
			proj_imp.killRoi();
		ring_roi = PolygonRoi(ring_roi_xs, ring_roi_ys, Roi.FREELINE);
		ring_rois.append(ring_roi);
		if generate_roi_stack:
			roi_stack.setZ(zidx+1);
			roi_stack.setRoi(ring_roi);
			IJ.run(roi_stack, "Line to Area", "");
			IJ.run(roi_stack, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");
		#egfp_imp.setRoi(ring_roi);
		projected_im_pix.append(projected_im_row);
		
#	for ch in split_chs:
#		ch.close();

	out_imp = ImagePlus("projected", FloatProcessor([list(x) for x in zip(*projected_im_pix)]));

	if generate_roi_stack:
		roi_stack.show();
		egfp_imp_disp.show();
		# merge?
	else:
		roi_stack = None;
	return out_imp, roi_stack, ring_rois, centres

def calculate_mean_r(imp, ring_rois, centres):
	"""calculate the average distance of the cell surface from the vessel axis"""
	w = imp.getWidth();
	h = imp.getHeight();
	rp = FloatProcessor(w, h)
	rpix = rp.getPixels();
	for lidx, (roi, centre) in enumerate(zip(ring_rois, centres)):
		for thetaidx, (x, y) in enumerate(zip(roi.getPolygon().xpoints, roi.getPolygon().xpoints)):
			rpix[lidx * w + thetaidx] = math.sqrt((x - centre[0])**2 + (y - centre[1])**2);
	rimp = ImagePlus("Radii", rp);
	IJ.setAutoThreshold(rimp, "Intermodes light");
	bp = rimp.createThresholdMask();
	bp.dilate();
	bp.erode();
	mask_imp = ImagePlus("Mask", bp);
	IJ.run(mask_imp, "Create Selection", "");
	roi = mask_imp.getRoi();
	rimp.setRoi(roi);
	mean_r = rimp.getStatistics().mean;
	rimp.close();
	mask_imp.close();
	return mean_r;

def generate_r_image(imp, ring_rois, centres, unwrap_axis, threshold_val):
	"""for each point in the projection, calculate the distance to the vessel axis and present as an image"""
	fp = imp.getProcessor()
	fp.setThreshold(threshold_val, fp.maxValue(), FloatProcessor.NO_LUT_UPDATE);
	bp = fp.createMask();
	bp.dilate();
	bp.erode();

	mask_imp = ImagePlus("Mask", bp);
#	mask_imp.show();
#	WaitForUserDialog("pasue - generated mask").show();
	mask_imp = do_unwrap(mask_imp, unwrap_axis);
#	WaitForUserDialog("pasue - unwrapped").show();
	IJ.run(mask_imp, "Fill Holes", "");
#	WaitForUserDialog("pasue - filled holes").show();
	IJ.run(mask_imp, "Divide...", "value=255");
#	WaitForUserDialog("pasue - scaled to 0-1").show();
#	mask_imp.show();

	r_list = [];
	for lidx, (roi, centre) in enumerate(zip(ring_rois, centres)):
		r_sublist = [math.sqrt((x - centre[0])**2 + (y - centre[1])**2) for x, y in zip(roi.getPolygon().xpoints, roi.getPolygon().ypoints)];
		r_list.append(r_sublist);
		
	r_imp = ImagePlus("Radii", FloatProcessor([list(x) for x in zip(*r_list)]));
	r_imp = do_unwrap(r_imp, unwrap_axis);
	r_imp = ImageCalculator().run("Multiply create", r_imp, mask_imp);
	IJ.run(r_imp, "Cyan Hot", "");
	
	return r_imp, mask_imp;

def calculate_area_and_aspect_ratio(r_imp, mask_imp, raw_voxel_side):
	"""return the area and aspect ratio of the cell based on the projected image and convert to real-world units"""
	mask_imp.show();
	bp = mask_imp.getProcessor();
	bp.setThreshold(0.5, bp.maxValue(), FloatProcessor.NO_LUT_UPDATE);
	IJ.run(mask_imp, "Create Selection", "");
	roi = mask_imp.getRoi();
	r_imp.setRoi(roi);
	roi = r_imp.getRoi();
	raw_height = roi.getBounds().height;
	min_y = roi.getBounds().y;
	widths = [0] * raw_height;
	fp = r_imp.getProcessor();
	for point in roi:
		print(fp.getPixelValue(point.x, point.y));
		widths[point.y - min_y] = widths[point.y - min_y] + fp.getPixelValue(point.x, point.y) * raw_voxel_side * 2 * math.pi/360;
	aspect_ratio_circumferential_to_axial = (sum(widths)/len(widths)) / (raw_height * raw_voxel_side);
	area = sum(widths) * raw_voxel_side;
	
	return area, aspect_ratio_circumferential_to_axial;
	
	

def do_slicewise_unwrap(imp):
	"""at each position along the lumen axis, get a shell from around the lumen and take the maximum in the cell channel"""
	split_chs = ChannelSplitter().split(imp);
	mch_imp = split_chs[0];
	egfp_imp = split_chs[1];
	cl_imp = split_chs[2];

	mch_imp.show();
	egfp_imp.show();

	IJ.setAutoThreshold(mch_imp, "IsoData light stack");
	
	centres = [];
	for zidx in range(200, 202):
		print(zidx)
		mch_imp.setZ(zidx + 1);
		cl_imp.setZ(zidx + 1);
		ip = cl_imp.getProcessor();
		maxpt = MaximumFinder().getMaxima(ip, 10, True)
		centre = (maxpt.xpoints[0], maxpt.ypoints[0]);
		centres.append(centre);
		bp = mch_imp.createThresholdMask();
		bp.erode();
		bp.dilate();
		mask_imp = ImagePlus("mask", bp);
		IJ.run(mask_imp, "Create Selection", "");
		roi = mask_imp.getRoi();
		print([(x, y) for x, y in zip(roi.getInterpolatedPolygon().xpoints, roi.getPolygon().ypoints)]);

	
	mask_imp.show();
	cl_imp.close()

Prefs.blackBackground = True;

out_imp, _, ring_rois, centres = do_angular_projection(imp, generate_roi_stack=True);
IJ.setAutoThreshold(out_imp, "Intermodes dark");
threshold_val = out_imp.getProcessor().getMinThreshold();
unwrapped_projection_imp, unwrap_axis = twist_and_unwrap(out_imp);
r_imp, mask_imp = generate_r_image(out_imp, ring_rois, centres, unwrap_axis, threshold_val);

r_imp.show();
out_imp.close();
unwrapped_projection_imp.show();
area, aspect_ratio = calculate_area_and_aspect_ratio(r_imp, mask_imp, 0.108333);
print("A = " + str(area));
print("Aspect ratio = " + str(aspect_ratio));

