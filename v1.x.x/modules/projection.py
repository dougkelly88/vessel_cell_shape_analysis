import math, sys, os
from ij import ImageStack, ImagePlus, IJ, Prefs
from ij.gui import PolygonRoi, Roi, WaitForUserDialog, Line, ProfilePlot
from ij.plugin import Duplicator, MontageMaker, ChannelSplitter, ImageCalculator
from ij.process import FloatProcessor, AutoThresholder

def make_tiled_imp(imp):
	"""generate a ImageProcessor that is the input tiled 3 time horizontally"""
	ip = imp.getProcessor();
	stack = ImageStack(imp.getWidth(), imp.getHeight());
	for idx in range(3):
		stack.addSlice(ip);
	temp_imp = ImagePlus("temp", stack);
	tile_imp = MontageMaker().makeMontage2(temp_imp, 3, 1, 1, 1, 3, 1, 0, False);
	temp_imp.close();
	return tile_imp;

def twist_and_unwrap(imp):
	"""from the output of angular projection, define an axis along which to unzip the vessel and return this unzipped image"""
	tile_imp = make_tiled_imp(imp);
	tile_imp.show();
	IJ.setTool("freeline");
	WaitForUserDialog("Input required...", "Please draw the line down which unwrapping should occur...").show();
	roi = tile_imp.getRoi();
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
	
	unwrapped_projection_imp = do_unwrap(tile_imp, unwrap_axis, imp_title=imp.getTitle());
	return unwrapped_projection_imp, unwrap_axis;

def do_unwrap(tile_imp, unwrap_axis, imp_title=None):
	"""given an unwrap axis extending from top to bottom of an image, generate an unwrapped image"""
	ip = tile_imp.getProcessor();
	ip.setValue(0);
	unwrap_poly_xs = [x for x, y in unwrap_axis];
	unwrap_poly_ys = [y for x, y in unwrap_axis];
	if sum(unwrap_poly_xs)/len(unwrap_poly_xs) > 1.5 * 360:
		rhs_poly_xs = unwrap_poly_xs;
	else:
		rhs_poly_xs = [x + 360 for x in unwrap_poly_xs];
	for _ in range(2):
		rhs_poly_xs.append(0);
	unwrap_poly_ys.append(unwrap_poly_ys[-1]);
	unwrap_poly_ys.append(unwrap_poly_ys[0]);
	crop_roi = PolygonRoi(rhs_poly_xs, unwrap_poly_ys, Roi.POLYGON);
	ip.fillOutside(crop_roi);
	ip.setRoi(crop_roi);
	ip = ip.crop();
	ip.setValue(0);
	lhs_poly_xs = [x - 360 for x in rhs_poly_xs[:-2]];
	for _ in range(2):
		lhs_poly_xs.append(ip.getWidth());
	crop_roi = PolygonRoi(lhs_poly_xs, unwrap_poly_ys, Roi.POLYGON);
	ip.fillOutside(crop_roi);
	ip.setRoi(crop_roi);
	ip = ip.crop();
	tile_imp.setProcessor(ip);
	tile_imp.updateAndRepaintWindow();
	if imp_title is not None:
		tile_imp.setTitle("{}, twisted and unwrapped".format(imp_title));
	else:
		tile_imp.setTitle("twisted and unwrapped")
	return tile_imp;
	
def do_angular_projection(imp, max_r_pix=60, min_r_pix=10, generate_roi_stack=True):
	"""perform ray-based projection of vessel wall, c.f. ICY TubeSkinner (Lancino 2018)"""
	Prefs.blackBackground = True;
	print("do angular projection input imp = " + str(imp));
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
	tile_mask = make_tiled_imp(mask_imp);
	#tile_mask.show();
	#WaitForUserDialog("pasue - generated mask").show();
	mask_imp = do_unwrap(tile_mask, unwrap_axis, imp_title=mask_imp.getTitle());
	#mask_imp.show();
	roi = PolygonRoi([x for (x, y) in unwrap_axis], [y for (x, y) in unwrap_axis], PolygonRoi.POLYLINE);
	mask_imp.setRoi(roi);
	#WaitForUserDialog("pasue - unwrapped").show();
	IJ.run(mask_imp, "Fill Holes", "");
	#WaitForUserDialog("pasue - filled holes").show();
	IJ.run(mask_imp, "Divide...", "value=255");
	#WaitForUserDialog("pasue - scaled to 0-1").show();
	#mask_imp.show();

	r_list = [];
	for lidx, (roi, centre) in enumerate(zip(ring_rois, centres)):
		r_sublist = [math.sqrt((x - centre[0])**2 + (y - centre[1])**2) for x, y in zip(roi.getPolygon().xpoints, roi.getPolygon().ypoints)];
		r_list.append(r_sublist);
		
	r_imp = ImagePlus("Radii", FloatProcessor([list(x) for x in zip(*r_list)]));
	tile_r_imp = make_tiled_imp(r_imp);
	r_imp = do_unwrap(tile_r_imp, unwrap_axis, imp_title=r_imp.getTitle());
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

#Prefs.blackBackground = True;

#out_imp, _, ring_rois, centres = do_angular_projection(imp, generate_roi_stack=True);
#IJ.setAutoThreshold(out_imp, "Intermodes dark");
#threshold_val = out_imp.getProcessor().getMinThreshold();
#unwrapped_projection_imp, unwrap_axis = twist_and_unwrap(out_imp);
#r_imp, mask_imp = generate_r_image(out_imp, ring_rois, centres, unwrap_axis, threshold_val);

#r_imp.show();
#out_imp.close();
#unwrapped_projection_imp.show();
#area, aspect_ratio = calculate_area_and_aspect_ratio(r_imp, mask_imp, 0.108333);
#print("A = " + str(area));
#print("Aspect ratio = " + str(aspect_ratio));
#
#test_path = "C:\\Users\\dougk\\Desktop\\tile_imp_test.tif";
#tile_imp = IJ.openImage(test_path);
#tile_imp.show();
#IJ.setTool("freeline");
#WaitForUserDialog("draw line").show();
#roi  = tile_imp.getRoi();
#unwrap_poly = roi.getPolygon();
#unwrap_poly_xs = [x for x in unwrap_poly.xpoints];
#unwrap_poly_ys = [y for y in unwrap_poly.ypoints];
## extend to the top and bottom of the image:
#unwrap_poly_xs.insert(0, unwrap_poly_xs[0]);
#unwrap_poly_xs.append(unwrap_poly_xs[-1])
#if unwrap_poly_ys[0] < unwrap_poly_ys[-1]:
#	unwrap_poly_ys.insert(0, 1);
#	unwrap_poly_ys.append(tile_imp.getHeight());
#else:
#	unwrap_poly_ys.insert(0, tile_imp.getHeight());
#	unwrap_poly_ys.append(1);
#unwrap_axis = [(x, y) for (x, y) in zip(unwrap_poly_xs, unwrap_poly_ys)];
#out_imp = do_unwrap(tile_imp, unwrap_axis)