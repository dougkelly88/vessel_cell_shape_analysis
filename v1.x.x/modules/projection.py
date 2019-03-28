import math, sys, os
from ij import ImageStack, ImagePlus, IJ, Prefs
from ij.gui import PolygonRoi, Roi, WaitForUserDialog, Line, ProfilePlot, NonBlockingGenericDialog
from ij.plugin import Duplicator, MontageMaker, ChannelSplitter, ImageCalculator, RGBStackMerge
from ij.plugin.filter import ParticleAnalyzer
from ij.process import FloatProcessor, AutoThresholder
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
#from ij.io import FileSaver

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
	MyWaitForUser("Input required...", "Please draw the line down which unwrapping should occur...");
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
	
def do_angular_projection(imp, output_path, vessel_axis, max_r_pix=60, min_r_pix=10, generate_roi_stack=False, smooth_radius_pix=1):
	"""perform ray-based projection of vessel wall, c.f. ICY TubeSkinner (Lancino 2018)"""
	input_n_ch = imp.getNChannels();
	Prefs.blackBackground = True;
	print("do angular projection input imp = " + str(imp));
	split_chs = ChannelSplitter().split(imp);
	mch_imp = split_chs[0];
	IJ.setAutoThreshold(mch_imp, "IsoData dark stack");
	egfp_imp = split_chs[1];
	proj_imp = Duplicator().run(egfp_imp);
	#if input_n_ch > 2:
	#	axis_imp = split_chs[2];
	if generate_roi_stack:
		egfp_imp_disp = Duplicator().run(egfp_imp);
		roi_stack = IJ.createImage("rois", egfp_imp.getWidth(), egfp_imp.getHeight(), egfp_imp.getNSlices(), 16);

	centres = [];
	projected_im_pix = [];
	ring_rois = [];
	centre = (proj_imp.getWidth()/2, proj_imp.getHeight()/2)
	for z in range(proj_imp.getNSlices()):
		z = z+1;
		print("z = {}".format(z));
		if ((z) % 100)==0:
			print("Progress = {}%".format(round(100*(float(z)/proj_imp.getNSlices()))));
		projected_im_row = [];
		proj_imp.setZ(z);
		centres.append(centre);
		#if input_n_ch < 3:
		#	mch_imp.setZ(z);
		#	bp = mch_imp.createThresholdMask();
		#	bp.dilate();
		#	bp.erode();
		#	bp.erode();
		#	bp.erode();
		#	mask_imp = ImagePlus("mask", bp);
		#	mask_imp.show();
		#	#MyWaitForUser("pause", "pause after mask");
		#	IJ.run(mask_imp, "Create Selection", "");
		#	roi = mask_imp.getRoi();
		#	proj_imp.setRoi(roi);
		#	IJ.run(proj_imp, "Set...", "value=0 slice");
		#	proj_imp.setRoi(roi);
		#	IJ.run(proj_imp, "Make Inverse", "");
		#	roi = proj_imp.getRoi();
		#	centre = (roi.getStatistics().xCentroid, roi.getStatistics().yCentroid);
		#	#print("centre = {}".format(centre));
		#	centres.append(centre);
		#else:
		#	axis_imp.setZ(z);
		#	ip = axis_imp.getProcessor();
		#	ip.setThreshold(255, ip.maxValue(), FloatProcessor.NO_LUT_UPDATE);
		#	bp = axis_imp.createThresholdMask();
		#	bp.dilate();
		#	mask_imp = ImagePlus("mask", bp);
		#	#mask_imp.show();
		#	IJ.run(mask_imp, "Create Selection", "");
		#	roi = mask_imp.getRoi();
		#	centre = (roi.getStatistics().xCentroid, roi.getStatistics().yCentroid);
		#	#print("centre = {}".format(centre));
		#	centres.append(centre);
		#	mask_imp.close();
		##proj_imp.show();
		#if roi is not None:
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
		#else:
		#	ring_roi = None;
		#	ring_rois.append(None);
		#	projected_im_row = [0 for _ in range(360)];
		if generate_roi_stack:
			roi_stack.setZ(z);
			#merged_radius_imp.setZ(zidx+1);
			if ring_roi is not None:
				roi_stack.setRoi(ring_roi);
				#merged_radius_imp.setRoi()
				IJ.run(roi_stack, "Line to Area", "");
				IJ.run(roi_stack, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");
				#IJ.run(merged_radius_imp, "Line to Area", "");
				#IJ.run(merged_radius_imp, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");
		#egfp_imp.setRoi(ring_roi);
		projected_im_pix.append(projected_im_row);
		
#	for ch in split_chs:
#		ch.close();

	out_imp = ImagePlus("projected", FloatProcessor([list(x) for x in zip(*projected_im_pix)]));
	IJ.run(out_imp, "Median...", "radius={}".format(smooth_radius_pix));

	if generate_roi_stack:
		roi_stack.show();
		#egfp_imp_disp.show();
		#mch_imp.show();
		merged_radius_imp = RGBStackMerge().mergeChannels([mch_imp, egfp_imp_disp, None, None, None, None, roi_stack], True)
		IJ.run(merged_radius_imp, "16-bit", "");
		merged_radius_imp.show();
		#IJ.saveAsTiff(merged_radius_imp, os.path.join(output_path, "radius ROIs.tif"))
		#FileSaver(merged_radius_imp).saveAsTiffStack(os.path.join(output_path, "identified vessel radius along vessel axis.tif"));
	else:
		roi_stack = None;
	roi_stack.changes = False;
	roi_stack.close();
	return out_imp, merged_radius_imp, ring_rois, centres

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

def generate_r_image(imp, ring_rois, centres, unwrap_axis, threshold_val, smooth_radius_pix=1):
	"""for each point in the projection, calculate the distance to the vessel axis and present as an image"""
	imp.show();
	fp = imp.getProcessor();
	fp.setThreshold(threshold_val, fp.maxValue(), FloatProcessor.RED_LUT);
	print("thresh val = {}".format(threshold_val));
	print("type fp = {}".format(type(fp)));
	bp = fp.createMask();
	bp.dilate();
	bp.erode();
	mask_imp = ImagePlus("Mask", bp);
	tile_mask = make_tiled_imp(mask_imp);
	IJ.run(tile_mask, "Fill Holes", "");
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
	bp = mask_imp.getProcessor();
	bp.setThreshold(0.5, bp.maxValue(), FloatProcessor.NO_LUT_UPDATE);
	keep_largest_blob(mask_imp);

	r_list = [];
	for lidx, (roi, centre) in enumerate(zip(ring_rois, centres)):
		if roi is None:
			r_sublist = [0 for _ in range(360)];
		else:
			r_sublist = [math.sqrt((x - centre[0])**2 + (y - centre[1])**2) for x, y in zip(roi.getPolygon().xpoints, roi.getPolygon().ypoints)];
		r_list.append(r_sublist);
		
	r_imp = ImagePlus("Radii", FloatProcessor([list(x) for x in zip(*r_list)]));
	tile_r_imp = make_tiled_imp(r_imp);
	r_imp = do_unwrap(tile_r_imp, unwrap_axis, imp_title=r_imp.getTitle());
	IJ.run(r_imp, "Median...", "radius={}".format(smooth_radius_pix));
	r_imp_masked = ImageCalculator().run("Multiply create", r_imp, mask_imp);
	IJ.run(r_imp_masked, "Cyan Hot", "");
	IJ.run(mask_imp, "Green", "");
	return r_imp_masked, mask_imp, r_imp;

def keep_largest_blob(imp):
	"""remove all blobs other than the largest by area"""
	rt = ResultsTable();
	mxsz = imp.width * imp.height;
	roim = RoiManager(False);
	pa = ParticleAnalyzer(ParticleAnalyzer.ADD_TO_MANAGER, ParticleAnalyzer.AREA | ParticleAnalyzer.SLICE, rt, 0, mxsz);
	pa.setRoiManager(roim);
	
	for idx in range(1, imp.getImageStackSize()+1):
		roim.reset();
		rt.reset();
		imp.setPosition(idx);
		pa.analyze(imp);
		rt_areas = rt.getColumn(rt.getColumnIndex("Area")).tolist();
		mx_ind = rt_areas.index(max(rt_areas))
		indices_to_remove = [a for a in range(0,len(rt_areas)) if a != mx_ind];
		indices_to_remove.reverse();
		for rem_idx in indices_to_remove:
			roim.select(imp, rem_idx);
			IJ.run(imp, "Set...", "value=0 slice");
	imp.killRoi();
	roim.reset();
	roim.close();

def calculate_area_and_aspect_ratio(r_imp, mask_imp, raw_voxel_side):
	"""return the area and aspect ratio of the cell based on the projected image and convert to real-world units"""
	mask_imp.show();
	MyWaitForUser("pause", "pause at mask imp before calculation");
	bp = mask_imp.getProcessor();
	bp.setThreshold(0.5, bp.maxValue(), FloatProcessor.NO_LUT_UPDATE);
	keep_largest_blob(mask_imp);
	MyWaitForUser("pause", "pause after keep largest blob");
	IJ.run(mask_imp, "Create Selection", "");
	roi = mask_imp.getRoi();
	r_imp.setRoi(roi);
	roi = r_imp.getRoi();
	MyWaitForUser("pause", "pause after getting roi from r_imp");
	# MyWaitForUser("pause", "pause to check mask has been filtered by size");
	stats = roi.getStatistics();
	print("stats = {}".format(stats));
	raw_height = roi.getBounds().height;
	min_y = roi.getBounds().y;
	widths = [0] * raw_height;
	fp = r_imp.getProcessor();
	for point in roi:
		widths[point.y - min_y] = widths[point.y - min_y] + fp.getPixelValue(point.x, point.y) * raw_voxel_side * 2 * math.pi/360;
		#widths[point.y - min_y] = widths[point.y - min_y] + 50 * raw_voxel_side * 2 * math.pi/360; # debug
	aspect_ratio_circumferential_to_axial = (sum(widths)/len(widths)) / (raw_height * raw_voxel_side);
	aspect_ratio_fitted_ellipse_IJ = stats.major/stats.minor;
	area = sum(widths) * raw_voxel_side;
	
	return area, aspect_ratio_fitted_ellipse_IJ, stats.min * raw_voxel_side, stats.max * raw_voxel_side;
	
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