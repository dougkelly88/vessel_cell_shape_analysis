# @ImagePlus imp
import math, sys, os
from ij import ImageStack, ImagePlus, IJ
from ij.gui import WaitForUserDialog, Roi, PolygonRoi, Line
from ij.plugin import ChannelSplitter, Straightener, Duplicator
from ij.plugin.filter import MaximumFinder
from ij.process import FloatProcessor

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "Unwrapping and shape analysis";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

from UpdateRoiImageListener import UpdateRoiImageListener

# angle method
max_r_pix = 60;
min_r_pix = 10;

split_chs = ChannelSplitter().split(imp);
mch_imp = split_chs[0];
egfp_imp = split_chs[1];
egfp_imp_disp = Duplicator().run(egfp_imp);
cl_imp = split_chs[2];
#egfp_imp.show();

centres = [];
projected_im_pix = [];
ring_rois = [];
unzip_axis = [];
roi_stack = IJ.createImage("rois", egfp_imp.getWidth(), egfp_imp.getHeight(), egfp_imp.getNSlices(), 16);
#roi_stack.show();

for zidx in range(cl_imp.getNSlices()):
#for zidx in range(200,350):
#for zidx in range(200,202):
	if ((zidx+1) % 100)==0:
		print("Z = " + str(zidx+1));
	projected_im_row = [];
	egfp_imp.setZ(zidx+1);
	cl_imp.setZ(zidx+1);
	ip = cl_imp.getProcessor();
	out = MaximumFinder().getMaxima(ip, 10, True)
	centres.append((out.xpoints[0], out.ypoints[0]));
	centre = (out.xpoints[0], out.ypoints[0]);
	ring_roi_xs = [];
	ring_roi_ys = [];
	for theta in range(359):
#	for theta_idx in range(6):
#		theta = theta_idx * 60 - 90;
		pt1 = (centre[0] + min_r_pix * math.cos(math.radians(theta)), 
				centre[1] + min_r_pix * math.sin(math.radians(theta)));
		pt2 = (centre[0] + max_r_pix * math.cos(math.radians(theta)), 
				centre[1] + max_r_pix * math.sin(math.radians(theta)));
		roi = Line(pt1[0], pt1[1], pt2[0], pt2[1]);
		egfp_imp.setRoi(roi);
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
		#print("Max val = " + str(max(profile)));
		egfp_imp.killRoi();
	ring_roi = PolygonRoi(ring_roi_xs, ring_roi_ys, Roi.FREELINE);
	ring_rois.append(ring_roi);
	roi_stack.setZ(zidx+1);
	roi_stack.setRoi(ring_roi);
	IJ.run(roi_stack, "Line to Area", "");
	IJ.run(roi_stack, "Set...", "value=" + str(roi_stack.getProcessor().maxValue()) + " slice");
	min_idx = projected_im_row.index(min(projected_im_row));
#	print("min idx = " + str(min_idx));
	unzip_axis.append(min_idx);
	egfp_imp.setRoi(ring_roi);
	projected_im_pix.append(projected_im_row);
#	WaitForUserDialog("pause").show();
		
#print(centres);
for ch in split_chs:
	ch.close();

w = len(projected_im_row);
h = len(projected_im_pix);
ip = FloatProcessor(w, h);
pix = ip.getPixels();
for x in range(w):
	for y in range(h):
		pix[y * w + x] = projected_im_pix[y][x];

out_imp = ImagePlus("projected", ip);
out_imp.show();

egfp_imp_disp.show();
roi_stack.show();
#egfp_imp_disp.addImageListener(UpdateRoiImageListener(ring_rois));

unwrapped_imp = IJ.createImage("twisted unwrap", 2*out_imp.getWidth(), out_imp.getHeight(), 1, 32);
unwrapped_imp.show();
IJ.setTool("freeline");
WaitForUserDialog("input unwrap axis").show();
unwrap_poly = out_imp.getRoi().getPolygon();
unwrap_poly_xs

#print("unzip axis = " + str(unzip_axis));
#right_side_xs = unzip_axis;
#right_side_xs.append(out_imp.getWidth());
#right_side_xs.append(out_imp.getWidth());
#right_side_ys = [y+1 for y in range(out_imp.getHeight())];
#right_side_ys.append(out_imp.getHeight());
#right_side_ys.append(1);
#right_side_roi = PolygonRoi(right_side_xs, right_side_ys, Roi.POLYGON);
#out_imp.setRoi(right_side_roi);


