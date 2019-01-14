 # @ImagePlus imp
from ij import ImageStack, ImagePlus, IJ
from ij.gui import WaitForUserDialog, Roi, PolygonRoi
from ij.plugin import ChannelSplitter, Straightener, ZProjector

# run on imp containing result of tubefitter + IJ.run(imp, "Reslice [/]...", "output=1.000 start=Top avoid") + spline-fitted 
# manual roi copied from z-projected version of this resliced image (add roi from z-proj to ROI manager, then add 
# to resliced image)

imp.killRoi();
IJ.run(imp, "Reslice [/]...", "output=1.000 start=Top avoid");
rot_imp = IJ.getImage();
crop_roi = Roi(10, 0, rot_imp.getWidth()-20, rot_imp.getHeight());
rot_imp.setRoi(crop_roi);
rot_imp.crop();
rot_imp.show();
WaitForUserDialog("pause").show();
split_ch = ChannelSplitter().split(rot_imp);
mch_imp = split_ch[0];
egfp_imp = split_ch[1];
roi_imp = split_ch[2];
zproj_imp = ZProjector.run(roi_imp, "max");
IJ.setRawThreshold(zproj_imp, 33153, 65535, None);
IJ.run(zproj_imp, "Make Binary", "")
zproj_imp.show();
raise error;
IJ.run(zproj_imp, "Skeletonize", "");
IJ.run(zproj_imp, "Create Selection", "");
roi = zproj_imp.getRoi();
floatpoly = roi.getContainedFloatPoints();
roi = PolygonRoi(floatpoly, Roi.FREELINE);
zproj_imp.setRoi(roi);
WaitForUserDialog("pause").show();

#IJ.setTool("freeline");
#WaitForUserDialog("Draw a line").show();
#IJ.run(zproj_imp, "Fit Spline", "");
#roi = zproj_imp.getRoi();
zproj_imp.changes = False;
zproj_imp.close();



z_planes = egfp_imp.getNSlices();
#straight_stack = Straightener().straightenStack(egfp_imp, roi, 100);
#out_imp = ImagePlus("Straightened", straight_stack);

for zidx in range(z_planes):
	for chidx in range(3):
		#print("Working on " + str(["egfp", "mCh", "roi"][chidx]));
		split_ch[chidx].setZ(zidx+1);
		split_ch[chidx].setRoi(roi);
		ip = Straightener().straightenLine(split_ch[chidx], 150);
		if chidx==1:
			if zidx==0:
				egfp_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
			egfp_straight_stack.addSlice(ip);
		elif chidx==0:
			if zidx==0:
				mch_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
			mch_straight_stack.addSlice(ip);
		else:
			if zidx==0:
				roi_straight_stack = ImageStack(ip.getWidth(), ip.getHeight());
			roi_straight_stack.addSlice(ip);

egfp_out_imp = ImagePlus("Straightened EGFP", egfp_straight_stack);
mch_out_imp = ImagePlus("Straightened mCh", mch_straight_stack);
roi_out_imp = ImagePlus("Straightened ROI", roi_straight_stack);
egfp_out_imp.show();
mch_out_imp.show();
roi_out_imp.show();
IJ.run("Merge Channels...", "c1=[" + mch_out_imp.getTitle() + 
									"] c2=[" + egfp_out_imp.getTitle() + 
									"] c7=[" + roi_out_imp.getTitle() + "] create keep");
