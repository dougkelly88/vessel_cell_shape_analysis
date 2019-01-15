# @ImagePlus imp
from ij import IJ
from ij.gui import PolygonRoi, Roi, WaitForUserDialog
from ij.plugin import Duplicator, MontageMaker

def twist_and_unwrap(imp):
	"""from the output of angular projection, define an axis along which to unzip the vessel and return this unzipped image"""
	IJ.setTool("freeline");
	WaitForUserDialog("Input required...", "Please draw the line down which unwrapping should occur...").show();
	roi = imp.getRoi();
	if roi is None:
		raise ValueError;
	unwrap_imp = IJ.createImage("unwrap", imp.getWidth(), imp.getHeight(), 2, 32);
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
	# extend to right hand corners...
	unwrap_poly_xs.append(imp.getWidth());
	unwrap_poly_xs.append(imp.getWidth());
	unwrap_poly_ys.append(unwrap_poly_ys[-1]);
	unwrap_poly_ys.append(unwrap_poly_ys[0]);
	new_roi = PolygonRoi(unwrap_poly_xs, unwrap_poly_ys, Roi.POLYGON);
	left_crop = new_roi.getBoundingRect().x;
	dummy = Duplicator().run(imp);
	dummy.setRoi(new_roi);
	IJ.run(dummy, "Make Inverse", "");
	IJ.run(dummy, "Set...", "value=0");
	ip = dummy.getProcessor();
	unwrap_imp.setProcessor(ip);
	dummy.close();
	unwrap_imp.setZ(2);
	dummy = Duplicator().run(imp);
	dummy.setRoi(new_roi);
	IJ.run(dummy, "Set...", "value=0");
	ip = dummy.getProcessor();
	unwrap_imp.setProcessor(ip);
	dummy.close();
	tile_imp = MontageMaker().makeMontage2(unwrap_imp, 2, 1, 1, 1, 2, 1, 0, False);
	crop_roi = Roi(left_crop, 1, tile_imp.getWidth() - left_crop, unwrap_imp.getHeight());
	unwrap_imp.close()
	tile_imp.setRoi(crop_roi);
	final_imp = tile_imp.crop();
	tile_imp.close()
	final_imp.setTitle("Twisted and unwrapped")
	final_imp.show();
	return final_imp, unwrap_axis;

_, unwrap_axis = twist_and_unwrap(imp);
print(unwrap_axis);