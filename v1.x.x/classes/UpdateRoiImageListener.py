# class to support updating ROI from list upon change of slice and modification of saved ROIs
#
# D. J. Kelly, 2018-12-05, douglas.kelly@riken.jp

from ij import ImageListener

class UpdateRoiImageListener(ImageListener):
	"""class to support updating ROI from list upon change of frame"""
	def __init__(self, roi_list):
		self.last_slice = 1;
		self.roi_list = roi_list;
		print("UpdateRoiImageListener started");

	def imageUpdated(self, imp):
		roi = imp.getRoi();
		if roi is not None and not roi.isArea():
			self.roi_list[self.last_slice - 1] = roi
		if imp.getNSlices() > imp.getNFrames():
			slc = imp.getZ();
		else:
			slc = imp.getT();
		print(slc);
		self.last_slice = slc;
		if slc < len(self.roi_list):
			imp.setRoi(self.roi_list[slc - 1]);
		
	def imageOpened(self, imp):
		print("UpdateRoiImageListener: image opened");
			
	def imageClosed(self, imp):
		print("UpdateRoiImageListener: image closed");
		imp.removeImageListener(self);

	def getRoiList(self):
		return self.roi_list;
