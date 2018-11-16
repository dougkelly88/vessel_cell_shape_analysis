import os
from ij import IJ
from ij import WindowManager as WM
from ij.gui import WaitForUserDialog
from ij.io import OpenDialog
from ij.plugin import HyperStackConverter, Duplicator, ZProjector
from ij.process import AutoThresholder, StackStatistics
from loci.plugins import BF as bf

def file_location_chooser(default_directory):
	"""choose folder locations and prepare output folder"""
	# input
	od = OpenDialog('Choose original file...', 
					default_directory, 
					'*.tif');
	file_path = od.getPath();
	if file_path is None:
		raise IOError('no input file chosen');
	return file_path;

def z_crop(imp):
	"""trim a z stack based on interactively-defined start and end points"""
	WaitForUserDialog("Choose first z plane and click OK...").show();
	start_z = imp.getZ();
	WaitForUserDialog("Now choose last z plane and click OK...").show();
	end_z = imp.getZ();
	frames = imp.getNFrames();
	channels = imp.getNChannels();
	#current_channel = imp.getC();
	dupimp = Duplicator().run(imp, 1, channels, start_z, end_z, 1, frames);
	imp.changes = False;
	imp.close();
	dupimp.show()
	#autoset_zoom(dupimp);
	return dupimp, (start_z, end_z);

def import_iq3_metadata(metadata_path):
	"""import basic image metadata based on the metadata saved by iQ3 software at acquisition time"""
	import re
	x_fmt_str = 'x \: (?P<x_pixels>\d+\.?\d*) \* (?P<x_physical_size>\d+\.?\d*) \: (?P<x_unit>\w+)';
	y_fmt_str = 'y \: (?P<y_pixels>\d+\.?\d*) \* (?P<y_physical_size>\d+\.?\d*) \: (?P<y_unit>\w+)';
	z_fmt_str = '\s*Repeat Z \- (?P<z_extent>[+-]?\d+\.?\d*) (?P<z_unit>\w+) in (?P<z_pixels>\d+\.?\d*) planes \(centre\)';
	t_fmt_str = '\s*Repeat T \- (?P<n_frames>\d+\.?\d*) times \((?P<frame_interval>\d+\.?\d*) sec\)'; # check if time is always in seconds (sec)?
	c_fmt_str = r"\s*Repeat \- Channel \((?P<raw_channels_str>\w.*)\)";
	format_strings = [x_fmt_str, y_fmt_str, z_fmt_str, t_fmt_str, c_fmt_str];
	
	meta_dict = {}
	metadata_file = open(metadata_path, 'r')
	try:
		for line in metadata_file.readlines():
			for fmt_str in format_strings:
				m = re.match(fmt_str, line)
				if (bool(m)):
					meta_dict.update(m.groupdict())
		p_num = re.compile('[+-]?\d+\.?\d*')
		for key, value in meta_dict.iteritems():
			if p_num.match(value):
				try:
					meta_dict[key] = float(value)
				except:
					#print("conversion error for key " + key);
					continue;
	finally:
		metadata_file.close();
		if 'raw_channels_str' in meta_dict:
			ch_list = meta_dict['raw_channels_str'].split(",")
			meta_dict['n_channels'] = len(ch_list);
			meta_dict['channel_list'] = ch_list;
		return meta_dict


default_path = "D:\\data\\Marcksl1 cell shape analysis\\zstacks\\2018-05-15 UAS-marcksl1b-delED-in-myrEGFP"
image_file_path = file_location_chooser(default_path)
imps = bf.openImagePlus(image_file_path);
imp = imps[0];
metadata_file_path = os.path.splitext(image_file_path)[0] + ".txt";
metadata = import_iq3_metadata(metadata_file_path);
#imp2 = HyperStackConverter.toHyperStack(imp, int(metadata['n_channels']), int(metadata['z_pixels']), 1, "Composite");
print(metadata);
IJ.run(imp, "Properties...", "channels=" + str(int(metadata['n_channels'])) + 
								" slices=" + str(int(metadata['z_pixels'])) + 
								" frames=1 unit=" + str(metadata['x_unit']) + 
								" pixel_width=" + str(metadata['x_physical_size']) + 
								" pixel_height=" + str(metadata['y_physical_size']) + 
								" voxel_depth=" + str(metadata['z_extent']/metadata['z_pixels']));
imp.show();
imp, z_lims = z_crop(imp);
use_ch = imp.getC();
print(use_ch);
IJ.setTool("rect");
imp.hide();
imp2 = ZProjector.run(imp, "max");
imp2.show();
imp2.setC(use_ch+1);
WaitForUserDialog("Select XY region to crop...").show();
crop_roi = imp2.getRoi();
imp2.close();
imp.show();
imp.setRoi(crop_roi);
IJ.run("Crop", "");
frames = imp.getNFrames();
slices = imp.getNSlices();
imp = Duplicator().run(imp, use_ch+1, use_ch+1, 1, slices, 1, frames);
imp.show();

# now do 3d segmentation....
#histo =  StackStatistics(imp).histogram;
#global_thr = AutoThresholder().getThreshold(AutoThresholder.Method.IJ_IsoData, histo) + 1;
#print(global_thr);
#IJ.run(imp, "3D Simple Segmentation", "low_threshold=" + str(global_thr) + " min_size=200 max_size=-1");
#segimp = WM.getImage("Seg");

