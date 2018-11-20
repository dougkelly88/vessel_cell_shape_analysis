import os, sys, re
from datetime import datetime
from ij import IJ
from ij import WindowManager as WM
from ij.gui import WaitForUserDialog, GenericDialog, YesNoCancelDialog
from ij.io import OpenDialog, DirectoryChooser, FileSaver
from ij.plugin import HyperStackConverter, Duplicator, ZProjector
from ij.process import AutoThresholder, StackStatistics
import ij.WindowManager as WM
from loci.plugins import BF as bf

release = False;

if not release:
	script_path = os.path.dirname(os.path.realpath(__file__));
else: 
	script_path = os.getcwd();
if "Fiji.app" in script_path:
	ss = script_path.split("Fiji.app");
	final_folder = "marcksl1 shape prescreener";
	script_path = os.path.join(ss[0], "Fiji.app", "plugins", "Scripts", "Plugins", final_folder);
sys.path.insert(0, os.path.join(script_path, 'modules'));
sys.path.insert(0, os.path.join(script_path, 'classes'));

from PrescreenInfo import PrescreenInfo

def file_location_chooser(default_directory):
	"""choose input data location"""
	od = OpenDialog('Choose original file...', 
					default_directory, 
					'*.tif');
	file_path = od.getPath();
	if file_path is None:
		raise IOError('no input file chosen');
	return file_path;

def output_folder_chooser(default_directory):
	"""choose where output data should be saved"""
	DirectoryChooser.setDefaultDirectory(default_directory);
	dc = DirectoryChooser('Select the folder for saving cropped images...');
	output_root = dc.getDirectory();
	if output_root is None:
		raise IOError('no output path chosen');
	timestamp = datetime.strftime(datetime.now(), '%Y-%m-%d %H-%M-%S');
	output_folder = os.path.join(output_root, (timestamp + ' output'));
	os.mkdir(output_folder);
	return output_folder;
	

def z_crop(imp):
	"""trim a z stack based on interactively-defined start and end points"""
	imp.setZ(1);
	WaitForUserDialog("Choose first z plane and click OK...").show();
	start_z = imp.getZ();
	WaitForUserDialog("Now choose last z plane and click OK...").show();
	end_z = imp.getZ();
	frames = imp.getNFrames();
	channels = imp.getNChannels();
	imp.killRoi();
	dupimp = Duplicator().run(imp, 1, channels, start_z, end_z, 1, frames);
	imp.hide();
	dupimp.show()
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

def parse_info_from_filename(info):
	"""attempt to populate channel label and embryo ID fields from filename"""
	fmt_str = r'(?P<exp_id>\w.*) e(?P<emb_id>\w{1,3})20(?P<date>\d+)_(?P<time>\d+) (?P<ampm>[AP])M.tif';
	m = re.match(fmt_str, os.path.basename(info.get_input_file_path()));
	if bool(m):
		d = m.groupdict();
		info.set_experiment_id(d['exp_id']);
		info.set_embryo_id(d['emb_id']);
	return info;

def save_etc(imp, info, output_folder):
	"""handle saving output, final labelling and UI of what to do next"""
	dialog = GenericDialog("Marksl1 cell shape prescreen");
	dialog.addChoice("Vessel type: ", 
						info.get_vessel_types(), 
						info.get_vessel_type());
	dialog.addStringField("Channel 1 label: ", info.get_ch1_label());
	dialog.addStringField("Channel 2 label: ", info.get_ch2_label());
	dialog.addStringField("Experiment identifier: ", info.get_experiment_id());
	dialog.addStringField("Embryo identifier: ", info.get_embryo_id());
	dialog.setOKLabel("Save preprocessing");
	dialog.showDialog();
	info.set_vessel_type(dialog.getNextChoice());
	info.set_ch1_label(dialog.getNextString());
	info.set_ch2_label(dialog.getNextString());
	info.set_experiment_id(dialog.getNextString());
	info.set_embryo_id(dialog.getNextString());
	if dialog.wasCanceled():
		return 0;
	if dialog.wasOKed():
		exsting_files = os.listdir(output_folder);
		r = re.compile(".*" + info.get_embryo_id() + " " + info.get_vessel_type() + ".*")
		fns = filter(r.match, exsting_files)
		numbers = list((int(s) for fn in fns for s in re.findall(r'\b\d+$', os.path.splitext(fn)[0])));
		append_digit = (max(numbers) + 1) if numbers else 1;
		ch1str = (info.get_ch1_label() + " ") if info.get_ch1_label() else "";
		ch2str = (info.get_ch2_label() + " ") if info.get_ch2_label() else "";
		expstr = (info.get_experiment_id() + " ") if info.get_experiment_id() else "";
		file_name = ("Cropped " + ch1str + ch2str + expstr + 
					"e" + str(info.get_embryo_id()) + " " + 
					info.get_vessel_type() + " " + 
					str(append_digit));
		FileSaver(imp).saveAsTiff(os.path.join(output_folder, (file_name + ".tif")));
		info.save_info_to_json(os.path.join(output_folder, (file_name + ".json")));
		continueDialog = YesNoCancelDialog(WM.getCurrentWindow(), 
											"Continue?", 
											"Continue with same input image or a new input image?", 
											"New image...", 
											"Same image");
		if continueDialog.cancelPressed():
			return 0;
		if continueDialog.yesPressed():
			return 2;
		return 1;
	return 0;

def perform_cropping(imp, info, output_folder, default_path):
	if imp is None:
		info.set_input_file_path(file_location_chooser(default_path));
		imps = bf.openImagePlus(info.get_input_file_path());
		imp = imps[0];
	
		info.set_metadata_file_path(os.path.splitext(info.get_input_file_path())[0] + ".txt");
		metadata = import_iq3_metadata(info.get_metadata_file_path());
		IJ.run(imp, "Properties...", "channels=" + str(int(metadata['n_channels'])) + 
										" slices=" + str(int(metadata['z_pixels'])) + 
										" frames=1 unit=" + str(metadata['x_unit']) + 
										" pixel_width=" + str(metadata['x_physical_size']) + 
										" pixel_height=" + str(metadata['y_physical_size']) + 
										" voxel_depth=" + str(metadata['z_extent']/metadata['z_pixels']));
		info.set_xy_pixel_size_um(metadata['x_physical_size']);
		info.set_z_plane_spacing_um(metadata['z_extent']/metadata['z_pixels']);
		info = parse_info_from_filename(info);

	imp.show();
	imp.setC(info.get_mosaic_labeled_ch());
	IJ.run("Enhance Contrast", "saturated=0.35");
	zcrop_imp, z_lims = z_crop(imp);
	info.set_z_crop_frames(z_lims);
	info.set_mosaic_labeled_ch(imp.getC());
	IJ.setTool("rect");
	zcrop_imp.hide();
	imp2 = ZProjector.run(zcrop_imp, "max");
	imp2.show();
	imp2.setC(info.get_mosaic_labeled_ch());
	crop_roi = None;
	while crop_roi is None:
		WaitForUserDialog("Select XY region to crop...").show();
		crop_roi = imp2.getRoi();
	info.set_xy_crop_rect(crop_roi.getBounds().toString());
	imp2.close();
	zcrop_imp.show();
	zcrop_imp.setRoi(crop_roi);
	IJ.run("Crop", "");
	return imp, zcrop_imp, info, info.get_input_file_path();
		
def main():
	info = PrescreenInfo();
	default_path = "D:\\data\\Marcksl1 cell shape analysis\\zstacks"
	output_folder = output_folder_chooser(default_path);
	answer = 1;
	imp = None;
	while answer > 0:
		answer = 0;
		imp, zcrop_imp, info, default_path = perform_cropping(imp, info, output_folder, default_path);
		answer = save_etc(zcrop_imp, info, output_folder);
		if answer == 2:
			imp = None;
		if zcrop_imp is not None:
			zcrop_imp.changes = False;
			zcrop_imp.close();
		if imp is not None:
			imp.changes = False;
			imp.close();	
	return;
	
# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main();