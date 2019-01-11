# module containing functions to handle loading/saving data
#
# D. J. Kelly, 2018-12-05, douglas.kelly@riken.jp

import os, re
from datetime import datetime
from ij import IJ
from ij import WindowManager as WM
from ij.io import OpenDialog, DirectoryChooser, FileSaver
from ij.gui import GenericDialog, YesNoCancelDialog
from loci.plugins import BF as bf
from loci.formats import ImageReader, MetadataTools
from loci.formats.gui import BufferedImageReader
from loci.plugins.in import ImporterOptions, ThumbLoader
from java.awt import Panel, Dimension, Checkbox, CheckboxGroup
from javax.swing import Box


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
		return "stop";
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
			return "stop";
		if continueDialog.yesPressed():
			return "continue_newimage";
		return "continue_sameimage";
	return "stop";

def get_image(info, default_path, used_files):
	"""handle getting image from file"""
	check_for_file_overlap_OK = False;
	while not check_for_file_overlap_OK:
		info.set_input_file_path(file_location_chooser(default_path));
		if info.get_input_file_path() in used_files:
			dlg = GenericDialog("Image already used...");
			dlg.addMessage("This image has already been used in this analysis run.  \n " + 
					"Continue with this image, or choose a new one?");
			dlg.setOKLabel("Continue");
			dlg.setCancelLabel("Choose again...");
			dlg.showDialog();
			check_for_file_overlap_OK = False if dlg.wasCanceled() else True;
		else:
			check_for_file_overlap_OK = True;
	import_opts, info = choose_series(info.get_input_file_path(), info)

	imps = bf.openImagePlus(import_opts);
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
	return imp, info;

def choose_series(filepath, info):
	"""if input file contains more than one image series (xy position), prompt user to choose which one to use"""
	# todo: if necessary (e.g. if lots of series), can improve thumbnail visuals based loosely on https://github.com/ome/bio-formats-imagej/blob/master/src/main/java/loci/plugins/in/SeriesDialog.java
	import_opts = ImporterOptions();
	import_opts.setId(filepath);
	
	reader = ImageReader();
	ome_meta = MetadataTools.createOMEXMLMetadata();
	reader.setMetadataStore(ome_meta);
	reader.setId(filepath);
	no_series = reader.getSeriesCount();
	if no_series == 1:
		return import_opts, info;
	else:
		series_names = [ome_meta.getImageName(idx) for idx in range(no_series)];
		dialog = GenericDialog("Select series to load...");
		dialog.addMessage("There are multiple series in this file! \n" + 
						"This is probably because there are multiple XY stage positions. \n " + 
						"Please choose which series to load: ");
		thumbreader = BufferedImageReader(reader);
		cbg = CheckboxGroup();
		for idx in range(no_series):
			p = Panel();
			p.add(Box.createRigidArea(Dimension(thumbreader.getThumbSizeX(), thumbreader.getThumbSizeY())));
			ThumbLoader.loadThumb(thumbreader, idx, p, True);
			dialog.addPanel(p);
			cb = Checkbox(series_names[idx], cbg, idx==0);
			p.add(cb);

		dialog.showDialog();
		if dialog.wasCanceled():
			raise KeyboardInterrupt("Run canceled");
		if dialog.wasOKed():
			selected_item = cbg.getSelectedCheckbox().getLabel();
			selected_index = series_names.index(selected_item);
			info.set_input_file_series(selected_index);
			for idx in range(0, no_series):
				import_opts.setSeriesOn(idx, True) if (idx==selected_index) else import_opts.setSeriesOn(idx, False);
	reader.close();
	return import_opts, info

def import_iq3_metadata(metadata_path):
	"""import basic image metadata based on the metadata saved by iQ3 software at acquisition time"""
	x_fmt_str = 'x \: (?P<x_pixels>\d+\.?\d*) \* (?P<x_physical_size>\d+\.?\d*) \: (?P<x_unit>\w+)';
	y_fmt_str = 'y \: (?P<y_pixels>\d+\.?\d*) \* (?P<y_physical_size>\d+\.?\d*) \: (?P<y_unit>\w+)';
	z_fmt_str = '\s*Repeat Z \- (?P<z_extent>[+-]?\d+\.?\d*) (?P<z_unit>\w+) in (?P<z_pixels>\d+\.?\d*) planes \(centre\)';
	t_fmt_str = '\s*Repeat T \- (?P<n_frames>\d+\.?\d*) times \((?P<frame_interval>\d+\.?\d*) (?P<time_unit>(min|sec))\)';
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