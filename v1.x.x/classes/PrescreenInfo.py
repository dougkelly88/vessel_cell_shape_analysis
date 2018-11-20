import json, os, platform

class PrescreenInfo(object):
	"""class containing information to link prescreened, cropped cell shape data with original dataset"""
	_version_string = "1.0.0";
	_vessel_types = {'xISV', 'aISV', 'vISV', 'DLAV', 'DA', 'PCV'};
	_persist_filename = "IJ_marksl1_analysis_last_info.json";
	_persist_folder = "IJ_marksl1_analysis";
	def __init__(self, 
					input_file_path=None,
					metadata_file_path=None,
					z_crop_frames=None, 
					xy_crop_rect=None, 
					vessel_type='xISV',
					ch1_label=None,
					ch2_label=None, 
					xy_pixel_size_um=None,
					z_plane_spacing_um=None, 
					embryo_id=None, 
					experiment_id=None, 
					mosaic_labeled_ch=2):
		"""constructor"""
		self.input_file_path = input_file_path;
		self.metadata_file_path = metadata_file_path;
		self.z_crop_frames = z_crop_frames;
		self.xy_crop_rect = xy_crop_rect;
		self.vessel_type = vessel_type;
		self.ch1_label = ch1_label;
		self.ch2_label = ch2_label;
		self.xy_pixel_size_um = xy_pixel_size_um;
		self.z_plane_spacing_um = z_plane_spacing_um;
		self.embryo_id = embryo_id;
		self.experiment_id = experiment_id;
		self.mosaic_labeled_ch = mosaic_labeled_ch;
		self.__prescreeninfo__ = True;

	# Setters
	def set_input_file_path(self, input_file_path):
		if os.path.isfile(input_file_path):
			self.input_file_path = input_file_path;
		else:
			raise IOError("File " + input_file_path + " not found!");

	def set_metadata_file_path(self, metadata_file_path):
		if os.path.isfile(metadata_file_path):
			self.metadata_file_path = metadata_file_path;
		else:
			raise IOError("File " + metadata_file_path + " not found!");

	def set_z_crop_frames(self, z_crop_frames):
		if len(z_crop_frames) == 2:
			self.z_crop_frames = z_crop_frames;
		else:
			raise ValueError("Crop frames should be defined by a length 2 list: [start, stop]");

	def set_xy_crop_rect(self, xy_crop_rect):
		self.xy_crop_rect = xy_crop_rect;

	def set_vessel_type(self, vessel_type):
		if vessel_type in self._vessel_types:
			self.vessel_type = vessel_type;
		else:
			raise ValueError("Requested vessel type not supported!")

	def set_ch1_label(self, ch1_label):
		try:
			self.ch1_label = str(ch1_label);
		except:
			raise ValueError("Can't convert requested label to string");

	def set_ch2_label(self, ch2_label):
		try:
		   self.ch2_label = str(ch2_label);
		except: 
			raise ValueError("Can't convert requested label to string");

	def set_xy_pixel_size_um(self, xy_pixel_size_um):
		try:
			self.xy_pixel_size_um = float(xy_pixel_size_um);
		except: 
			raise ValueError("Can't convert request pixel size " + str(xy_pixel_size_um) + "to float");

	def set_z_plane_spacing_um(self, z_plane_spacing_um):
		try:
			self.z_plane_spacing_um = float(z_plane_spacing_um);
		except: 
			raise ValueError("Can't convert request pixel size " + str(z_plane_spacing_um) + "to float");

	def set_embryo_id(self, embryo_id):
		self.embryo_id = embryo_id;

	def set_experiment_id(self, experiment_id):
		self.experiment_id = experiment_id;

	def set_mosaic_labeled_ch(self, mosaic_labeled_ch):
		self.mosaic_labeled_ch = mosaic_labeled_ch;

	# Getters
	def get_input_file_path(self):
		return self.input_file_path;

	def get_metadata_file_path(self):
		return self.metadata_file_path;

	def get_z_crop_frames(self):
		return self.z_crop_frames;

	def get_xy_crop_rect(self):
		return self.xy_crop_rect;

	def get_vessel_type(self):
		return self.vessel_type;

	def get_ch1_label(self):
		return self.ch1_label;

	def get_ch2_label(self):
		return self.ch2_label;

	def get_xy_pixel_size_um(self):
		return self.xy_pixel_size_um;

	def get_z_plane_spacing_um(self):
		return self.z_plane_spacing_um;

	def get_embryo_id(self):
		ret = self.embryo_id if self.embryo_id else "undefined_embryo";
		return ret;

	def get_experiment_id(self):
		return self.experiment_id;

	def get_vessel_types(self):
		return list(PrescreenInfo._vessel_types);

	def get_mosaic_labeled_ch(self):
		return self.mosaic_labeled_ch;

	# persistence functions
	def save_info_to_json(self, file_path):
		try:
			f = open(file_path, 'w');
			json.dump(self.__dict__, f);
		finally:
			f.close();
		return;

	def load_info_from_json(self, file_path):
		try:
			f = open(file_path, 'r');
			dct = json.loads(f.read());
			if "__prescreeninfo__" in dct:
				self.set_input_file_path(dct["input_file_path"]);
				self.set_z_crop_frames(dct["z_crop_frames"]);
				self.set_xy_crop_rect(dct["xy_crop_rect"]);
				self.set_vessel_type(dct["vessel_type"]);
				self.set_ch1_label(dct["ch1_label"]);
				self.set_ch2_label(dct["ch2_label"]);
				self.set_xy_pixel_size_um(dct["xy_pixel_size_um"]);
				self.set_z_plane_spacing_um(dct["z_plane_spacing_um"]);
				self.set_embryo_id(dct["embryo_id"]);
				self.set_experiment_id(dct["experiment_id"]);
				self.set_metadata_file_path(dct["metadata_file_path"]);
				self.set_mosaic_labeled_ch(dct["mosaic_labeled_ch"]);
			else:
				raise ValueError("JSON file doesn't translate to prescreeninfo format")
		except IOError:
			print("IOError reading from JSON file");
			return False;
		except: 
			return False;
		finally:
			f.close();
		return True;

	def load_last(self):
		success = True;
		try:
			temp_path = self.get_persistence_file_location();
			if temp_path:
				temp_params_path = os.path.join(temp_path, PrescreenInfo._persist_filename);
				if os.path.isfile(temp_params_path):
					success = self.load_info_from_json(temp_params_path);
				else:
					success = False;
			else:
				success = False;
		except Exception as e:
			print("Warning: Error loading previous settings, reverting to default...");
			raise e;
			return False;
		if not success:
			print("Warning: Error loading previous settings, reverting to default...");
		return success;

	def persist(self):
		temp_path = self.get_persistence_file_location();
		if temp_path:
			temp_params_path = os.path.join(temp_path, PrescreenInfo._persist_filename);
			self.save_info_to_json(temp_params_path);
		return;

	def get_persistence_file_location(self):
		try:
			st = platform.mac_ver()[0];
			if not st:
				# windows
				temp_path = os.path.join(os.getenv('APPDATA'), PrescreenInfo._persist_folder);
			else:
				# mac
				temp_path = os.path.join(os.path.expanduser("~"), "Library", PrescreenInfo._persist_folder);
			if not os.path.isdir(temp_path):
				os.mkdir(temp_path);
		except Exception as e:
			print("Error: " + e.message);
			return "";
		return temp_path;
