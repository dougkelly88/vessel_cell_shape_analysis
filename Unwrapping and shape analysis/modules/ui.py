# module containing functions to handle user interface
#
# D. J. Kelly, 2019-01-11, douglas.kelly@riken.jp

from ij.gui import NonBlockingGenericDialog

def choose_segmentation_and_projection_channels(info):
	dialog = NonBlockingGenericDialog("Select channels...");
	channels = [info.ch1_label, info.ch2_label];
	dialog.addRadioButtonGroup("Segmentation channel: ", 
							channels,
							1, len(channels), 
							channels[0]);
	dialog.addRadioButtonGroup("Projection channel: ", 
							channels,
							1, len(channels), 
							channels[1]);
	dialog.showDialog();
	if dialog.wasCanceled():
		return None;
	seg_ch = dialog.getNextRadioButton();
	proj_ch = dialog.getNextRadioButton();
	return channels.index(seg_ch), channels.index(proj_ch);