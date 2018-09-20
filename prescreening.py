# perform first pass screening of maximum projections to populate analysis metadata
# D.J. Kelly, 2018-09-19, douglas.kelly@riken.jp

# python (jython) imports
import os, sys

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing
import javax.swing.table.TableModel

# imagej imports
from ij.gui import NonBlockingGenericDialog
from ij.io import DirectoryChooser

# custom imports
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_path)

def quickSetGridBagConstraints(gridbagconstraints, gridx, gridy, weightx, weighty):
	gridbagconstraints.weightx = weightx;
	gridbagconstraints.weighty = weighty;
	gridbagconstraints.gridx = gridx;
	gridbagconstraints.gridy = gridy;
	return gridbagconstraints;

def generateUI(file_list):
	# window
	frame = swing.JFrame("Marcksl1 prescreen tool")
	frame.setDefaultCloseOperation(swing.JFrame.DISPOSE_ON_CLOSE)
	frame.setPreferredSize(Dimension(700, 500))
	frame.setSize(Dimension(700, 500))
	frame.setLocationByPlatform(True)
	
	layout = GridBagLayout()
	frame.setLayout(layout)

	# image list
	lst_constraints = quickSetGridBagConstraints(GridBagConstraints(), 0, 0, 0.5, 0.33);
	lst_constraints.fill = GridBagConstraints.BOTH;
	lst = swing.JList(file_list)
	lst_scrollpane = swing.JScrollPane(lst)
	layout.setConstraints(lst_scrollpane, lst_constraints)
	frame.add(lst_scrollpane)
	
	# radio buttons
	radio_panel_constraints = quickSetGridBagConstraints(GridBagConstraints(), 1, 0, 0.25, 0.33)
	radio_panel_constraints.fill = GridBagConstraints.HORIZONTAL
	
	isvButton = swing.JRadioButton("ISV", True)
	dlavButton = swing.JRadioButton("DLAV")

	radio_group = swing.ButtonGroup()
	radio_group.add(isvButton)
	radio_group.add(dlavButton)
	
	radio_panel = swing.JPanel()
	radio_panel_layout = GridLayout(2, 1)
	radio_panel.setLayout(radio_panel_layout)
	radio_panel.add(isvButton)
	radio_panel.add(dlavButton)
	radio_panel.setBorder(swing.BorderFactory.createTitledBorder(
		swing.BorderFactory.createEtchedBorder(), "Vessel type"))
	layout.setConstraints(radio_panel, radio_panel_constraints)
	frame.add(radio_panel)
	
	# checkboxes
	chk_constraints = quickSetGridBagConstraints(GridBagConstraints(), 2, 0, 0.25, 0.33);
	chk_constraints.fill = GridBagConstraints.HORIZONTAL;
	chkbox = swing.JCheckBox("Is lumenised?", True);
	layout.setConstraints(chkbox, chk_constraints);
	frame.add(chkbox);

	# add ROI button
	roi_but_constraints = quickSetGridBagConstraints(GridBagConstraints(), 0, 1, 1, 0.16);
	roi_but_constraints.gridwidth = 3;
	roi_but_constraints.fill = GridBagConstraints.BOTH;
	roi_but = swing.JButton("Add current ROI")
	layout.setConstraints(roi_but, roi_but_constraints);
	frame.add(roi_but);

    # metadata display table
	mdtbl_constraints = quickSetGridBagConstraints(GridBagConstraints(), 0, 2, 1, 0.33);
	mdtbl_constraints.gridwidth = 3;
	mdtbl_constraints.fill = GridBagConstraints.BOTH;

	mdtbl_colnames = ['Experiment', 'Embryo', 'Imaged region index', 'ROI index', 'Vessel type', 'Lumenised?'];
	mdtbl_model = swing.table.DefaultTableModel(mdtbl_colnames, 10);
	mdtbl = swing.JTable();
	
	mdtbl_scrollpane = swing.JScrollPane(mdtbl);
	layout.setConstraints(mdtbl_scrollpane, mdtbl_constraints);
	frame.add(mdtbl_scrollpane);

	# buttons
	del_but_constraints = quickSetGridBagConstraints(GridBagConstraints(), 1, 3, 0.25, 0.16);
	del_but_constraints.fill = GridBagConstraints.BOTH;
	del_but = swing.JButton("Delete selected ROI");
	layout.setConstraints(del_but, del_but_constraints);
	frame.add(del_but);

	save_but_constraints = quickSetGridBagConstraints(GridBagConstraints(), 2, 3, 0.25, 0.16);
	save_but_constraints.fill = GridBagConstraints.BOTH;
	save_but = swing.JButton("Save...");
	layout.setConstraints(save_but, save_but_constraints);
	frame.add(save_but);

	# show window
	frame.setVisible(True)
	frame.pack()
	

def getFileList():
	default_root_path = r'D:\\data\Marcksl1 cell shape analysis\\MaximumProjections'
	file_list = []
	directory_list = []
	dc = DirectoryChooser('Select the root folder containing all data...')
	root_path = dc.getDirectory()
	if root_path is None:
		return file_list
	for directory, dir_names, file_names in os.walk(root_path):
		if (directory is not root_path):
			#file_list = [file_list.append(os.path.join(directory.replace(root_path, ''), file_name)) for file_name in file_names]
			file_list = [os.path.join(directory.replace(root_path, ''), file_name) for file_name in file_names if '.tif' in file_name]
	return root_path, file_list
	
def main():
	#print (sys.version_info) # debug
	#print(sys.path) # debug
	root_path, file_list = getFileList()
	generateUI(file_list)
	# get list of images by walking over root directory and populate FOV list
	# on selection of an image, load and display
	# on draw of ROI, add ROI and populate table of ROIs with defaults
	# on click of table, load that FOV, overlay existing ROIs, highlight selected ROI and allow modification/deletion
	# 
	

# It's best practice to create a function that contains the code that is executed when running the script.
# This enables us to stop the script by just calling return.
if __name__ in ['__builtin__','__main__']:
    main()