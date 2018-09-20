# perform first pass screening of maximum projections to populate analysis metadata
# D.J. Kelly, 2018-09-19, douglas.kelly@riken.jp

# python (jython) imports
import os, sys

# java imports - aim to have UI components entirely swing, listeners and layouts awt
from java.awt import Dimension, GridBagLayout, GridBagConstraints, GridLayout
import javax.swing as swing

# imagej imports
from ij.gui import NonBlockingGenericDialog
from ij.io import DirectoryChooser

def generateUI(file_list):
	# window
	frame = swing.JFrame("ROI vessel tool")
	frame.setDefaultCloseOperation(swing.JFrame.DISPOSE_ON_CLOSE)
	frame.setPreferredSize(Dimension(500, 500))
	frame.setSize(Dimension(500, 500))
	frame.setLocationByPlatform(True)
	
	#dialog = NonBlockingGenericDialog("ROI vessel tool")
	#dialog.setPreferredSize(awt.Dimension(500, 500))
	layout = GridBagLayout()
	frame.setLayout(layout)

	# image list
	lst_constraints = GridBagConstraints()
	lst_constraints.fill = GridBagConstraints.BOTH
	lst_constraints.weightx = 0.5
	lst_constraints.weighty = 0.8
	lst_constraints.gridx = 0
	lst_constraints.gridy = 0
	lst = swing.JList(file_list)
	lst_scrollpane = swing.JScrollPane(lst)
	layout.setConstraints(lst_scrollpane, lst_constraints)
	frame.add(lst_scrollpane)
	
	
	#lst = awt.List(10, False)
	#[lst.add(item) for item in file_list]
	#layout.setConstraints(lst, lst_constraints)
	#dialog.add(lst)

	# radio buttons
	radio_panel_constraints = GridBagConstraints()
	radio_panel_constraints.weightx = 0.5
	radio_panel_constraints.weighty = 0.8
	radio_panel_constraints.gridx = 1
	radio_panel_constraints.gridy = 0
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
	#radio_frame = awt.ButtonGroup();
	#radio_frame_constraints = awt.GridBagConstraints()
	#radio_frame_constraints.weightx = 0.5
	#radio_frame_constraints.weighty = 0.8
	#radio_frame_constraints.gridx = 1
	#radio_frame_constraints.gridy = 0

	# buttons
	#buts = dialog.getButtons()
	#print(buts)

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