import os, shutil
from ij import IJ
from ij.plugin import FolderOpener, HyperStackConverter, Memory, Concatenator
from ij.io import FileSaver 

max_memory = Memory().maxMemory();
no_paths_to_skip = 0;
root = "D:\\data\\2019-02-22 Lumen stained cell shape analysis\\Control LynEGFP";
#walkout = [(root, subdir) for root,subdir,_ in os.walk(root)];
#paths_to_work_on=[w[0] for w in walkout if len(w[1])==0];
paths_to_work_on = ['e1\\e1b'];
paths_to_work_on = [os.path.join(root, p) for p in paths_to_work_on];
paths_with_problems = [];
for idx, p in enumerate(paths_to_work_on):
	# skip first path...
	if idx>(no_paths_to_skip-1):
		try:
			file_list = os.listdir(p);
			metadata_file = [f for f in file_list if os.path.splitext(f)[1]=='.txt'][0]
			fname = os.path.splitext(metadata_file)[0];
			print("Working on file {} of {}: {}".format(idx+1, len(paths_to_work_on), fname));
			# check memory requirements
			imfilenames = [f for f in file_list if os.path.splitext(f)[1]=='.tif'];
			imfilename = imfilenames[0];
			fsize = os.path.getsize(os.path.join(p, imfilename));
			if (len(file_list)-1) > max_memory//fsize:
				paths_with_problems.append(p);
				print("Too much data in path {}; moving to next path...".format(p));
				continue;
				
			print("Loading data...");
			if len(imfilenames)>2:
				fo = FolderOpener();
				imp = fo.open(p);
				imp.show();
				print("N_z = {}".format(imp.getNSlices()));
				print("N_c = {}".format(imp.getNChannels()));
				print("N_t = {}".format(imp.getNFrames()));
				
			else:
				ch1_imp = IJ.openImage(os.path.join(p, imfilenames[0]));
				ch2_imp = IJ.openImage(os.path.join(p, imfilenames[1]));
				imp = Concatenator().concatenate(ch1_imp, ch2_imp, False);
			imp2 = HyperStackConverter.toHyperStack(imp, 2, int(imp.getNSlices()/2), 1, "xyzct", "Composite");
			imp2.show();
			imp.close();
			print("Saving data...");
			out_dir = os.path.dirname(os.path.dirname(p))
			save_path = os.path.join(out_dir, fname + ".tif");
			FileSaver(imp2).saveAsTiffStack(save_path);
			print("Copying metadata...");
			shutil.copyfile(os.path.join(p, fname + ".txt"), os.path.join(out_dir, metadata_file));
			imp2.close();
			
		except Exception as e:
			print(e);
			print("Problem with {}".format(p));
			paths_with_problems.append(p);
print("Done");
if len(paths_with_problems)>0:
	print("Problems with:");
	for pwp in paths_with_problems:
		print(pwp);
