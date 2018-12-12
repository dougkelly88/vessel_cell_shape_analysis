import ij.WindowManager as WM
titles = WM.getImageTitles();
for title in titles:
	WM.getImage(title).close();