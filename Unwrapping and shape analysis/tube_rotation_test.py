# generate a tilted cylinder sliced through z - n.b. rotation around z-axis doesn't matter (here, because of symmetry, 
# IRL because unwrapping will occur around z

import math
from ij import IJ
from ij.gui import EllipseRoi

w = 100;
h = 100;
d = 100;
r = 20;
zang = 0;
yang = -30;
xang = -45;
rotation_centre = (50, 50, 50);

xsemiax = r / math.cos(math.radians(xang));
ysemiax = r / math.cos(math.radians(yang));

imp = IJ.createImage("tilted tube", w, h, d, 8);

for zidx in range(d):
	imp.setZ(zidx+1);
	cx = (zidx + 1 - rotation_centre[2]) * math.tan(math.radians(yang));
	cy = (zidx + 1 - rotation_centre[2]) * math.tan(math.radians(xang));
	
	aspect_ratio = 
	roi = EllipseRoi(cx - xsemiax, xy - ysemiax, cx + xsemiax, cy + ysemiax, xsemiax/ysemiax

