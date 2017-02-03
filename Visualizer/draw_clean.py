#  drawclean.py uses VPython to draw neutrons in a simulation geometry.
#  Geometry is defined by A.T.Holley.
#  ARR 1/11/12, 2/06/12, 2012-04-09

from __future__ import division
from visual import *
import os
from ucn_import_lib import *
from ucn_render_lib import *

scene.width =1200
scene.height = 800
scene.autocenter = true
# #Locate scene in 2nd monitor. Hard-coded, because it feels good to be dirty.
# scene.x = 1440 
# scene.y = -457

#os.chdir(os.getenv("HOME")+'/Dropbox/NCSU/Research/Simulation/AdamUCN/Visualizer')


regions = read_regions();
connex = read_connex();
pieces = draw_geom(regions, connex);

pieces[-1].color=color.blue

hitplaces = read_simfile('noscattercleantraj.txt');
trails = draw_simfile(hitplaces,1);
