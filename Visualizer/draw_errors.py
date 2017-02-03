# drawerrors.py uses VPython to draw unphysically lost neutrons in sim geometry.
# Geometry is defined by A.T. Holley's sim standards: regions_* and conn*
# This exists to verify a text based file.
# ARR 1/11/12, 2/06/12

from __future__ import division
from visual import *
import os
from ucn_import_lib import *
from ucn_render_lib import *
#import Carbon.Snd

scene.width =1200
scene.height = 800
scene.autocenter = true
# # Locate scene in 2nd monitor. Hard-coded
# scene.x = 1440 
# scene.y = -457

axes = false

#os.chdir(os.getenv("HOME")+'/Dropbox/Development/Visualizer')


regions = read_regions();
connex = read_connex();
pieces = draw_geom(regions, connex);
# scene.autocenter = false

hitplaces = read_simfile('errortraj.txt');

problems = [];
#errcolors is a dictionary that matches an error code to a corresponding color
errcolors = {-1.0 : color.red,#intersection
             -2.0 : color.green, #propagation in bulk
             -3.0 : color.blue, #escaped geometry
             -4.0 : color.yellow, #determining if T-junction
             -5.0 : color.orange, # bounce without velocity into wall
             -6.0 : color.cyan, # error determining wall loss
             -7.0 : color.magenta, # error determining bounce type
             -8.0 : (1,0.7,0.2), #error determining polarization state
             -9.0 : (1,1,1), # error determining result of nonspecular bounce
             -10.0: (.2,0.7,0.2) # propagation through change in Fermi potential
}
errorballs = [];#last recorded location for lost neutrons gets a sphere.

for places in hitplaces:
  if places[4]<0: #all error codes are negative
    if len(problems)>1 and problems[-1][0]==places[0]: #double-counted...replace
      problems[-1]=[places[0],places[4]];
      errorballs[-1].pos = vector(places[1],places[2],places[3])
    else: #new element.
      problems.append([places[0],places[4]])
      errorballs.append(sphere(pos = vector(places[1],places[2],places[3]),
                               radius = .008, color=errcolors[places[4]]))
      errorballs[-1].material = materials.emissive;#make those balls glow
errortrails = [];
prev = -1;
for i in problems:
  if i[0]!=prev: #same neutron was flagged multiple times.
    prev = i[0];
    errortrails.append(curve(color=errcolors[i[1]]))
    for event in hitplaces:
      if event[0]==i[0]:
        errortrails[-1].append(pos=vector(event[1],event[2],event[3]));
del hitplaces
Carbon.Snd.SysBeep(1) #Done rendering- play a quick beep.
