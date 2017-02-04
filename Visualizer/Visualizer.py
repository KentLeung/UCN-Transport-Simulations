from __future__ import division
import visual as vs
import os
import ucn_import_lib as ucni
import ucn_render_lib as ucnr
from numpy import savetxt
import sys

regionfile = ""
eventsfile = ""
trajyn = 0
trajfrac = 1

print sys.argv

if len(sys.argv) == 2:
    regionfile = sys.argv[1]
else:
    regionfile = sys.argv[1]
    eventsfile = sys.argv[2]
    trajyn = int(sys.argv[3])
    trajfrac = float(sys.argv[4])

print trajyn
print trajfrac


def visualize(regionfile, eventsfile, trajyn, trajfrac):#run the visualizer
    vs.scene.width =1200
    vs.scene.height = 800
    vs.scene.autocenter = vs.true
    
    regions = ucni.read_regions(regionfile);
    pieces = ucnr.draw_geom(regions);
    
    pieces[-1].color=vs.color.blue
    
    if trajyn == 1:
        eventstotrajectories(eventsfile)
        hitplaces = ucni.read_simfile('cleantraj.txt');
        os.remove("cleantraj.txt")
        trails = ucnr.draw_simfile(hitplaces,min(trajfrac,max(trajfrac,0)));

def eventstotrajectories(filename): #Convert events.sim file to an acceptable trajectory file
    
    eventsfile = open(filename,'r');
    rawtext = eventsfile.read(); #extract text
    eventsfile.close();
    
    eventlines = rawtext.splitlines()[2:];
    stringevents = [];
    
    for lines in eventlines:
        stringevents.append([lines.split()[0],lines.split()[5],lines.split()[6],lines.split()[7]]);
    
    #Convert read text values to floats
    events = [map(float,x) for x in stringevents]

    savetxt('cleantraj.txt', events, fmt='%i,%f,%f,%f')

visualize(regionfile, eventsfile, trajyn, trajfrac)
