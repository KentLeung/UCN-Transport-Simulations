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

if len(sys.argv) == 2:
    regionfile = sys.argv[1]
else:
    regionfile = sys.argv[1]
    eventsfile = sys.argv[2]
    trajfrac = float(sys.argv[3])
    cleanorerrors = sys.argv[4]


def visualize(regionfile, eventsfile, trajyn, trajfrac):#run the visualizer
    vs.scene.width =1200
    vs.scene.height = 800
    vs.scene.autocenter = vs.true
    
    regions = ucni.read_regions(regionfile);
    pieces = ucnr.draw_geom(regions);
    
    pieces[-1].color=vs.color.blue
    
    if len(sys.argv) == 5:
        eventstotrajectories(eventsfile)
        hitplaces = ucni.read_simfile('traj.txt');
        trails = ucnr.draw_simfile(hitplaces,min(trajfrac,max(trajfrac,0)));
#       if cleanorerrors == 'errors': draw error points
#           errorpoints = ucnr.draw_errorpoints()
def eventstotrajectories(filename): #Convert events.sim file to an acceptable trajectory file
    
    eventsfile = open(filename,'r');
    rawtext = eventsfile.read(); #extract text
    eventsfile.close();
    
    eventlines = rawtext.splitlines()[2:];
    stringevents = [];
    errors = [];
    errorpoints = [];
    
    for lines in eventlines:
        stringevents.append([lines.split()[0],lines.split()[1],lines.split()[5],lines.split()[6],lines.split()[7]]);
    
    #Convert read text values to floats
    events = [map(float,x) for x in stringevents]

    #log error particles
    for line in events:
        if line[1] == 19:
            errors.append(line[0])
            errorpoints.append([line[2], line[3], line[4]])

    #remove event code
    for event in events:
        event = event.pop(1)

    print 'There are {:d} errors.'.format(len(errors))

    finalevents = [];

    if cleanorerrors == 'clean':
        for i in events:
            if i[0] not in errors:
                finalevents.append(i)

    if cleanorerrors == 'errors':
        for i in events:
            if i[0] in errors:
                finalevents.append(i)

    savetxt('traj.txt', finalevents, fmt='%i,%f,%f,%f')
    if not (len(errorpoints) == 0):
        savetxt('errorpoints.txt', errorpoints, fmt='%f,%f,%f')

visualize(regionfile, eventsfile, trajyn, trajfrac)
