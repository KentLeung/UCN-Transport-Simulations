#This code reads the events.sim output file from the Transport Code and
#returns a trajectory text file where each line is of the form n,x,y,z

import os
from numpy import savetxt

eventsfile = open('events.sim','r');
rawtext = eventsfile.read(); #extract text
  
eventsfile.close();
  
eventlines = rawtext.splitlines()[2:];
stringevents = [];
  
for lines in eventlines:
    stringevents.append([lines.split()[0],lines.split()[5],lines.split()[6],lines.split()[7]]);
  
#Convert read text values to floats
events = [map(float,x) for x in stringevents]

savetxt('cleantraj.txt', events, fmt='%i,%f,%f,%f')
