#Bottling Effects - Erik Lutz
#Use this code to extract bottling data from events.sim files
#to run: python Bottling.py events.sim

from pylab import *
from sys import *

# OPTIONS
bottlingregion = 5

# import data
file = open(argv[1],'r')
raw = file.read().splitlines()[2:]
file.close()

bottlingdata = []
for line in raw:
    splitline = line.split()
    if int(splitline[2]) == bottlingregion:
        bottlingdata.append(splitline)
