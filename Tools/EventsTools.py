#Various tools for handling large events.sim files
#to run: python EventsTools.py events.sim [args]

from pylab import *
from sys import *

# Parameters
headerlength = "short"

# functions
def printEvent(line):
    if headerlength == "short":
        print "  {:11s}{:8s}{:7s}{:8s}{:12s}{:12s}{:12s}{:12s}".format(str(line[0]),str(line[1]),str(line[2]),str(line[3]),str(line[4]),str(line[5]),str(line[6]),str(line[7]))

# import data
file = open(argv[1],'r')
raw = [map(lambda x: int(x),line.split()[:4]) + map(lambda x: float(x),line.split()[4:]) for line in file.read().splitlines()[2:]]
file.close()

# resources
shortheader = "#Neutron# | ecode | Reg# | xcode |     t     |     x     |     y     |     z     |"
longheader  = "#Neutron# | ecode | Reg# | xcode |     t     |     x     |     y     |     z     |     vx    |     vy    |     vz    |   Spinz   |   Data1   |   Data2   |   Data3   |   Data4   |"

# print particle events
if argv[2] == "printparticle":
    if headerlength == "short": print shortheader
    i = 1
    while True:
        if raw[i][0] < int(argv[3]):
            i += 1
        elif raw[i][0] == int(argv[3]):
            printEvent(raw[i])
            i += 1
        else: break
