from pylab import *
from sys import *

def narrate(story):
    print "\nNeutron #{:d}:\n".format(story[0][0])
    for event in story: describe(event)

def describe(event):

    reg = event[2]
    if event[3] == -1: xcode = reg
    else: xcode = event[3]
    
    print "         - t = {:f}".format(event[4])
    
    if event[1] == -1:
        print "   - Error in Reg. {:d} interacting with Reg. {:d}.".format(reg,xcode)
    elif event[1] == 0:
        print "      - r = ({:f},{:f},{:f}), v = ({:f},{:f},{:f})".format(event[5],event[6],event[7],event[8],event[9],event[10])
    elif event[1] == 1:
        print "\n   - Intersection with Reg. {:d}.".format(event[2])
    elif event[1] == 2:
        print "\n   - Intersection with cutplane {:d} in Reg. {:d}.".format(xcode,reg)
    elif event[1] == 3:
        print "   - Left specular bounce.\n"
    elif event[1] == 4:
        print "   - Left non-specular bounce.\n"
    elif event[1] == 5:
        print "   - Left specular bounce, depolarized.\n"
    elif event[1] == 6:
        print "   - Left non-specular bounce, depolarized.\n"
    elif event[1] == 7:
        print "   - Lost to wall due to non-penetration effects in Reg. {:d} interacting with Reg. {:d}.".format(reg,xcode)
    elif event[1] == 8:
        print "   - Absorbed in Reg. {:d}".format(reg)
    elif event[1] == 9:
        print "   - Absorbed and detected in Reg. {:d}".format(reg)
    elif event[1] == 10:
        print "   - Spin-flip."
    elif event[1] == 11:
        print "   - Beta decay."
    elif event[1] == 12:
        print "      - Checking for bounce or penetration."
    elif event[1] == 13:
        print "   - Penetrates wall or Fermi gradient.\n"
    elif event[1] == 14:
        print "   - Particle created with E = {:f} neV.".format(event[14])
    elif event[1] == 15:
        print "   - Energy shifted from Fermi gradient."
    elif event[1] == 16:
        print "   - Interaction with T-junction."
    elif event[1] == 17:
        print "   - Scattered in bulk."
    elif event[1] == 18:
        print "   - Detected by cutplane detector."
    elif event[1] == 19:
        print "   - Bounced off of cutplate {:d}.".format(xcode)
    elif event[1] == 20:
        print "   - Penetrated cutplate {:d}.".format(xcode)
    else: print "\n   - Unknown event code: {}\n".format(event[1])
    
with open(argv[1],'r') as file:
    raw = file.read().splitlines()

N = int(raw[0].split()[3])

events = [[] for i in range(N)]

for line in raw[2:]:
    temp = [int(x) for x in line.split()[:4]]+[float(x) for x in line.split()[4:]]
    events[temp[0]].append(temp)

for event in events: narrate(event)
print "\n"
