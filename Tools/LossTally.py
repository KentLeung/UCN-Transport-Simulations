#Loss Tally - Erik Lutz
#Use this code to find where losses happen in events.sim files
#to run: python LossTally.py events.sim

from pylab import *
from sys import *

# OPTIONS
detailed = False

# import data
file = open(argv[1],'r')
if detailed:
    raw = [map(lambda x: int(x),line.split()[:4]) + map(lambda x: float(x),line.split()[4:]) for line in file.read().splitlines()[2:]]
else:
    raw = [map(lambda x: int(x),line.split()[:4]) for line in file.read().splitlines()[2:]]
file.close()

Unphysical = []                 # -1
NonPenetrationWallLoss = []     # 7
RegionAbsorption = []           # 8
DetectorAbsorbtion = []         # 9
BetaDecay = []                  # 11
WallOrBulkPenetration = []      # 13
FermiGradientPenetration = []   # 20

codes = {   -1:Unphysical,
            7:NonPenetrationWallLoss,
            8:RegionAbsorption,
            9:DetectorAbsorbtion,
            11:BetaDecay,
            13:WallOrBulkPenetration,
            20:FermiGradientPenetration
        }

for line in raw:
    try: codes[line[1]].append(line)
    except: continue
print "\n###### Loss Data ######\n"
print "Unphysical Losses          =",len(Unphysical)
print "Non-Penetration Wall Loss  =",len(NonPenetrationWallLoss)
print "Region Absorption          =",len(RegionAbsorption)
print "Detector Absorption        =",len(DetectorAbsorbtion)
print "Beta Decay                 =",len(BetaDecay)
print "Wall or Bulk Penetration   =",len(WallOrBulkPenetration) - len(FermiGradientPenetration),"\n"

tot = len(Unphysical) + len(NonPenetrationWallLoss) + len(RegionAbsorption) + len(DetectorAbsorbtion) + len(BetaDecay) + len(WallOrBulkPenetration) - len(FermiGradientPenetration)
print "Total                      =",tot,"\n"

print "###### Excel Data ######\n"

print len(Unphysical)
print len(NonPenetrationWallLoss)
print len(RegionAbsorption)
print len(DetectorAbsorbtion)
print len(BetaDecay)
print len(WallOrBulkPenetration)- len(FermiGradientPenetration)
print tot,"\n"



#Event Codes:    -1 => Error
#                0 => Record a trajectory point
#                1 => Intersection with a region
#                2 => Intersection with a cut-plane
#                3 => Leaving a specular bounce
#                4 => Leaving a non-specular bounce
#                5 => Leaving a specular bounce where particle was depolarized
#                6 => Leaving a non-specular bounce where particle was depolarized
#                7 => Loss of particle upon wall interaction other than due to penetration.
#                8 => Absorption of particle in a region
#                9 => Absorption of particle in a detector region.
#                10 => Spin-flip
#                11 => Beta-decay of a neutron
#                12 => Entering bounce or checking for reflection off an increased Fermi potential.
#                13 => Particle lost due to penetration into wall or penetrates a bulk medium.
#                14 => A created particle's dynamical information.
#                15 => Particle gets an energy shift due to an interaction with a bulk medium.
#                16 => Intersection with a T-junction.
#                17 => Particle scattered inside a bulk medium.
#                18 => Detection of a particle by a cut-plane
#                19 => Particle bounced off a cutplane with defined surface roughness
#                20 => Particle penetrated a cutplane with defined surface roughness
#                21 => Particle penetrated a cutplane between 2 regions with differing potentials
