#This file varies the density of a regionfile with layers
#python SD2ConditionSpread.py regionfile densityfraction

from pylab import *
from sys import *
from UCNToolsLib import *

#inputs
regionfile = argv[1]
densityfraction = float(argv[2])
numberoflayers = 10.
totalthickness = 2.e-3
vacuumthickness = totalthickness*(1-densityfraction)/numberoflayers
SD2thickness = totalthickness*(densityfraction)/numberoflayers

#write parameter strings
vacuum = '0.17,{:.2e}'.format(vacuumthickness)
SD2 = '0.17,{:.2e}'.format(SD2thickness)

#vacuum layers
editRegion(regionfile,range(2,int(numberoflayers)*2+1,2),'Dim',vacuum)

#sd2 layers
editRegion(regionfile,range(3,int(numberoflayers)*2+2,2),'Dim',SD2)
