#creates regionfiles and connexfiles of a certain number of layers

from pylab import *
from sys import *
from UCNToolsLib import *

#input
numoflayers = int(argv[2])
shape = argv[1]

#parameters (also go in and change ending variable to match if you make changes)
SD2BulkPot = '106.0'
WallPot = '350.0'
SD2Dim = '0.17,0.03900'
FrostSpec = '1.000'
WallSpec = '0.970'
SD2mfp = '0.04'
SD2absorb = '25.0'
LayerDim = '0.17,1.00e-5'
TravelDim = '0.17,0.30000'
startDim = '0.17,0.01000'
SD2Det = '2'
TransDet = '1'
FrostDet = '3'
roughnesscode = '2(1)'

#vars
lineTemplate = ['0', '2', '*', '0.0,0.0', '0.0,-90.0,0.0', '0.0,0.0,0.0', '0.0', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '0.0', '0', '0', '0', '0', '0']
header = "Reg#  RType  BP(x,y,z)           Dim[m]          Orient          Grad B [T/m]  Spec   Loss    Depol   WPot[neV]  BPot[neV]  Scat   Absorb[1/s]   Det  PM  SM  LM  DM\n"
options = {'Reg':0, 'RType':1, 'BP':2, 'Dim':3, 'Orient':4, 'Grad B':5, 'Spec':6, 'Loss':7, 'Depol':8, 'WPot':9, 'BPot':10, 'Scat':11, 'Absorb':12, 'Det':13, 'PM':14, 'SM':15, 'LM':16, 'DM':17}

#create file template
output = [header]
if shape == 'straight':
    ending = [['0', '2', '*', '0.17,0.10000', '0.0,-90.0,0.0', '0.0,0.0,0.0', '0.970', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '9999.', '1', '1', '0', '0', '0']]
    ending[0][0] = str(numoflayers*2+4)
if shape == 'bent':
    ending =[['0', '2', '*', '0.17,0.15000', '0.0,-60.0,0.0', '0.0,0.0,0.0', '0.970', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '0.0', '0', '0', '0', '0', '0'], \
             ['0', '2', '*', '0.17,0.15000', '0.0,-30.0,0.0', '0.0,0.0,0.0', '0.970', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '0.0', '0', '0', '0', '0', '0'], \
             ['0', '2', '*', '0.17,0.30000', '0.0,0.0,0.0', '0.0,0.0,0.0', '0.970', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '0.0', '0', '0', '0', '0', '0'], \
             ['0', '2', '*', '0.17,0.10000', '0.0,0.0,0.0', '0.0,0.0,0.0', '0.970', '2.0e-4', '1.0e-6', '350.0', '0.0', '0.0', '9999.', '1', '1', '0', '0', '0']]
    ending[0][0] = str(numoflayers*2+4)
    ending[1][0] = str(numoflayers*2+5)
    ending[2][0] = str(numoflayers*2+6)
    ending[3][0] = str(numoflayers*2+7)


for i in range(0,4+numoflayers*2):
    line = lineTemplate
    line[0] = str(i)
    output.append(formatRegion(lineTemplate))
for i in ending:
    output.append(formatRegion(i))
output.append("/")

#write out
Regfile = open("Regionfile", "w")
Regfile.writelines(output)
Regfile.close()

#####Begin editing lines#####
#0th region
editRegion("Regionfile",[0],'BP','0.0,0.0,0.0')
editRegion("Regionfile",[0],'Dim',startDim)
editRegion("Regionfile",[0],'Spec',FrostSpec)
editRegion("Regionfile",[0],'BPot',WallPot)
editRegion("Regionfile",[0],'Absorb','9999.')

#SD2 Crystal Region
editRegion("Regionfile",[1],'Dim',SD2Dim)
editRegion("Regionfile",[1],'Spec',FrostSpec)
editRegion("Regionfile",[1],'BPot',SD2BulkPot)
editRegion("Regionfile",[1],'Scat',SD2mfp)
editRegion("Regionfile",[1],'Absorb',SD2absorb)
editRegion("Regionfile",[1],'PM','4')
editRegion("Regionfile",[1],'Det',SD2Det)

#Guide Region
editRegion("Regionfile",[numoflayers*2+3],'Dim',TravelDim)
editRegion("Regionfile",[numoflayers*2+3],'Spec',WallSpec)

#vacuum layers
for i in range(0,numoflayers+1):
    editRegion("Regionfile",[2+2*i],'Dim',LayerDim)
    editRegion("Regionfile",[2+2*i],'Spec',FrostSpec)

#SD2 layers
for i in range(0,numoflayers):
    editRegion("Regionfile",[3+2*i],'Dim',LayerDim)
    editRegion("Regionfile",[3+2*i],'Spec',FrostSpec)
    editRegion("Regionfile",[3+2*i],'Absorb',SD2absorb)
    editRegion("Regionfile",[3+2*i],'Scat',SD2mfp)
    editRegion("Regionfile",[3+2*i],'Det',FrostDet)
    editRegion("Regionfile",[3+2*i],'PM','4')
    editRegion("Regionfile",[3+2*i],'BPot',SD2BulkPot)

#Connexfile stuff

if shape == 'straight':
    blankConnex("Connexfile",numoflayers*2 + 4)
    #last line
    editHandling("Connexfile",[numoflayers*2 + 4],'4')



if shape == 'bent':
    blankConnex("Connexfile",numoflayers*2 + 7)
    #last line
    editHandling("Connexfile",[numoflayers*2 + 7],'4')


#first line
editHandling("Connexfile",[0],'1')

#layers
for i in range(2,numoflayers*2+3):
    editHandling("Connexfile",[i],roughnesscode)


#Print info for spectral script use
if shape == 'straight':
    print "Detector at region {:d}.".format(numoflayers*2+4)
if shape == 'bent':
    print "Detector at region {:d}.".format(numoflayers*2+7)










