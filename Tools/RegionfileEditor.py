import sys

def format(line):
    a = ""
    a = a + "{:6s}".format(line[0])
    a = a + "{:7s}".format(line[1])
    a = a + "{:25s}".format(line[2])
    a = a + "{:16s}".format(line[3])
    a = a + "{:16s}".format(line[4])
    a = a + "{:14s}".format(line[5])
    a = a + "{:7s}".format(line[6])
    a = a + "{:8s}".format(line[7])
    a = a + "{:8s}".format(line[8])
    a = a + "{:11s}".format(line[9])
    a = a + "{:11s}".format(line[10])
    a = a + "{:7s}".format(line[11])
    a = a + "{:14s}".format(line[12])
    a = a + "{:5s}".format(line[13])
    a = a + "{:4s}".format(line[14])
    a = a + "{:4s}".format(line[15])
    a = a + "{:4s}".format(line[16])
    a = a + "{:1s}\n".format(line[17])
    return a



file = open(sys.argv[1], "r")
raw = file.read().splitlines()
file.close()

header = "Reg#  RType       BP(x,y,z)           Dim[m]          Orient          Grad B [T/m]  Spec   Loss    Depol   WPot[neV]  BPot[neV]  Scat   Absorb[1/s]   Det  PM  SM  LM  DM\n"
raw = [i.split() for i in raw[1:-1]]
options = {'Reg':0, 'RType':1, 'BP':2, 'Dim':3, 'Orient':4, 'Grad B':5, 'Spec':6, 'Loss':7, 'Depol':8, 'WPot':9, 'BPot':10, 'Scat':11, 'Absorb':12, 'Det':13, 'PM':14, 'SM':15, 'LM':16, 'DM':17}


#########################################################################################################
#Editor Logic#
if sys.argv[2] == 'num':
    for i in range(0, len(raw)):
        raw[i][0] = str(i)

if sys.argv[2] == 'param': #Change a specific parameter on several or all lines
    #get values to change
    regs  = sys.argv[3]
    param = sys.argv[4]
    value = sys.argv[5]

    #change all values in a column to the same thing
    if regs == "-1":
        for i in raw:
            i[options[param]] = value
    #Change a value in specified regions to something
    else:
        regs = regs.split(',')
        for i in raw:
            if i[0] in regs:
                i[options[param]] = value

#########################################################################################################


output = [header]
for i in raw:
    output.append(format(i))
output.append("/")

file = open(sys.argv[1], "w")
file.writelines(output)
file.close()
