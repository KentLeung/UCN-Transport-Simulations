#Analyzer for particle energy spectra, produces lists of energy
#to run: python SpectrumAnalyzer.py events.sim events.sim events.sim...

from pylab import *
import sys

def dataGrab(directory):
    #grab events.sim file
    file = open(directory,"r")
    raw = file.read().splitlines()[2:]
    file.close()

    #create array with all events
    events = []
    for line in raw:
        events.append([float(x) for x in line.split()])

    #separate particle instantiation events and angular detector cutplane events, log detectors
    poofs = []
    hits = []
    detectors = set()
    for event in events:
        if event[1] == 14:
            poofs.append([int(event[0]),int(event[3]),event[14],event[15]])
            # data = [Neutron #, xcode, energy, forward directed energy]
        elif event[1] == 18:
            hits.append([int(event[0]),int(event[3]),event[14],event[15]])
            detectors.add(event[3])
            #if event[14] <=80: print "N#{:d}: E = {:f}".format(int(event[0]),event[14])
    detectors = list(detectors)

    #sort cutplanes that appear in hits
    detid = {}
    for i in range(0,len(detectors)): #dictionary to hold cp number vs cphits index
        detid[detectors[i]] = i

    cphits = [[] for item in detid]
    for i in range(0,len(detectors)):
        for hit in hits:
            if hit[1] == detectors[i]:
                cphits[i].append(hit)
    return poofs, cphits, detid, detectors

def writefile(name,list,num):
    with open(name, 'w') as file:
        file.write("Neutron Energy (neV) - {:d}\n".format(num))
        for item in list:
            file.write(str(item[2]) + "\n")

for eventfile in sys.argv[1:]:
    data = dataGrab(eventfile)
    filename = "_".join(eventfile.replace("\ "," ").split('/')[-4:-1]) + "_"
    writefile(filename + "Poofs" + ".txt",data[0],len(data[0]))
    for cutplane in data[3]:
        writefile(filename + "CP" + str(int(cutplane)) + ".txt",data[1][data[2][cutplane]],len(data[0]))

