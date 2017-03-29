#Plotter for particle energy spectra

from pylab import *
import sys

nev2J = 1.602177e-28
vcutoff = 8/nev2J

#polynomial fit function
def polynomial(x,y,n,domain):
    p = polyfit(x,y,n)
    f = 0
    for i in p:
        f = f + p[i]*domain**(n-i)
        print f
    return f


#grab events.sim file
file = open(sys.argv[1],"r")
raw = file.read().splitlines()[2:]
file.close()

#create array with all events
events = []
for line in raw:
    events.append([float(x) for x in line.split()])

#separate particle instantiation events
poofs = []
for event in events:
    if event[1] == 14:
        poofs.append([int(event[0]),event[14],event[15]])

#separate angular detector cutplane events
hits = []
for event in events:
    if event[1] == 18:
        hits.append([int(event[0]),int(event[3]),event[14],event[15]])

#find all cutplanes that appear in hits
detectors  = list({item[1] for item in hits})
detid = {}
for i in range(0,len(detectors)): #dictionary to hold cp number vs cphits index
    detid[detectors[i]] = i

cphits = [[] for item in detid]
for i in range(0,len(detectors)):
    for hit in hits:
        if hit[1] == detectors[i]:
            cphits[i].append(hit)

close()
#plot starting distribution
if sys.argv[2] == 'start':
    hist([item[1] for item in poofs],25)
    xlabel("Energy (neV)")
    ylabel("Frequency")
    xlim((0,500))
    title("Initial Energy Distribution")

#plot distribution at a cutplane
elif int(sys.argv[2]) in detectors:
    hist([item[2] for item in cphits[detid[int(sys.argv[2])]]],25)
    xlabel("Energy (neV)")
    ylabel("Frequency")
    title("Energy Distribution at Cut Plane {:d}".format(int(sys.argv[2])))
    xlim((0,500))

show()
