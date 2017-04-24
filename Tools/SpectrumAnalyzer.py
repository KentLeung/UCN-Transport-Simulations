#Plotter for particle energy spectra
#to run: python SpectrumAnalyzer.py events.sim

from pylab import *
import sys

nev2J = 1.602177e-28
nMASS = 1.67493e-27
vcutoff = (0.5)*nMASS*(8**2)/nev2J

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
        poofs.append([int(event[0]),int(event[3]),event[14],event[15]])

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

n1, binsinitial, patches = hist([i[3] for i in poofs],50,normed='true')
clf()
bins1 = []
for i in range(0,len(binsinitial)-1):
    bins1.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])

n2, binsinitial, patches = hist([i[3] for i in cphits[detid[3]]],50,normed='true')
clf()
bins2 = []
for i in range(0,len(binsinitial)-1):
    bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])

plot(bins1,n1*100,'k--',bins2,n2*100,'k-')
show()


#ylim((0,2))
#yticks([0,1,2])
#xlim((0,500))
#xticks([0,100,200,300,400,500])

#def plotinit(typeofgraph, destination):
#    close()
#    type = {"energy":1 , "direction":2}
#    lims = {"energy":0 , "direction":-100}
#    titles = {"energy":"Energy", "direction":"Vx"}
#    hist([item[type[typeofgraph]] for item in poofs],25)
#    xlabel("Energy (neV)")
#    ylabel("Frequency")
#    xlim((lims[typeofgraph],500))
#    title("Initial {:s} Distribution".format(titles[typeofgraph]))
#    savefig(destination + "/Initial {:s} Distribution".format(titles[typeofgraph]))
#    close()
#
#def plotcut(cutplane, typeofgraph, destination):
#    close()
#    type = {"energy":2 , "direction":3}
#    lims = {"energy":0 , "direction":-100}
#    titles = {"energy":"Energy", "direction":"Vx"}
#    hist([item[type[typeofgraph]] for item in cphits[detid[cutplane]]],25)
#    xlabel("Energy (neV)")
#    ylabel("Frequency")
#    title("{:s} Distribution at Cut Plane {:d}".format(titles[typeofgraph], cutplane))
#    xlim((lims[typeofgraph],500))
#    savefig(destination + "/{:d} {:s}".format(cutplane, titles[typeofgraph]))
#    close()

#destination = sys.argv[1][:-11]
#
#plotinit("energy",destination)
#plotinit("direction", destination)
#
#for i in detectors:
#    plotcut(i,"energy",destination)
#    plotcut(i,"direction", destination)


#close()


#plot starting distribution
#if sys.argv[2] == 'start':
#    type = {"energy":1 , "direction":2}
#    lims = {"energy":0 , "direction":-100}
#    titles = {"energy":"Energy", "direction":"Vx"}
#    hist([item[type[sys.argv[3]]] for item in poofs],25)
#    xlabel("Energy (neV)")
#    ylabel("Frequency")
#    xlim((lims[sys.argv[3]],500))
#    title("Initial {:s} Distribution".format(titles[sys.argv[3]]))
#    savefig(sys.argv[4] + "/Initial {:s} Distribution".format(titles[sys.argv[3]]))

#plot distribution at a cutplane
#elif int(sys.argv[2]) in detectors:
#    type = {"energy":2 , "direction":3}
#    lims = {"energy":0 , "direction":-100}
#    titles = {"energy":"Energy", "direction":"Vx"}
#    hist([item[type[sys.argv[3]]] for item in cphits[detid[int(sys.argv[2])]]],25)
#    xlabel("Energy (neV)")
#    ylabel("Frequency")
#    title("{:s} Distribution at Cut Plane {:d}".format(titles[sys.argv[3]], int(sys.argv[2])))
#    xlim((lims[sys.argv[3]],500))

#show()
