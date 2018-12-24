#Plotter for particle energy spectra
#to run: python SpectrumAnalyzer.py initialevents.sim (comparator.sim) events.sim events.sim ...

from pylab import *
import sys

nev2J = 1.602177e-28
nMASS = 1.67493e-27
vcutoff = (0.5)*nMASS*(8**2)/nev2J

colordict = {1:'b',2:'g',3:'c',4:'y',5:'m',6:'r'}
colorindex = 1

def dataGrab(directory):
    #grab events.sim file
    file = open(directory,"r")
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
    return poofs, cphits, detid

def plotDistBaseline((poofs,cphits,detid),plottype,nbins):
    if plottype == "energy":
        n1, binsinitial, patches = hist([i[2] for i in poofs],nbins,normed=normalized)
        clf()
        bins1 = []
        for i in range(0,len(binsinitial)-1):
            bins1.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        n2, binsinitial, patches = hist([i[2] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized)
        clf()
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins1,n1,linewidth=3,color='k',ls='-')
        plot(bins2,n2,linewidth=3,color='k',ls='--')

    if plottype == "angular":
        n1, binsinitial, patches = hist([i[3] for i in poofs],nbins,normed=normalized)
        clf()
        bins1 = []
        for i in range(0,len(binsinitial)-1):
            bins1.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        n2, binsinitial, patches = hist([i[3] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized)
        clf()
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins1,n1,linewidth=3,color='k',ls='-')
        plot(bins2,n2,linewidth=3,color='k',ls='--')

    if dataPrint == "ON":
        print "\nBaseline Initial"
        for i in range(0,len(bins1)):
            print bins1[i],n1[i]
        print "\nBaseline Final"
        for i in range(0,len(bins2)):
            print bins2[i],n2[i]

def plotDistComparator((poofs,cphits,detid),plottype,nbins):
    if plottype == "energy":
        n2, binsinitial, patches = hist([i[2] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized,alpha=0)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2,linewidth=3,color='r',ls=':')
    
    if plottype == "angular":
        n2, binsinitial, patches = hist([i[3] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized,alpha=0)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2,linewidth=3,color='r',ls=':')

    if dataPrint == "ON":
        for i in range(0,len(bins2)):
            print bins2[i],n2[i]


def plotDist((poofs,cphits,detid),plottype,nbins):

    global colorindex
    
    if plottype == "energy":
        n2, binsinitial, patches = hist([i[2] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized,alpha=0)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2,color=colordict[colorindex],ls='-')
        colorindex += 1
    
    if plottype == "angular":
        n2, binsinitial, patches = hist([i[3] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized,alpha=0)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2,color=colordict[colorindex],ls='-')
        colorindex += 1

    if dataPrint == "ON":
        for i in range(0,len(bins2)):
            print bins2[i],n2[i]


#Plot stuff
plottype = "energy"
nbins = 30
dataPrint = "ON"
normalized = False
comparator = False

close()

plotDistBaseline(dataGrab(sys.argv[1]),plottype,nbins)

if comparator:
    if dataPrint == "ON":
        print "\n{:s}".format(sys.argv[2])
    plotDistComparator(dataGrab(sys.argv[2]),plottype,nbins)

    for file in sys.argv[3:]:
        if dataPrint == "ON":
            print "\n{:s}".format(file)
        plotDist(dataGrab(file),plottype,nbins)
else:
    for file in sys.argv[2:]:
        if dataPrint == "ON":
            print "\n{:s}".format(file)
        plotDist(dataGrab(file),plottype,nbins)

ylabel("Frequency")

if plottype == "energy":
    xlabel("Energy (neV)")
    title("Energy Distribution")
    if comparator:
        legend(["Baseline Initial","Baseline Final","Single Interface","5 Layers","10 Layers","15 Layers","20 Layers"],loc=2)
    else:
        legend(["Baseline Initial","Baseline Final","5 Layers","10 Layers","15 Layers","20 Layers"],loc=2)


if plottype == "angular":
    xlabel("Energy Directed Downstream (neV)")
    title("Angular Distribution")
    if comparator:
        legend(["Baseline Initial","Baseline Final","Single Interface","5 Layers","10 Layers","15 Layers","20 Layers"],loc=1)
    else:
        legend(["Baseline Initial","Baseline Final","5 Layers","10 Layers","15 Layers","20 Layers"],loc=1)

show()
