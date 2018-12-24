#Bottling Effects - Erik Lutz
#Use this code to extract bottling data from events.sim files
#to run: python Bottling.py events.sim

from pylab import *
from sys import *

# OPTIONS
configuration = "Close Foil"
N = 10000
plot = False

# import data
file = open(argv[1],'r')
raw = [map(lambda x: int(x),line.split()[:4]) + map(lambda x: float(x),line.split()[4:]) for line in file.read().splitlines()[2:]]
file.close()

neutrons = [0 for i in range(0,N)]

top = []
bottom = []
front = []
back = []
left = []
right = []

if configuration == "Close Foil":
    bottlingregion = 5
    for i in range(0,len(raw)-1):
        if raw[i][2] == bottlingregion and raw[i+1][2] == bottlingregion: # find residency times
            neutrons[raw[i][0]] += raw[i+1][4] - raw[i][4]
        if raw[i][2] == bottlingregion and raw[i][1] in [1,2]:
            if raw[i][5] == -2.0: front.append(raw[i])
            if raw[i][5] == 2.0: back.append(raw[i])
            if raw[i][6] == -0.27: bottom.append(raw[i])
            if raw[i][6] == 2.23: top.append(raw[i])
            if raw[i][7] == 9.3001: left.append(raw[i])
            if raw[i][7] == 13.3001: right.append(raw[i])

if configuration == "Far Foil":
    bottlingregion = 18
    for i in range(0,len(raw)-1):
        if raw[i][2] == bottlingregion and raw[i+1][2] == bottlingregion: # find residency times
            neutrons[raw[i][0]] += raw[i+1][4] - raw[i][4]
        if raw[i][2] == bottlingregion and raw[i][1] in [1,2]:
            if raw[i][5] == -2.0: front.append(raw[i])
            if raw[i][5] == 2.0: back.append(raw[i])
            if raw[i][6] == -0.954533: bottom.append(raw[i])
            if raw[i][6] == 1.545467: top.append(raw[i])
            if raw[i][7] == 9.994513: left.append(raw[i])
            if raw[i][7] == 13.994513: right.append(raw[i])

if configuration == "Downward Guide":
    bottlingregion = 13
    for i in range(0,len(raw)-1):
        if raw[i][2] == bottlingregion and raw[i+1][2] == bottlingregion: # find residency times
            neutrons[raw[i][0]] += raw[i+1][4] - raw[i][4]
        if raw[i][2] == bottlingregion and raw[i][1] in [1,2]:
            if raw[i][5] == -2.0: front.append(raw[i])
            if raw[i][5] == 2.0: back.append(raw[i])
            if raw[i][6] == -2.712266: bottom.append(raw[i])
            if raw[i][6] == -0.212266: top.append(raw[i])
            if raw[i][7] == 9.300031: left.append(raw[i])
            if raw[i][7] == 13.300031: right.append(raw[i])

if configuration == "Sloped Guide":
    bottlingregion = 7
    for i in range(0,len(raw)-1):
        if raw[i][2] == bottlingregion and raw[i+1][2] == bottlingregion: # find residency times
            neutrons[raw[i][0]] += raw[i+1][4] - raw[i][4]
        if raw[i][2] == bottlingregion and raw[i][1] in [1,2]:
            if raw[i][5] == -2.0: front.append(raw[i])
            if raw[i][5] == 2.0: back.append(raw[i])
            if raw[i][6] == -0.77: bottom.append(raw[i])
            if raw[i][6] == 1.73: top.append(raw[i])
            if raw[i][7] == 8.652256: left.append(raw[i])
            if raw[i][7] == 12.652256: right.append(raw[i])

# find mean (deltaT)^2
bouncetime = []
for i in range(len(raw)-2): # look for bounce events
    
    if raw[i][2] == bottlingregion - 1 and raw[i+1][2] == bottlingregion: # neutron enters bottle
        bouncetime.append("Enter")
    
    if raw[i][2] == raw[i+1][2] == raw[i+1][2] == bottlingregion: # event in bottling region
        if raw[i][1] in [1,2]: # event is an intersection with region or cutplane
            if (raw[i+1][1] == 12) and (raw[i+2][1] == 3): # particle enters and exits a bounce
                bouncetime.append(raw[i][4])

    if raw[i][2] == bottlingregion and raw[i+1][2] == bottlingregion - 1: # neutron exits bottle
        bouncetime.append("Exit")

deltaTsq = []
deltaT = []

i = 0
while i < len(bouncetime) - 1: # talley real bounce times
    
    # exclude time between neutrons and time individual neutrons exit and enter bottle
    if bouncetime[i] == "Exit" or bouncetime[i] == "Enter":
        i += 1
        continue
    if bouncetime[i+1] == "Exit" or bouncetime[i+1] == "Enter":
        i += 1
        continue

    # now we know we have a single instance of the neutron being in the bottle
    deltaTsq.append((bouncetime[i+1]-bouncetime[i])**2)
    deltaT.append(bouncetime[i+1]-bouncetime[i])
    i += 1

bottledneutrons = [x for x in neutrons if x != 0]

print "\n###### Geometry Bottling Stats ######\n"

print "*   {:.1f}% of the neutrons entered the bottle.\n".format(float(len(bottledneutrons))*100/N)

if len(bottledneutrons):
    print "*   The average residency time for a neutron that enters the bottle is {:d} seconds.\n*   The average residency time overall is {:d} seconds.\n".format(int(mean(bottledneutrons)), int(mean(neutrons)))

    print "*   {:.1f}% of bottle interactions were with the top of the box.\n".format(float(len(top))/(len(top)+len(bottom)+len(back)+len(front)+len(left)+len(right))*100)

    print "*   Neutrons entering the bottle averaged {:d} bounces.\n".format(int(len(deltaTsq)/len(bottledneutrons)))

    print "*   deltaT^2 Stats:"
    print "*      mean(deltaT^2) =",mean(deltaTsq)
    print "*      median(deltaT^2) =",median(deltaTsq)
    print "*      min(deltaT^2) =",min(deltaTsq)
    print "*      max(deltaT^2) =",max(deltaTsq),"\n"

    print "###### Excel Output ######\n"

    print "{:.1f}".format(float(len(bottledneutrons))*100/N)
    print "{:d}".format(int(mean(bottledneutrons)))
    print "{:d}".format(int(mean(neutrons)))
    print "{:.1f}".format(float(len(top))/(len(top)+len(bottom)+len(back)+len(front)+len(left)+len(right))*100)
    print "{:d}".format(int(len(deltaTsq)/len(bottledneutrons)))
    print "{:.3f}".format(mean(deltaTsq))
    print "{:.3f}".format(median(deltaTsq))
    print "{:.3f}".format(min(deltaTsq))
    print "{:.3f}\n".format(max(deltaTsq))

if plot:
    close()


    fig = figure(2)

    ax = fig.add_subplot(121)
    ax.hist(bottledneutrons,25)
    ax.set_title("Residency Time Frequency Plot")
    ax.set_xlabel("Residency Time (s)")
    ax.set_ylabel("Frequecny")

    ax = fig.add_subplot(122)
    nbins = 30
    normalized = False

    def eExtract(data):
        return [sqrt(x[8]**2 + x[9]**2 + x[10]**2) for x in data]

    data = [eExtract(top),eExtract(bottom),eExtract(back),eExtract(front),eExtract(left),eExtract(right)]


    for dataset in data:
        n, binsinitial, patches = ax.hist(dataset,nbins,normed=normalized,alpha=0)
        bins = []
        for i in range(0,len(binsinitial)-1):
            bins.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])

        ax.plot(bins,n)

    ax.legend(["Top","Bottom","Back","Front","Left","Right"],loc=2)
    ax.set_title("Energy Spectra of Wall Interactions")
    ax.set_xlabel("Velocity (m/s)")
    ax.set_ylabel("Frequency")

    show()





