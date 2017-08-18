#Bottling Effects - Erik Lutz
#Use this code to extract bottling data from events.sim files
#to run: python Bottling.py events.sim

from pylab import *
from sys import *

# OPTIONS
bottlingregion = 5
N = 10000

# import data
file = open(argv[1],'r')
raw = [map(lambda x: int(x),line.split()[:4]) + map(lambda x: float(x),line.split()[4:]) for line in file.read().splitlines()[2:]]
file.close()

#for line in raw[:5]: print line

neutrons = [0 for i in range(0,N)]

top = []
bottom = []
front = []
back = []
left = []
right = []

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


bottledneutrons = [x for x in neutrons if x != 0]

print "The average residency time for a neutron that enters the bottle is {:d} seconds. The average residency time overall is {:d} seconds.".format(int(mean(bottledneutrons)), int(mean(neutrons)))

print "{:.1f}% of bottle interactions were with the top of the box.".format(float(len(top))/(len(top)+len(bottom)+len(back)+len(front)+len(left)+len(right))*100)

close()

fig, axs = plt.subplots(nrows=1,ncols=2)

ax = axs[0]
ax.hist(bottledneutrons,25)
ax.set_title("Residency Time Frequency Plot")
ax.set_xlabel("Residency Time (s)")
ax.set_ylabel("Frequecny")

ax = axs[1]
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





