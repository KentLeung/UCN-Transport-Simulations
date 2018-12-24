# Spectrum Histogram Plotter
# plots histograms of energy spectra based upon lists of energies created by SpectrumAnalyzer.py
# To use: python SpectrumHistogramPlotter.py list1.txt list2.txt ...

from pylab import *
from sys import *

legend = []

fig = figure()
ax = fig.add_subplot(111)

for datafile in argv[1:]:
    with open(datafile,'r') as file:
        data = [float(item) for item in file.read().splitlines()[1:]]
    legend.append(datafile.replace("\ "," ").split('/')[-1][:-4].replace("_"," "))
    energy,binEdges=histogram(data,bins=25,density=True)
    binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
    ax.plot(binCenters,energy)

ax.legend(legend)
ax.set_xlabel("Energy (neV)")
ax.set_ylabel("Normalized Frequency")
show()
