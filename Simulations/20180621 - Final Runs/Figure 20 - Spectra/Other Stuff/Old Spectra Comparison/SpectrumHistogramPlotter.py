# Spectrum Histogram Plotter
# plots histograms of energy spectra based upon lists of energies created by SpectrumAnalyzer.py
# To use: python SpectrumHistogramPlotter.py list1.txt list2.txt ...

from pylab import *
from sys import *

legend = []

fig = figure()
ax = fig.add_subplot(111)

headers = []
filedata = []

for datafile in argv[1:]:
    with open(datafile,'r') as file:
        data = [float(item) for item in file.read().splitlines()[1:]]
    text = datafile.replace("\ "," ").split('/')[-1][:-4].replace("_"," ")
    legend.append(text)
    headers.append(text + "_Bin")
    headers.append(text + "_Frequency")
    frequency,binEdges=histogram(data,bins=25,density=True)
    binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
    ax.plot(binCenters,frequency)
    filedata.append(binCenters)
    filedata.append(frequency)

ax.legend(legend)
ax.set_xlabel("Energy (neV)")
ax.set_ylabel("Normalized Frequency")
show()

if True:
    with open("SpectrumHistogramData.txt","w") as file:
        file.write(",".join(headers)+"\n")
        for i in range(0,len(filedata[0])):
            for list in filedata[:-1]: file.write(str(list[i]) + ",")
            file.write(str(filedata[-1][i]))
            file.write("\n")

