# Spectrum Histogram Plotter
# plots histograms of energy spectra based upon lists of energies created by SpectrumAnalyzer.py
# To use: python SpectrumHistogramPlotter.py list1.txt list2.txt ...

from pylab import *
from sys import *

fig = figure()
ax = fig.add_subplot(111)

headers = []
filedata = []

def analyzeData(datafile):
    with open(datafile,'r') as file:
        data = [float(item) for item in file.read().splitlines()[1:]]
    text = datafile.replace("\ "," ").split('/')[-1][:-4].replace("_"," ")
    headers.append(text + "_Bin")
    headers.append(text + "_Frequency")
    frequency,binEdges=histogram(data,bins=25,density=True)
    frequency = array(frequency)*1000
    binCenters = 0.5*(binEdges[1:]+binEdges[:-1])
    filedata.append(binCenters)
    filedata.append(frequency)
    return binCenters,frequency

first_legend_handles=[]
second_legend_handles=[]

if len(argv) >= 2:
    baselineinitialx,baselineinitialy = analyzeData(argv[1])
    l1, = ax.plot(baselineinitialx,baselineinitialy,'k-',linewidth=2,label="Baseline Initial")
    first_legend_handles.append(l1)

if len(argv) >= 3:
    baselinefinalx, baselinefinaly = analyzeData(argv[2])
    l2, = ax.plot(baselinefinalx,baselinefinaly,'k--',linewidth=2,label="Baseline Final")
    first_legend_handles.append(l2)

if len(argv) >= 4:
    singleinterfacex, singleinterfacey = analyzeData(argv[3])
    l3, = ax.plot(singleinterfacex,singleinterfacey,'r:',linewidth=2,label="Single Interface")
    first_legend_handles.append(l3)

if len(argv) >= 5:
    fivelayersx, fivelayersy = analyzeData(argv[4])
    l4, = ax.plot(fivelayersx,fivelayersy,'b-',linewidth=1,label="5 Layers")
    second_legend_handles.append(l4)

if len(argv) >= 6:
    tenlayersx, tenlayersy = analyzeData(argv[5])
    l5, = ax.plot(tenlayersx,tenlayersy,'g-',linewidth=1,label="10 Layers")
    second_legend_handles.append(l5)

if len(argv) >= 7:
    fifteenlayersx, fifteenlayersy = analyzeData(argv[6])
    l6, = ax.plot(fifteenlayersx,fifteenlayersy,'c-',linewidth=1,label="15 Layers")
    second_legend_handles.append(l6)

if len(argv) >= 8:
    twentylayersx, twentylayersy = analyzeData(argv[7])
    l7, = ax.plot(twentylayersx,twentylayersy,'y-',linewidth=1,label="20 Layers")
    second_legend_handles.append(l7)

if len(argv) >= 9:
    thirtylayersx, thirtylayersy = analyzeData(argv[8])
    l8, = ax.plot(thirtylayersx,thirtylayersy,'-',color='darkorange',linewidth=1,label="30 Layers")
    second_legend_handles.append(l8)

ax.set_xlabel("Energy (neV)",fontsize=13)
ax.set_ylabel("Detected UCN Spectrum (arbitrary units)",fontsize=13)
xticks(fontsize=12)
yticks(fontsize=12)

first_legend = ax.legend(handles=first_legend_handles,loc=2,fontsize=12)
ax = gca().add_artist(first_legend)
legend(handles=second_legend_handles,loc=4,fontsize=12)

savefig("Spectral_Hardening_Plot.png",dpi=1000)

if True:
    with open("Spectrum_Hardening_Plot.txt","w") as file:
        file.write(",".join(headers)+"\n")
        for i in range(0,len(filedata[0])):
            for list in filedata[:-1]: file.write(str(list[i]) + ",")
            file.write(str(filedata[-1][i]))
            file.write("\n")

