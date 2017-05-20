#Library for common functions used in tool scripts

from pylab import *
from subprocess import *
from re import *
from sys import *
from time import *
import os

##Regionfile Editor Functions##
def formatRegion(line): #creates a regionfile line with parameters listed in argument
    a = ""
    a = a + "{:6s}".format(line[0])
    a = a + "{:7s}".format(line[1])
    a = a + "{:20s}".format(line[2])
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

def editRegion(file, regs, param, value): #change a value in a regionfile
    #get original file
    Regfile = open(file, "r")
    raw = Regfile.read().splitlines()
    Regfile.close()

    #vars
    header = "Reg#  RType  BP(x,y,z)           Dim[m]          Orient          Grad B [T/m]  Spec   Loss    Depol   WPot[neV]  BPot[neV]  Scat   Absorb[1/s]   Det  PM  SM  LM  DM\n"
    raw = [i.split() for i in raw[1:-1]]
    options = {'Reg':0, 'RType':1, 'BP':2, 'Dim':3, 'Orient':4, 'Grad B':5, 'Spec':6, 'Loss':7, 'Depol':8, 'WPot':9, 'BPot':10, 'Scat':11, 'Absorb':12, 'Det':13, 'PM':14, 'SM':15, 'LM':16, 'DM':17}

    #change line/lines
    for i in raw:
        if int(i[0]) in regs:
            i[options[param]] = value

    #construct new file
    output = [header]
    for i in raw:
        output.append(formatRegion(i))
    output.append("/")

    #write out
    Regfile = open(file, "w")
    Regfile.writelines(output)
    Regfile.close()

##Connexfile Editor Functions##
def blankConnex(file,num): #creates file for contiguous connexfile with "num+1" regions
    #create lines
    header1 = "connex: The first region should be region 0, and a special-handling code must be used to specify how to treat its cut-plane.\n"
    header2 = "-----------------------------------------------------------------------------------------------------------------------------\n"
    header3 = "Region#     Connects through its own cut-plane to:     AND also to:     AND also to:     AND also to:     AND also to:     AND also to:\n"
    line1 = "0           0                                             1/\n"
    output = [header1,header2,header3,line1]

    for i in range(1,min(10,num)):
        output.append("{:d}{:12d}{:46d}/\n".format(i,i-1,i+1))

    if num >= 10:
        for i in range(10,num):
            output.append("{:d}{:11d}{:46d}/\n".format(i,i-1,i+1))
        output.append("{:d}{:11d}/\n".format(num,num-1))
    else:
        output.append("{:d}{:12d}/\n".format(num,num-1))

    #write file
    Confile = open(file, "w")
    Confile.writelines(output)
    Confile.close()

def editHandling(file,regs,code): #adds or modifies special handling codes in connexfiles
    #get lines
    Confile = open(file, "r")
    raw = [i.split() for i in Confile.read().splitlines()[3:]]
    Confile.close()
    
    #get existing codes
    sCodes = ["" for i in zeros(len(raw))]
    for i in range(0,len(raw)):
        if len(raw[i][0].split(',')) == 2:
            sCodes[i] = raw[i][0].split(',')[1]

    #edit codes
    for i in range(0,len(raw)):
        if int(raw[i][0].split(',')[0]) in regs:
            sCodes[i] = str(code)

    lines = ["connex: The first region should be region 0, and a special-handling code must be used to specify how to treat its cut-plane.\n","-----------------------------------------------------------------------------------------------------------------------------\n","Region#     Connects through its own cut-plane to:     AND also to:     AND also to:     AND also to:     AND also to:     AND also to:\n"]

    for i in range(0,len(raw)):
        if sCodes[i] != "":
            raw[i][0] = raw[i][0].split(',')[0] + "," + sCodes[i]

    for i in raw[:-1]:
        lines.append("{:12s}{:46s}{:s}\n".format(i[0],i[1],i[2]))

    lines.append("{:12s}{:s}".format(raw[-1][0],raw[-1][1]))

    #write out
    Confile = open(file, "w")
    Confile.writelines(lines)
    Confile.close()

##Spectrum Analyzer Functions##
def spectralData(directory):
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

def plotSpectralBaselineDist((poofs,cphits,detid),plottype,nbins,normalized,dataPrint = "OFF"):
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
        
        plot(bins1,n1*100,linewidth=3,color='k',ls='-')
        plot(bins2,n2*100,linewidth=3,color='k',ls='--')

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
        
        plot(bins1,n1*100,linewidth=3,color='k',ls='-')
        plot(bins2,n2*100,linewidth=3,color='k',ls='--')

    if dataPrint == "ON":
        print "\nBaseline Initial"
        for i in range(0,len(bins1)):
            print bins1[i],n1[i]
        print "\nBaseline Final"
        for i in range(0,len(bins2)):
            print bins2[i],n2[i]


def plotSpectralDist((poofs,cphits,detid),plottype,nbins,normalized,dataPrint = "OFF"):
    if plottype == "energy":
        n2, binsinitial, patches = hist([i[2] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized,alpha=0)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2*100,'-')
    
    if plottype == "angular":
        n2, binsinitial, patches = hist([i[3] for i in cphits[detid[detid.keys()[-1]]]],nbins,normed=normalized)
        bins2 = []
        for i in range(0,len(binsinitial)-1):
            bins2.append((1./2.)*(binsinitial[i+1]-binsinitial[i])+binsinitial[i])
        
        plot(bins2,n2*100,'-')

    if dataPrint == "ON":
        for i in range(0,len(bins2)):
            print bins2[i],n2[i]

def plotLayersSpectral(baseline,eventsfiles,plottype,nbins,normalized,dataPrint = "OFF"):
    close()
    plotSpectralBaselineDist(spectralData(baseline),plottype,nbins,normalized)
    for file in eventsfiles:
        if dataPrint == "ON":
            print "\n{:s}".format(file)
        plotSpectralDist(spectralData(file),plottype,nbins,normalized)

    ylabel("Frequency (%)")

    if plottype == "energy":
        xlabel("Energy (neV)")
        title("Energy Distribution")
        legend(["Baseline Initial","Baseline Final","5 Layers","10 Layers","15 Layers","20 Layers"],loc=2)


    if plottype == "angular":
        xlabel("Energy Directed Downstream (neV)")
        title("Angular Distribution")
        legend(["Baseline Initial","Baseline Final","5 Layers","10 Layers","15 Layers","20 Layers"],loc=1)

    show()

##Batch Data Aquisition Functions##
def simRun(binary,directory,detectors=[1,2,3,4,5,6,7,8,9,10]):
    sim = Popen(binary + " " + directory, shell=True, stdout=PIPE)
    output, error = sim.communicate()
    uloss = [int(i) for i in findall("\d", findall("\(\d*\) of the particles",output)[0])]
    dets = [int(i[1]) for i in [findall("\d+",i) for i in findall("Detector \d* -> \d*",output)]]
    return [dets[i-1] for i in detectors] + uloss

def sortFiles(run,directory):
    call("mkdir -p " + directory + "Run\ {:d}".format(run),shell=True)
    call("mv " + directory + "events.sim " + directory + "geomout.sim " + directory + "detectors.sim " + directory + "Run\ {:d}".format(run),shell=True)

def finalRuns(binary,directory,detectors,runs=5):
    output = []
    for i in range(1,runs+1):
        output.append(simRun(binary,directory,detectors))
        sortFiles(i,directory)
    fileoutput.append(directory)
    for i in output:
        fileoutput.append(str(i)[1:-1].replace(",",""))
    fileoutput.append("")
    fileoutput.append("")
    return output

def makePathsFromBinary(binary,dirs):
    paths = []
    for dir in dirs:
        direc = []
        for i in dir.split('/'):
            if i not in binary.split('/'):
                direc.append(i)
        path = os.path.join(*direc) + '/'
        paths.append(path.replace(" ","\ "))
    return paths

def makePaths(dirs):
    paths = []
    for dir in dirs:
        direc = []
        for i in dir.split('/'):
            direc.append(i)
        path = "/" + os.path.join(*direc)
        paths.append(path.replace(" ","\ "))
    return paths

def batchDataAq(pathtobinary,dirs,DET,N):
    tic = time()

    binary = './' + pathtobinary.split('/')[-1]

    print "################################################################"
    print "\nRunning Batch Final Simulations.\n"

    for path in makePaths(dirs):
        print "Simulating runs in {:s}.\nData:\n".format(path)
        data =  finalRuns(binary,path,DET,N)
        for i in data:
            print i
        print "\n\n"

    file = open("Batch Output.txt","w")
    file.write('\n'.join(fileoutput))
    file.close()

    toc = time()
    elapsedtime = (toc - tic)/60.

    print "\a\a\a\a"

    print "All simulations finished in {:.2f} minutes.\n".format(elapsedtime)
    print "################################################################"


##More Functions##
