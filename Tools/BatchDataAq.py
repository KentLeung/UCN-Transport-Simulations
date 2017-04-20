#Batch Data Aquisition Code - Erik Lutz
#Use this code to complete final data collection for simulations
#to run: python BatchDataAq.py binary directory1 directory2 directory3 ...

from pylab import *
from subprocess import *
from re import *
from sys import *
from time import *
import os

#PARAMETERS#
N = 3 #Number of runs
DET = [1,2,3]



fileoutput = []

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
        fileoutput.append(str(i)[1:-1])
    fileoutput.append("")
    fileoutput.append("")
    return output

def makePaths():
    paths = []
    for dir in argv[2:]:
        direc = []
        for i in dir.split('/'):
            if i not in argv[1].split('/'):
                direc.append(i)
        path = os.path.join(*direc) + '/'
        paths.append(path.replace(" ","\ "))
    return paths
tic = time()

binary = './' + argv[1].split('/')[-1]

print "################################################################"
print "\nRunning Batch Final Simulations.\n"

for path in makePaths():
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
