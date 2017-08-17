#This file runs a spread of simulations over a variety of SD2 conditions
#python SD2ConditionSpread.py binary directory regionfile
from __future__ import print_function

#from pylab import *
from sys import *
from UCNToolsLib import *
from time import *
import numpy

#Temporary Variables
binary = "./a.out"
directory = "./"
Regionfile = "Connexfile"


#files
#Regionfile = argv[3]

#parameters to vary over
param1 = ["7(0.005)","7(0.010)","7(0.025)","7(0.050)","7(0.060)"]

#Other parameters
regs = [14]
DET = [1,2,3,4]
totalout = [[] for i in param1]
count = 1
numofruns = 10


tic = time()
storage = []
stddeviation = []
stedevTemp = []

for i in range(0,len(param1)):
    #print("Run {:d} of {:d}.".format(count,len(param1)*len(param2)))
    print("Run {:d}".format(i+1))

    output = zeros(len(DET)+1)

    #change regionfile
    editHandling(Regionfile,regs,param1[i])

    #run simulation
    for ij in range(0,numofruns):
        tempoutput = simRun(binary,directory,DET)
        storage.append(tempoutput)
        for ijk in range(0,len(output)):
            output[ijk] += tempoutput[ijk]
    for ijkl in range(0,len(output)):
        output[ijkl] = output[ijkl]/float(numofruns)
    for stddev in range(0,len(output)):
        stedevTemp.append(numpy.std([row[stddev] for row in storage],ddof=1))
    stddeviation.append(stedevTemp)

    totalout[i].append(output)
    toc = time()
    elapsedtime = (toc - tic)/3600.
    estimatedtime = (elapsedtime/count)*(len(param1)-count)
    print("Elapsed Time: {:.2f} hours. Estimated Time Remaining: {:.2f} hours.".format(elapsedtime,estimatedtime))
    count += 1
    storage = []
    stedevTemp = []


file = open("spreadoutput.txt","w")

for condition in totalout:
    string = ""
    for det in condition[0]:
        string += "{:f},".format(det)
    string = string[:-1]
    print(string)
    file.write(string)
    file.write("\n")

print("\n")
file.write("\n")

for condition in stddeviation:
    string = ""
    for det in condition:
        string += "{:f},".format(det)
    string = string[:-1]
    print(string)
    file.write(string)
    file.write("\n")


file.close()


#Change regionfile back
editHandling(Regionfile,regs,"7(0.025)")

print("\a\a\a")

