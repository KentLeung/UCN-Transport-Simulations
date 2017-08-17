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
Regionfile = "Regionfile"


#files
#Regionfile = argv[3]

#parameters to vary over
label1 = 'Spec'
param1 = ["0.0","0.5","0.97"]
label2 = 'Absorb'
param2 = ["0.287"]

#Other parameters
regs = [8]
DET = [1,2,3,4]
totalout = [[] for i in param1]
count = 1
numofruns = 10

tic = time()
storage = []
stddeviation = []
stedevTemp = []

for i in range(0,len(param1)):
    for j in range(0,len(param2)): #Loop over parameter space
    
        #print("Run {:d} of {:d}.".format(count,len(param1)*len(param2)))
        print("Run {:d}.{:d}.".format(i+1,j+1))
    
        output = zeros(len(DET)+1)

        #change regionfile
        editRegion(Regionfile,regs,label1,param1[i])
        editRegion(Regionfile,regs,label2,param2[j])

        #run simulation
        for run in range(0,numofruns):
            tempoutput = simRun(binary,directory,DET)
            storage.append(tempoutput)
            for det in range(0,len(output)):
                output[det] += tempoutput[det]
        for average in range(0,len(output)):
            output[average] = output[average]/float(numofruns)
        for stddev in range(0,len(output)):
            stedevTemp.append(numpy.std([row[stddev] for row in storage],ddof=1))
        stddeviation.append(stedevTemp)


        totalout[i].append(output)
        toc = time()
        elapsedtime = (toc - tic)/3600.
        estimatedtime = (elapsedtime/count)*(len(param1)*len(param2)-count)
        print("Elapsed Time: {:.2f} hours. Estimated Time Remaining: {:.2f} hours.".format(elapsedtime,estimatedtime))
        count += 1
        storage = []
        stedevTemp = []

print("\n\n")

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
editRegion(Regionfile,regs,label1,"0.97")
editRegion(Regionfile,regs,label2,"0.287")

print("\a\a\a")

