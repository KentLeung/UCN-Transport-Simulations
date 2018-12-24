#This file runs a spread of simulations over a variety of SD2 conditions
#python ConditionSpread.py binary directory regionfile
from __future__ import print_function

#from pylab import *
from sys import *
from UCNToolsLib import *
from time import *
import numpy

#Temporary Variables
binary = "./transmission.out"
directory = "./"
Regionfile = "Regionfile"


#files
#Regionfile = argv[3]

#parameters to vary over
label1 = 'Spec'
param1 = ["0.700"]
param1endvalue = "1.000"

label2 = 'WPot'
param2 = ["350.0"]
param2endvalue = "350.0"

concurrentvariation = False
if concurrentvariation: #concurrent variation,
    label1b = 'BP'
    param1b = [">0.0,0.01,1.47"] # must be same length as param1
    param1bendvalue = ">0.0,0.01,1.47"
    regs1b = [6]

#Other parameters
regs = [0]
DET = [1,2,3]
totalout = [[[] for i in param1] for i in range(len(DET)+1)]
count = 1
numofruns = 5

tic = time()
storage = []
stddeviation = [[[] for i in param1] for i in range(len(DET)+1)]
stedevTemp = []

for i in range(0,len(param1)):
    for j in range(0,len(param2)): #Loop over parameter space
    
        #print("Run {:d} of {:d}.".format(count,len(param1)*len(param2)))
        print("Run {:d}.{:d}.".format(i+1,j+1))
    
        output = zeros(len(DET)+1)

        #change regionfile
        editRegion(Regionfile,regs,label1,param1[i])
        editRegion(Regionfile,regs,label2,param2[j])
        if concurrentvariation: editRegion(Regionfile,regs1b,label1b,param1b[i])

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
    

        for k in range(len(DET)+1):
            totalout[k][i].append(output[k])
            stddeviation[k][i].append(stedevTemp[k])
        
        
        toc = time()
        elapsedtime = (toc - tic)/3600.
        estimatedtime = (elapsedtime/count)*(len(param1)*len(param2)-count)
        print("Elapsed Time: {:.2f} hours. Estimated Time Remaining: {:.2f} hours.".format(elapsedtime,estimatedtime))
        count += 1
        storage = []
        stedevTemp = []

print("\n\n")

file = open("spreadoutput.txt","w")

print("OUTPUT: (each row is from param1, each column is from param2)")
file.write("OUTPUT: (each row is from param1, each column is from param2)\n")

print("\nparam1: Varied \"{:s}\" over".format(label1),param1)
print("param2: Varied \"{:s}\" over".format(label2),param2)
print("Varied regions:", regs)
print("Ran each simulation {:d} times.".format(numofruns))

file.write("\nparam1: Varied \"{:s}\" over ".format(label1))
file.write(str(param1))
file.write("\nparam2: Varied \"{:s}\" over".format(label2))
file.write(str(param2))
file.write("\nVaried regions:")
file.write(str(regs))
file.write("\nRan each simulation {:d} times.\n".format(numofruns))

for k in range(len(DET)): # print out detector data
    print("\nDetector {:d} Means".format(DET[k]))
    file.write("\nDetector {:d} Means\n".format(DET[k]))
    for line in totalout[k]:
        string = ""
        for entry in line:
            string += "{:f},".format(entry)
        string = string[:-1]
        print(string)
        file.write(string)
        file.write("\n")
    print("\nDetector {:d} Standard Deviations".format(DET[k]))
    file.write("\nDetector {:d} Standard Deviations\n".format(DET[k]))
    for line in stddeviation[k]:
        string = ""
        for entry in line:
            string += "{:f},".format(entry)
        string = string[:-1]
        print(string)
        file.write(string)
        file.write("\n")

print("\nUnphysical Loss Means")
file.write("\nUnphysical Loss Means\n")
for line in totalout[-1]:
    string = ""
    for entry in line:
        string += "{:f},".format(entry)
    string = string[:-1]
    print(string)
    file.write(string)
    file.write("\n")
print("\nUnphysical Loss Standard Deviations")
file.write("\nUnphysical Loss Standard Deviations\n")
for line in stddeviation[-1]:
    string = ""
    for entry in line:
        string += "{:f},".format(entry)
    string = string[:-1]
    print(string)
    file.write(string)
    file.write("\n")

file.close()


#Change regionfile back
editRegion(Regionfile,regs,label1,param1endvalue)
editRegion(Regionfile,regs,label2,param2endvalue)
if concurrentvariation: editRegion(Regionfile,regs1b,label1b,param1bendvalue)

print("\a\a\a")

