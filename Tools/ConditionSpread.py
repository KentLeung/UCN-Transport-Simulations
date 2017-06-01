#This file runs a spread of simulations over a variety of SD2 conditions
#python SD2ConditionSpread.py binary directory regionfile
from __future__ import print_function

#from pylab import *
from sys import *
from UCNToolsLib import *
from time import *

#Temporary Variables
binary = "./a.out"
directory = "./"
Regionfile = "Regionfile"


#files
#Regionfile = argv[3]

#parameters to vary over
label1 = 'Scat'
param1 = ["0.01","0.02","0.03","0.04"]
label2 = 'Absorb'
param2 = ["50.00","25.00","16.66","12.50"]

#Other parameters
regs = [1]
DET = [1,2,3]
totalout = [[] for i in param1]
count = 1

tic = time()

for i in range(0,len(param1)):
    for j in range(0,len(param2)): #Loop over parameter space
    
        #print("Run {:d} of {:d}.".format(count,len(param1)*len(param2)))
        print("Run {:d}.{:d}.".format(i+1,j+1))
    
        output = zeros(len(DET)+1)

        #change regionfile
        editRegion(Regionfile,regs,label1,param1[i])
        editRegion(Regionfile,regs,label2,param2[j])

        #run simulation
        for ij in range(0,5):
            tempoutput = simRun(binary,directory,DET)
            for ijk in range(0,len(output)):
                output[ijk] += tempoutput[ijk]
        for ijkl in range(0,len(output)):
            output[ijkl] = output[ijkl]/5.

        totalout[i].append(output)
        toc = time()
        elapsedtime = (toc - tic)/3600.
        estimatedtime = (elapsedtime/count)*(len(param1)*len(param2)-count)
        print("Elapsed Time: {:.2f} hours. Estimated Time Remaining: {:.2f} hours.".format(elapsedtime,estimatedtime))
        count += 1

file = open("spreadoutput.txt","w")
for d in range(0,len(DET)+1):
    print("\n\n")
    file.write("\n\n")
    for i in range(0,len(param1)):
        print("")
        file.write("\n")
        for j in range(0,len(param2)): #construct output
            print(totalout[i][j][d],end=",")
            file.write(str(totalout[i][j][d]) + ",")
print("\n")


#Change regionfile back
editRegion(Regionfile,regs,label1,"0.04")
editRegion(Regionfile,regs,label2,"25.0")

print("\a\a\a")

