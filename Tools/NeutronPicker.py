from pylab import *
from sys import *

def readline(line):
    return [float(item) for item in line.split()]

def neutronNum(line):
    return int(readline(line)[0])

def neutronE(line):
    data = readline(line)
    return (1/2.)*1.67493e-27*(data[8]**2 + data[9]**2 + data[10]**2)/1.602177e-28

target = int(argv[2])

t = []
y = []

with open(argv[1],'r') as file:
    for i,line in enumerate(file):
        if i<=1: continue
        if neutronNum(line) == target:
            print [int(item) for item in readline(line)[:4]],neutronE(line)
            t.append(float(readline(line)[4]))
            y.append(float(readline(line)[6]))

plot(t,y)
show()
