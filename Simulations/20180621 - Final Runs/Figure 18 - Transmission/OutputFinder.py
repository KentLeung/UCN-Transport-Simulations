from pylab import *
from sys import *

number = []
error = []

for item in argv[1:]:
    with open(item,'r') as file:
        raw = file.read().splitlines()
        number.append(float(raw[8])/100)
        error.append(float(raw[11])/100)

print number
print error

