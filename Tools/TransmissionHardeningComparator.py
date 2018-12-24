# Plotter for spectral hardening vs transmission, uses energy lists from SpectralAnalyzer.py
# to use:
#    python TransmissionHardeningComparator.py list1.txt list2.txt ...


from pylab import *
from sys import *

eCutoff = 190

lists = []
totals = []
for filepath in argv[1:]:
    with open(filepath,'r') as file:
        raw = file.read().splitlines()
        totals.append(float(raw[0][23:]))
        lists.append([float(entry) for entry in raw[1:]])

hardening = []
losses = []

for i in range(len(totals)):
    losses.append((totals[i]-len(lists[i]))/totals[i])
    hardening.append(len(lists[i])/float(sum(i>eCutoff for i in lists[i])))

plot(losses,hardening,'o')
xlabel("Loss Fraction")
ylabel("Number of Neutrons/Number of Neutrons over {:d} neV".format(eCutoff))
show()

