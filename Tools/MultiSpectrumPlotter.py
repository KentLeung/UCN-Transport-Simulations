#Plotter for particle energy spectra

from pylab import *
import sys

file = open(sys.argv[1],"r")
raw = [map(float,i.split(',')) for i in [i[1:-1] for i in file.read().splitlines()]]
file.close()

#For Spectral Plot
fiveibin       = raw[0]
fivein         = raw[1]
fivefbin       = raw[2]
fivefn         = raw[3]
tenibin        = raw[4]
tenin          = raw[5]
tenfbin        = raw[6]
tenfn          = raw[7]
fifteenibin    = raw[8]
fifteenin      = raw[9]
fifteenfbin    = raw[10]
fifteenfn      = raw[11]
twentyibin     = raw[12]
twentyin       = raw[13]
twentyfbin     = raw[14]
twentyfn       = raw[15]
initialibin     = raw[16]
initialin       = raw[17]
initialfbin     = raw[18]
initialfn       = raw[19]

plot(fiveibin,[100*i for i in fivein],'k-',linewidth=3)
plot(initialfbin,[100*i for i in initialfn],'k--',linewidth=3)
plot(fivefbin,[100*i for i in fivefn],'r-',markersize=4)
plot(tenfbin,[100*i for i in tenfn],'b-')
plot(fifteenfbin,[100*i for i in fifteenfn],'g-',markersize=4)
plot(twentyfbin,[100*i for i in twentyfn],'m-',markersize=3)
legend(["Initial Distribution","Baseline Final Distribution","5 layers - at detector","10 layers - at detector","15 layers - at detector","20 layers - at detector"])
xlabel("Energy Directed Downstream (neV)")
#xlabel("Energy (neV)")
#ylabel("Frequency (%)")
yticks([])
xlim((0,500))
ylim((0,2))
show()
