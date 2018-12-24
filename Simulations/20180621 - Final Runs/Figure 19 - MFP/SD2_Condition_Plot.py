#Plotting of SD2 Condition plot

from pylab import *
from sys import *

mfp = [1,2,3,4,8]

max = 7772.6

with open(argv[1],'r') as file:
    raw = [i.split(",") for i in file.read().splitlines()[8:13]]
    zerolayershundredms = [float(i[0])/max for i in raw]
    zerolayerstwentyms = [float(i[1])/max for i in raw]

with open(argv[2],'r') as file:
    raw = [i.split(",") for i in file.read().splitlines()[8:13]]
    fourlayershundredms = [float(i[0])/max for i in raw]
    fourlayerstwentyms = [float(i[1])/max for i in raw]

with open(argv[3],'r') as file:
    raw = [i.split(",") for i in file.read().splitlines()[8:13]]
    eightlayershundredms = [float(i[0])/max for i in raw]
    eightlayerstwentyms = [float(i[1])/max for i in raw]

close()

marksize = 4

l1, = plot(mfp,zerolayershundredms,'ko-',ms=marksize,label="100 ms, 0 layers")
l2, = plot(mfp,fourlayershundredms,'bs-',ms=marksize,label="100 ms, 4 layers")
l3, = plot(mfp,eightlayershundredms,'gx-',ms=marksize,label="100 ms, 8 layers")

l4, = plot(mfp,zerolayerstwentyms,'ko--',ms=marksize,label="20 ms, 0 layers")
l5, = plot(mfp,fourlayerstwentyms,'bs--',ms=marksize,label="20 ms, 4 layers")
l6, = plot(mfp,eightlayerstwentyms,'gx--',ms=marksize,label="20 ms, 8 layers")

ylim([0,1.2])

xlabel("Mean Free Path (cm)",fontsize=13)
xticks([1,2,3,4,8],fontsize=12)
yticks(fontsize=12)

ylabel("UCN Yield (relative units)",fontsize=13)

first_legend = legend(handles=[l1,l2,l3],loc=2,fontsize=10)
ax = gca().add_artist(first_legend)
legend(handles=[l4,l5,l6],loc=4,fontsize=10)

with open("SD2_Condition_Plot.txt","w") as file:
    file.write("mfp,zero_hundred,four_hundred,eight_hundred,zero_twenty,four_twenty,eight_twenty\n")
    for i in range(len(mfp)):
        for list in [mfp,zerolayershundredms,fourlayershundredms,eightlayershundredms,zerolayerstwentyms,fourlayerstwentyms,eightlayerstwentyms]:
            file.write(str(list[i]))
            file.write(",")
        file.write("\n")

savefig("SD2_Condition_Plot.png",dpi=1000)
