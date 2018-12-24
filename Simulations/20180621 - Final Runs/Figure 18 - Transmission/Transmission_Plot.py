#Plotting of correlation plots

from pylab import *
from matplotlib.ticker import FixedLocator

#DATA
layers = [0,1,2,3,5,10,15,20,30]
straight_07_layers = [56.288,42.438,34.716,29.844,26.404,16.602,12.010,9.348,6.15]
straight_07_layerserror = [0.1375,0.5109,0.5060,0.4227,0.3827,0.067,0.3097,0.1035,0.2909]
straight_10_layers = [58.076, 50.254, 45.943999999999996, 42.093999999999994, 39.596, 30.87, 25.548000000000002, 21.506, 15.92]
straight_10_layerserror = [0.5952142500000001,0.54284436,0.41740867000000004,0.6866076,0.21007142,0.24909836999999999,0.62247088,0.49333559,0.49295029999999995]
bent_10_layers = [49.198,41.55,38.99,35.525999999999996,33.488,26.741999999999997,22.158,19.186,14.59]
bent_10_layerserror = [0.71117508,0.33466400999999996,0.35637059,0.35788266,0.44007954,0.40763955,0.35322797,0.42009523000000004,0.22737634]

close()


#top plot
fig, axs = plt.subplots(nrows=1, ncols=1)

ax = axs
ax.errorbar(layers, array(straight_10_layers)/58.076, yerr=array(straight_10_layerserror)/58.076, fmt='b.-')
ax.errorbar(layers, array(straight_07_layers)/58.076, yerr=array(straight_07_layerserror)/58.076, fmt='r.--')

#ax.errorbar(layers, array(bent_10_layers)/58.076, yerr=array(bent_10_layerserror)/58.076, fmt='r.-')
ax.xaxis.set_ticks([0,5,10,15,20,25,30])

ax.set_ylabel("UCN Yield (arbitrary units)",fontsize=13)
ax.set_xlabel("Layers",fontsize=13)
ax.legend(["1.0 Specularity","0.7 Specularity"],fontsize=12)
xticks(fontsize=12)
yticks(fontsize=12)

with open("Transmission_Plot.txt","w") as file:
    file.write("Layers,Straight_0.7_Transmission_Percentage,Straight_0.7_Transmission_Error,Straight_1.0_Transmission_Percentage,Straight_1.0_Transmission_Error,Bent_1.0_Transmission_Percentage,Bent_1.0_Transmission_Error\n")
    for i in range(len(layers)):
        for list in [layers,straight_07_layers,straight_07_layerserror,straight_10_layers,straight_10_layerserror,bent_10_layers,bent_10_layerserror]:
            file.write(str(list[i]))
            file.write(",")
        file.write("\n")

savefig("Transmission_Plot.png",dpi=1000)
