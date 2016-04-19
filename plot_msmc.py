#!/usr/bin/python3.4

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plot_msmc_utils as pmu


plt.style.use('ggplot')

# convert msmc output to pop size and time in years
converted_input = pmu.popSizeStepPlot(sys.argv[1].rstrip("\n"),mu=1e-8,gen=5.0)

# convert y axis to 10^4
log10_y = [x / 10000 for x in converted_input[1]]

#plt the stepwise graph
plt.step(converted_input[0], log10_y)

#impose log scale on x axis
plt.xscale("log")


plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#label the axis
plt.xlabel("Time in Years (yr)")
plt.ylabel("Effective Population Size x10^4")

#give a tite
name = "_".join(sys.argv[1].rstrip("\n").split("_")[0:2])
plt.title(name)

#plot the figure
plt.savefig(name + ".png")
