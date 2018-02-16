from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


GenPath = mydir + "/results/simulated_data/"
#dat = open(GenPath + 'partitions-1000.csv', 'r')
dat = open(GenPath + 'partitions-1000-40.csv', 'r')


ranks, abundances = [], []
for line in dat:
    partition = eval(line)
    for rank, ab in enumerate(partition):
        ranks.append(rank+1)
        abundances.append(np.log10(ab))


fig = plt.figure()
xlab = 'Rank in abundance'
ylab = 'Abundance, '+r'$log_{10}$'

plt.style.use('dark_background')
plt.hexbin(ranks, abundances, mincnt=1, gridsize = 100, bins='log', cmap=plt.cm.jet)
plt.xlabel(xlab, fontsize=14)
plt.ylabel(ylab, fontsize=14)
plt.tick_params(axis='both', labelsize=8)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/Partitions-N=1000-S=40.png',
    dpi=400, bbox_inches = "tight")
plt.close()
