from __future__ import division
from os.path import expanduser
import partitions as parts
import numpy as np
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/Emergence/model")
GenPath = mydir + "GitHub/ScaleEmerge/results/simulated_data/"

OUT = open(GenPath + 'partitions.csv', 'w+')
OUT.close()


for i in range(10000):
    N = int(10**np.random.uniform(1, 3))
    print N
    S = int(N)
    sample_size = 1
    D = {}
    name = 'best'

    zeros = False
    exact = False
    x = parts.rand_partitions(N, S, sample_size, name, D, zeros, exact)

    OUT = open(GenPath + 'partitions.csv', 'a')
    for i in x:
        print>>OUT, i
    OUT.close()
