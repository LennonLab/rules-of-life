from __future__ import division
from os.path import expanduser
import partitions as parts
import numpy as np
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/Emergence/model")
GenPath = mydir + "GitHub/ScaleEmerge/results/simulated_data/"

N = 1000
S = 40 #int(N)
sample_size = 1
D = {}

name = 'best'
zeros = False
exact = False

#OUT = open(GenPath + 'partitions-'+str(N)+'.csv', 'w+')
#OUT.close()

for i in range(10000):
    x = parts.rand_partitions(N, S, sample_size, name, D, zeros, exact)
    OUT = open(GenPath + 'partitions-'+str(N)+'-'+str(S)+'.csv', 'a')
    print>>OUT, x[0]
    print i
    OUT.close()
