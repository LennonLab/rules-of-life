from __future__ import division
from random import shuffle, choice, randint, randrange, seed
from os.path import expanduser
from numpy import log10
from scipy import stats
import numpy as np
import sys
import os


mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/Emergence/model")
GenPath = mydir + "GitHub/ScaleEmerge/results/simulated_data/"

col_headers = 'sim,total.abundance,species.richness,simpson.e,N.max,logmod.skew,clr,choose1,choose2,choose3,p'
OUT = open(GenPath + 'SimData-Random-simple.csv', 'w+')
print>>OUT, col_headers
OUT.close()


def GetRAD(vector):
    RAD = []
    unique = list(set(vector))
    for val in unique: RAD.append(vector.count(val))
    return RAD, unique


def e_simpson(sad):
    sad = filter(lambda a: a != 0, sad)
    D = 0.0
    N = sum(sad)
    S = len(sad)
    for x in sad: D += (x*x) / (N*N)
    E = round((1.0/D)/S, 4)
    return E


cDict = {}
for sim in range(10**6):
    p = choice([0.2, 0.4, 0.6, 0.8])
    choose1 = choice(['0'])
    choose2 = choice(['0'])
    choose3 = choice(['0'])
    cp1 = choice(['0'])
    cp2 = choice(['0'])
    cp3 = choice(['0'])
    ID = choose1+'-'+choose2+'-'+choose3+'-'+cp1+'-'+cp2+'-'+cp3+'-'+str(p)

    if ID not in cDict:
        if p == 0.2: clr = 'c'
        if p == 0.4: clr = 'DodgerBlue'
        if p == 0.6: clr = 'Blue'
        if p == 0.8: clr = 'DarkBlue'
        #r = lambda: randint(0,255)
        #clr = '#%02X%02X%02X' % (r(),r(),r())
        cDict[ID] = clr

    clr = cDict[ID]
    sample_size = 10
    ct = 0

    while ct < sample_size:
        ct += 1
        Ns = []
        Ss = []
        ESs = []
        Ds = []
        Rs = []

        N = int(10**np.random.uniform(0, 4))
        RAC = [N]
        for i in range(2000):

            p1 = p * (len(RAC) - RAC.count(1))/len(RAC)
            if len(RAC) < N and np.random.binomial(1, p1) == 1:
                RAC1 = filter(lambda a: a == 1, RAC)
                RAC2 = filter(lambda a: a > 1, RAC)

                RAC2.sort()
                v = RAC2.pop(0)
                v1 = randint(1, v)
                v2 = v - v1
                RAC2.extend([v1, v2])
                RAC = RAC1 + RAC2
                RAC = filter(lambda a: a != 0, RAC)

            p1 = p * (RAC.count(1))/len(RAC)
            if len(RAC) > 1 and np.random.binomial(1, p1) == 1:

                RAC.sort()
                v = randint(1, RAC[0])
                RAC[-1] += v
                RAC[0] -= v
                RAC = filter(lambda a: a != 0, RAC)


            if i > 1000 and i%100 == 0:
                N1 = sum(RAC)
                S = len(RAC)
                ES = e_simpson(RAC)
                Nm = max(RAC)

                skew = stats.skew(RAC)
                lms = log10(abs(float(skew)) + 1)
                if skew < 0: lms = lms * -1

                Ns.append(N1)
                Ss.append(S)
                ESs.append(ES)
                Ds.append(Nm)
                Rs.append(lms)

        N = np.mean(Ns)
        S = np.mean(Ss)
        ES = np.mean(ESs)
        D = np.mean(Ds)
        R = np.mean(Rs)

        OUT = open(GenPath + 'SimData-Random-simple.csv', 'a')
        outlist = [sim, N, S, ES, D, R, clr, choose1, choose2, choose3, p]
        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

        print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  int(N),
        print '  S:', '%4s' %  int(S), '  E:', '%4s' %  round(ES,2)
