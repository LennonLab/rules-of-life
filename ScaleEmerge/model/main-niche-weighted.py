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
OUT = open(GenPath + 'SimData-Random.csv', 'w+')
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
    #p = choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    p = 1.0
    choose1 = choice(['1','2'])
    choose2 = choice(['1','2'])
    choose3 = choice(['1','2'])
    cp1 = choice(['1','2'])
    cp2 = choice(['1','2'])
    cp3 = choice(['1','2'])
    ID = choose1+'-'+choose2+'-'+choose3+'-'+cp1+'-'+cp2+'-'+cp3+'-'+str(p)

    if ID not in cDict:
        r = lambda: randint(0,255)
        clr = '#%02X%02X%02X' % (r(),r(),r())
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

            if cp1 == '1': p1 = 1/len(RAC) * (e_simpson(RAC))
            else: p1 = 1/len(RAC) * (1-len(RAC)/sum(RAC))
            if np.random.binomial(1, p1) == 1:
                ''' differentiate '''
                RAC.sort()
                ranks = range(len(RAC))
                ps = np.array(RAC)/sum(RAC)
                if choose1 == '1': i1 = np.random.choice(ranks, 1, replace=True, p=ps)
                else: i1, i2 = np.random.choice(ranks, 2, replace=True)
                v = RAC.pop(i1)
                v1 = randint(1, v)
                v2 = v - v1
                RAC.extend([v1, v2])
                RAC = filter(lambda a: a != 0, RAC)

            if cp2 == '1': p1 = 1/len(RAC) * (e_simpson(RAC))
            else: p1 = 1/len(RAC) * (1-len(RAC)/sum(RAC))
            if np.random.binomial(1, p1) == 1:
                ''' rich to poor '''
                RAC.sort()
                ranks = range(len(RAC))
                ps = np.array(RAC)/sum(RAC)
                if choose2 == '1': i1, i2 = np.random.choice(ranks, 2, replace=True, p=ps)
                else: i1, i2 = np.random.choice(ranks, 2, replace=True)
                ii = [i1, i2]
                ii.sort()
                v = randint(1, RAC[ii[1]])
                RAC[ii[1]] -= v
                RAC[ii[0]] += v
                RAC = filter(lambda a: a != 0, RAC)

            if cp3 == '1': p1 = 1/len(RAC) * (e_simpson(RAC))
            else: p1 = 1/len(RAC) * (1-len(RAC)/sum(RAC))
            if np.random.binomial(1, p1) == 1:
                ''' poor to rich '''
                RAC.sort()
                ranks = range(len(RAC))
                ps = np.array(RAC)/sum(RAC)
                if choose3 == '1': i1, i2 = np.random.choice(ranks, 2, replace=True, p=ps)
                else: i1, i2 = np.random.choice(ranks, 2, replace=True)
                ii = [i1, i2]
                ii.sort(reverse=True)
                v = randint(1, RAC[ii[1]])
                RAC[ii[1]] -= v
                RAC[ii[0]] += v
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

        OUT = open(GenPath + 'SimData-Random.csv', 'a')
        outlist = [sim, N, S, ES, D, R, clr, choose1, choose2, choose3, p]
        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

        print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  int(N),
        print '  S:', '%4s' %  int(S), '  E:', '%4s' %  round(ES,2)
