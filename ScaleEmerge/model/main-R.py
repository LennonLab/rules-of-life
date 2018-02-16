from __future__ import division
from random import shuffle, choice, randint, sample
from os.path import expanduser
from numpy import log10
from scipy import stats
import numpy as np
import time
import copy
import sys

mydir = expanduser("~/")
sys.path.append(mydir + "GitHub/Emergence/model")
GenPath = mydir + "GitHub/ScaleEmerge/results/simulated_data/"

col_headers = 'sim,u,r,gr,mt,q,ct,total.abundance,species.richness,simpson.e,N.max,logmod.skew'
OUT = open(GenPath + 'SimData-mainR.csv', 'w+')
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


def output(iD, sD, rL, sim, ct, ps):
    r, u, gr, mt, q = ps
    IndIDs, SpIDs = [], []
    for k, v in iD.items():
            IndIDs.append(k)
            SpIDs.append(v['sp'])

    N = len(IndIDs)
    R = sum(rL)
    S = len(list(set(SpIDs)))

    if N > 0:
        RAD, splist = GetRAD(SpIDs)
        ES = e_simpson(RAD)
        Nm = max(RAD)

        skew = stats.skew(RAD)
        lms = log10(abs(float(skew)) + 1)
        if skew < 0: lms = lms * -1

        OUT = open(GenPath + 'SimData-MainR.csv', 'a')
        outlist = [sim, u, r, gr, mt, q, ct, N, S, ES, Nm, lms]
        outlist = str(outlist).strip('[]')
        outlist = outlist.replace(" ", "")
        print>>OUT, outlist
        OUT.close()

    print 'sim:', '%3s' % sim, 'ct:', '%3s' % ct,'  N:', '%4s' %  N, '  S:', '%4s' %  S,  '  R:', '%4s' %  R
    return


def immigration(sD, iD, rL, ps, sd=1):
    r, u, gr, mt, q = ps

    for j in range(sd):
        if sd == 1 and np.random.binomial(1, u) == 0: continue
        p = np.random.randint(1, 1000)
        if p not in sD:
            sD[p] = {'gr' : np.random.uniform(gr, 0)}
            sD[p]['mt'] = np.random.uniform(mt, 0)
            es = np.random.uniform(1, 100, 3)
            sD[p]['ef'] = es/sum(es)
            sD[p]['et'] = sample(range(len(rL)), 3)

        ID = time.time()
        iD[ID] = copy.copy(sD[p])
        iD[ID]['sp'] = p
        iD[ID]['q'] = np.random.uniform(0, q)

    return [sD, iD]


def ResIn(rL, ps):
    r, u, gr, mt, q = ps
    if np.random.binomial(1, u) == 1:
        rL += np.random.uniform(0, r, len(rL))
    return rL


def consume(iD, rL, ps):
    r, u, gr, mt, q = ps
    keys = list(iD)
    shuffle(keys)

    for k in keys:
        rl = iD[k]['et']
        c = choice(rl)
        i = rl.index(c)
        if rL[c] > 0:
            e = iD[k]['ef'][i] * iD[k]['q']
            iD[k]['q'] += min([rL[c], e])
            rL[c] -= min([rL[c], e])
    return [iD, rL]


def maintenance(iD):
    for k, v in iD.items():
        iD[k]['q'] -= v['mt']
        if v['q'] < v['mt']: del iD[k]
    return iD


def reproduce(sD, iD, ps, p = 0):
    for k, v in iD.items():
        if v['q'] > v['mt']*2 and np.random.binomial(1, v['gr']) == 1:
            iD[k]['q'] = v['q']/2.0
            i = time.time()
            iD[i] = copy.copy(iD[k])
    return [sD, iD]


def iter_procs(iD, sD, rL, ps, ct):
    procs = range(6)
    shuffle(procs)
    for p in procs:
        if p == 0: rL = ResIn(rL, ps)
        if p == 1: sD, iD = immigration(sD, iD, rL, ps)
        elif p == 2: iD, rL = consume(iD, rL, ps)
        elif p == 3: iD = maintenance(iD)
        elif p == 4: sD, iD = reproduce(sD, iD, ps)
    N = len(list(iD))
    return [iD, sD, rL, N, ct+1]


for sim in range(10**5):
    print '\n'
    gr = 10**np.random.uniform(-2, -1)
    mt = 10**np.random.uniform(-2, -1)
    q = 10**np.random.uniform(-2, -1)
    r = 10**np.random.uniform(-1, 1)
    u = 10**np.random.uniform(-2, 0)

    sD, iD, ct = {}, {}, 0
    ps = r, u, gr, mt, q
    rL = np.zeros(randint(3,100))

    sD, iD = immigration(sD, iD, rL, ps, 1000)
    while ct < 300:
        iD, sD, rL, N, ct = iter_procs(iD, sD, rL, ps, ct)
        if ct > 200 and ct%10 == 0: output(iD, sD, rL, sim, ct, ps)
