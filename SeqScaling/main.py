from __future__ import division
from random import choice, randint
import numpy as np
import sklearn.cluster
#from sklearn.cluster import Birch
#import sys



def GetRAD(vector):
    vector = vector.tolist()
    RAD, unique = [], list(set(vector))
    for val in unique: RAD.append(vector.count(val))
    return RAD, unique


def log_modulo_skew(sad):
    if len(sad) == 1: return float('NaN')
    x = np.array(sad)
    s = len(sad)
    skw = sum((x - np.mean(x))**3)/((s-1)*np.std(x)**3)
    lms = np.log10(abs(float(skw)) + 1)
    if skw < 0: lms = lms * -1
    return skw


def e_simpson(sad):
    if len(sad) == 0: return float('NaN')
    sad = filter(lambda a: a != 0, sad)
    D, N, S = [0, sum(sad), len(sad)]
    for x in sad: D += (x*x) / (N*N)
    E = round((1.0/D)/S, 4)
    return E



def iterative_levenshtein(s, t):
    """
        iterative_levenshtein(s, t) -> ldist
        ldist is the Levenshtein distance between the strings
        s and t.
        For all i and j, dist[i,j] will contain the Levenshtein
        distance between the first i characters of s and the
        first j characters of t

        This code was obtained from:
        https://www.python-course.eu/levenshtein_distance.php
        on 2018-02-08
    """
    rows = len(s)+1
    cols = len(t)+1
    dist = [[0 for x in range(cols)] for x in range(rows)]

    for i in range(1, rows):
        dist[i][0] = i

    for i in range(1, cols):
        dist[0][i] = i

    for col in range(1, cols):
        for row in range(1, rows):
            if s[row-1] == t[col-1]:
                cost = 0
            else:
                cost = 1
            dist[row][col] = min(dist[row-1][col] + 1,      # deletion
                                 dist[row][col-1] + 1,      # insertion
                                 dist[row-1][col-1] + cost) # substitution
    return dist[row][col]



def random_seqs(samplesize, length, m, generations):
    chars = ['A', 'T', 'G', 'C']
    o_seq = []
    for i in range(length): o_seq.append(choice(chars))

    seqs = [list(o_seq) for _ in xrange(samplesize)]
    for g in range(generations):

        #i1 = seqs[choice(range(samplesize))]
        #i2 = choice(range(samplesize))
        #seqs[i2] = i1

        for i, seq in enumerate(seqs):
            for j, p in enumerate(seq):
                if np.random.binomial(1, m) == 1:
                    seq[j] = choice(chars)

    for i, seq in enumerate(seqs):
        seqs[i] = "".join(str(x) for x in seq)

    return seqs



OUT = open('SimData.csv','w+')
h = 'rep,length,sample_size,mutation_rate,generations,N,S,E,Nmax,R'
print>>OUT, h
OUT.close()

for rep in range(10000):
    for samplesize in [5, 6, 7, 8, 9, 10, 12, 15, 17, 20, 25, 30, 35, 40, 60, 80, 120, 160, 240, 320, 480, 720, 960]:

        #samplesize = randint(4, 1000)
        samplesize = int(round(10**np.random.uniform(1, 3.5), 0))

        for length in [5]:
            for m in [0.00001]:
                for generations in [200]:

                    sequences = random_seqs(samplesize, length, m, generations)
                    sequences = np.asarray(sequences)
                    lev_similarity = -1*np.array([[iterative_levenshtein(w1,w2) for w1 in sequences] for w2 in sequences])

                    affprop = sklearn.cluster.AffinityPropagation(damping=0.5)
                    affprop.fit(lev_similarity)

                    vector = affprop.predict(lev_similarity)
                    sad, unique = GetRAD(vector)

                    OUT = open('SimData.csv','a')
                    N = sum(sad)
                    S = len(sad)
                    E = e_simpson(sad)
                    R = log_modulo_skew(sad)
                    Nmax = max(sad)

                    print rep, samplesize, length, m, generations, N, S

                    outlist = [rep, length, samplesize, m, generations, N, S, E, Nmax, R]
                    outlist = str(outlist).strip('[]')
                    outlist = outlist.replace(" ", "")
                    print>>OUT, outlist
                    OUT.close()
