from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
import scipy as sc
from scipy import stats
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table
from scipy.stats.kde import gaussian_kde


mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")


def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list), max(_list), n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


def figplot(df3, Ys, xlabs, ylab, fig):

    clrs = list(df3['color'])
    unique_colors = list(set(clrs))
    Rmlist = [list([]) for _ in xrange(len(unique_colors))]
    Dmlist = [list([]) for _ in xrange(len(unique_colors))]
    Emlist = [list([]) for _ in xrange(len(unique_colors))]
    Smlist = [list([]) for _ in xrange(len(unique_colors))]

    for i, clr in enumerate(unique_colors):
        print i
        df4 = df3[df3['color'] == clr]
        sims = list(df4['sim'])
        unique_sims = list(set(sims))

        for sim in unique_sims:
            for Y in Ys:
                df5 = df4[df4['sim'] == sim]
                y2 = list(df5[Y])
                x2 = list(df5['N'])

                m, b, r, p, std_err = stats.linregress(x2, y2)
                if Y == 'R': Rmlist[i].append(m)
                elif Y == 'D': Dmlist[i].append(m)
                elif Y == 'E': Emlist[i].append(m)
                elif Y == 'S': Smlist[i].append(m)

    clist = list(unique_colors)
    alist = [0] * len(unique_colors)
    llist = [0] * len(unique_colors)

    for i, clr in enumerate(clist):
        kernel = 0.5
        D = get_kdens_choose_kernel(Rmlist[i], kernel)
        x = D[0].tolist()
        y = D[1].tolist()
        minx, maxx = [0.8*0.11, 0.14*1.2]
        md = max(y)
        md = y.index(md)
        md = x[md]
        if minx <= md and md <= maxx: clist[i], alist[i], llist[i] = clr, 1, 1
        else: clist[i], alist[i], llist[i] = '0.6', 0.3, 0.4

        D = get_kdens_choose_kernel(Dmlist[i], kernel)
        x = D[0].tolist()
        y = D[1].tolist()
        minx, maxx = [0.8*0.92, 0.99*1.2]
        md = max(y)
        md = y.index(md)
        md = x[md]
        if minx <= md and md <= maxx and clist[i] == clr: clist[i], alist[i], llist[i] = clr, 1, 1
        else: clist[i], alist[i], llist[i] = '0.6', 0.3, 0.4

        D = get_kdens_choose_kernel(Emlist[i], kernel)
        x = D[0].tolist()
        y = D[1].tolist()
        minx, maxx = [-0.23*1.2, -0.22*0.8]
        md = max(y)
        md = y.index(md)
        md = x[md]
        if minx <= md and md <= maxx and clist[i] == clr: clist[i], alist[i], llist[i] = clr, 1, 1
        else: clist[i], alist[i], llist[i] = '0.6', 0.3, 0.4

        D = get_kdens_choose_kernel(Smlist[i], kernel)
        x = D[0].tolist()
        y = D[1].tolist()
        minx, maxx = [0.8*0.24, 0.38*1.2]
        md = max(y)
        md = y.index(md)
        md = x[md]
        if minx <= md and md <= maxx and clist[i] == clr: clist[i], alist[i], llist[i] = clr, 1, 1
        else: clist[i], alist[i], llist[i] = '0.6', 0.3, 0.4

    ns = [1,2,3,4]
    for i, n in enumerate(ns):
        if n == 1:
            fig.add_subplot(2, 2, n)
            for j, ls in enumerate(Rmlist):
                D = get_kdens_choose_kernel(Rmlist[j], kernel)
                plt.plot(D[0], D[1], color = clist[j], lw=llist[j], alpha=alist[j])

            plt.axvline(0.14, color='k', ls='--', lw = 1)
            plt.axvline(0.11, color='0.5', ls='--', lw = 1)
            #plt.xlim(0.05, 0.2)

        elif n == 2:
            fig.add_subplot(2, 2, n)
            for j, ls in enumerate(Dmlist):
                D = get_kdens_choose_kernel(Dmlist[j], kernel)
                plt.plot(D[0], D[1], color = clist[j], lw=llist[j], alpha=alist[j])

            plt.axvline(0.92, color='k', ls='--', lw = 1)
            plt.axvline(0.99, color='0.5', ls='--', lw = 1)
            #plt.xlim(0.9, 1.075)
            #plt.ylim(0, 40)

        elif n == 3:
            fig.add_subplot(2, 2, n)
            for j, ls in enumerate(Emlist):
                D = get_kdens_choose_kernel(Emlist[j], kernel)
                plt.plot(D[0], D[1], color = clist[j], lw=llist[j], alpha=alist[j])

            plt.axvline(-0.23, color='k', ls='--', lw = 1)
            plt.axvline(-0.21, color='0.5', ls='--', lw = 1)
            #plt.ylim(0, 30)
            #plt.xlim(-0.3, -0.125)

        elif n == 4:
            fig.add_subplot(2, 2, n)
            for j, ls in enumerate(Smlist):
                D = get_kdens_choose_kernel(Smlist[j], kernel)
                plt.plot(D[0], D[1], color = clist[j], lw=llist[j], alpha=alist[j])

            plt.axvline(0.38, color='k', ls='--', lw = 1)
            plt.axvline(0.24, color='0.5', ls='--', lw = 1)
            #plt.xlim(0.1, 0.4)

        plt.xlabel(xlabs[n-1], fontsize=8)
        plt.ylabel(ylab, fontsize=8)
        plt.tick_params(axis='both', labelsize=5)

    return fig

df_w = pd.read_csv(mydir + '/results/simulated_data/SimData-Random.csv')
df_s = pd.read_csv(mydir + '/results/simulated_data/SimData-Random-simple.csv')

df = pd.concat([df_w,df_s], axis=0)
#print list(df)
#sys.exit()

df2 = pd.DataFrame({'sim' : df['sim']})
df2['N'] = np.log10(df['total.abundance'])
df2['D'] = np.log10(df['N.max'])
df2['S'] = np.log10(df['species.richness'])
df2['E'] = np.log10(df['simpson.e'])
df2['R'] = df['logmod.skew']
df2['p'] = df['p']
df2['choose1'] = df['choose1']
df2['choose2'] = df['choose2']
df2['choose3'] = df['choose3']
df2['sim'] = df['sim']
df2['color'] = [s.replace('\'', '') for s in df['clr']]

fig = plt.figure()
xlabs = ['Scaling exponent for rarity', 'Scaling exponent for dominance',
    'Scaling exponent for evenness', 'Scaling exponent for $S$']
ylab = 'kernel density'
fig = figplot(df2, ['R','D', 'E', 'S'], xlabs, ylab, fig)

#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)

#plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Random.png', dpi=400, bbox_inches = "tight")
#plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Simple.png', dpi=400, bbox_inches = "tight")
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Combo.png', dpi=400, bbox_inches = "tight")

plt.close()
