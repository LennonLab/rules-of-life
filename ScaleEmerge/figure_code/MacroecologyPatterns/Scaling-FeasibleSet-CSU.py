from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser('~/GitHub/ScaleEmerge')
sys.path.append(mydir+'/tools')
mydir2 = os.path.expanduser("~/")

def e_simpson(sad):
    sad = filter(lambda a: a != 0, sad)
    D = 0.0
    N = sum(sad)
    S = len(sad)
    for x in sad: D += (x*x) / (N*N)
    E = round((1.0/D)/S, 4)
    return E


def figplot(x, y, xlab, ylab, fig, n):

    '''main figure plotting function'''

    fig.add_subplot(2, 2, n)
    plt.style.use('dark_background')

    d = pd.DataFrame({'x': list(x)})
    d['y'] = list(y)
    f = smf.ols('y ~ x', d).fit()

    m, b, r, p, std_err = stats.linregress(x, y)

    st, data, ss2 = summary_table(f, alpha=0.05)
    fitted = data[:,2]
    mean_ci_low, mean_ci_upp = data[:,4:6].T
    ci_low, ci_upp = data[:,6:8].T

    x, y, fitted, ci_low, ci_upp = zip(*sorted(zip(x, y, fitted, ci_low, ci_upp)))

    if n == 1:
        lab = r'$Rarity$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        plt.text(1.1, 0.9, lab, fontsize=8)

    elif n == 2:
        lab = r'$Dominance$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        plt.text(1.1, 2.2, lab, fontsize=8)

    elif n == 3:
        lab = r'$Evenness$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        plt.text(1.1, -1.2, lab, fontsize=8)

    elif n == 4:
        lab = r'$S$'+ ' = '+str(round(10**b,2))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'+'\n'
        plt.text(1.1, 2.1, lab, fontsize=8)

    #plt.hexbin(x, y, mincnt=1, gridsize = 100, bins='log', cmap=plt.cm.jet)
    plt.scatter(x, y, color = 'SkyBlue', alpha= 1 , s = 3, linewidths=0.25, edgecolor='Steelblue')

    if n == 3: plt.legend(loc='best', fontsize=6, frameon=False)

    plt.plot(x, fitted,  color='w', ls='--', lw=1.0, alpha=0.9)
    plt.xlabel(xlab, fontsize=9)
    plt.ylabel(ylab, fontsize=10)
    plt.tick_params(axis='both', labelsize=5)
    plt.xlim(1, 1.05*max(x))

    if n == 1: plt.ylim(0.0, 1.1)
    elif n == 2: plt.ylim(0.2, 2.7)
    elif n == 3: plt.ylim(-1.2, 0.05)
    elif n == 4: plt.ylim(0.4, 2.5)

    return fig


GenPath = mydir + "/results/simulated_data/"
dat = open(GenPath + 'partitions.csv', 'r')

Ns, Rs, Es, Ds, Ss = [], [], [], [], []
for line in dat:
    partition = eval(line)

    N = np.log10(sum(partition))
    D = np.log10(max(partition))
    S = np.log10(len(partition))
    E = np.log10(e_simpson(partition))
    R = stats.skew(partition)
    R = np.log10(abs(float(R)) + 1)
    if R < 0: R = R * -1

    Ns.append(N)
    Ds.append(D)
    Ss.append(S)
    Es.append(E)
    Rs.append(R)


fig = plt.figure()
xlab = '$N$, $log$'+r'$_{10}$'
ylab = 'Rarity, '+r'$log_{10}$'
fig = figplot(Ns, Rs, xlab, ylab, fig, 1)

ylab = 'Dominance, '+r'$log_{10}$'
fig = figplot(Ns, Ds, xlab, ylab, fig, 2)

ylab = 'Evenness, ' +r'$log_{10}$'
fig = figplot(Ns, Es, xlab, ylab, fig, 3)

ylab = 'Richness, ' +r'$log_{10}$'
fig = figplot(Ns, Ss, xlab, ylab, fig, 4)


#### Final Format and Save #####################################################
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/results/figures/DiversityAbundanceScaling-Partitions.png',
    dpi=400, bbox_inches = "tight")
plt.close()
