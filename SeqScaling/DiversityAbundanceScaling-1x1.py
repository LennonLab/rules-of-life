from __future__ import division
import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table


p, _lw, w, fs, sz = 2, 0.25, 1, 6, 8

mydir = os.path.expanduser('~/GitHub/SeqScaling')
df = pd.read_csv(mydir + '/SimData.csv')
#df = df[df['S'] > 1]

'rep,length,sample_size,mutation_rate,generations,N,S,E,Nmax,R'

fig = plt.figure()
plt.style.use('classic')
xlab = 'Abundance of sequences, ' + '$log$'+r'$_{10}$'
ylab = 'Diversity of sequences, ' +r'$log_{10}$'
x, y = np.log10(df['N']), np.log10(df['S'])
y2 = list(y)
x2 = list(x)

d = pd.DataFrame({'x': list(x2)})
d['y'] = list(y2)
f = smf.ols('y ~ x', d).fit()

m, b, r, p, std_err = stats.linregress(x2, y2)
st, data, ss2 = summary_table(f, alpha=0.05)
fitted = data[:,2]
mean_ci_low, mean_ci_upp = data[:,4:6].T
ci_low, ci_upp = data[:,6:8].T

x2, y2, fitted, ci_low, ci_upp = zip(*sorted(zip(x2, y2, fitted, ci_low, ci_upp)))

lbl = r'$S$'+ ' = '+str(round(10**b,1))+'*'+r'$N$'+'$^{'+str(round(m,2))+'}$'
plt.fill_between(x2, np.array(ci_upp)+0.05, np.array(ci_low)+0.05, color='0.8', lw=0.5, alpha=0.4)

plt.hexbin(x2, y2, mincnt=1, gridsize = 20, bins='log', cmap=plt.cm.jet)
#plt.scatter(x2, y2, s = 40, color='blue', linewidths=0.5, edgecolor='w', alpha=0.5)
plt.plot(x2, np.array(fitted)+0.05,  color='k', ls='--', lw=2, alpha=0.9)

plt.text(1.1, 1.0, lbl, fontsize=26)
plt.xlabel(xlab, fontsize=20)
plt.ylabel(ylab, fontsize=20)
plt.tick_params(axis='both', labelsize=14)
plt.xlim(1, 3.5)
plt.ylim(min(ci_low), max(ci_upp))
plt.ylim(0.21, 1.5)

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig(mydir + '/DiversityAbundanceScaling-1x1.png', dpi=200, bbox_inches = "tight")
plt.close()
