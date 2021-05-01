import numpy as np
import matplotlib.pyplot as plt
import tellurium as te #
import roadrunner #
from rrplugins import Plugin
auto = Plugin("tel_auto2000")
import re
import seaborn as sns
import os
import te_bifurcation as bf
import matplotlib as mpl
mpl.use('GTK3Agg')
from sys import exit
from matplotlib.ticker import NullFormatter

model_mmi3 = {
    'pars':{
        'mu'    : 1,
        'b1'    : 1,
        'b2'    : 1,
        'b3'    : 1,
        'kR'    : 1,
        'kr'    : 1,
        'K1'    : 10000,
        'K2'    : 10000,
        'K3'    : 10000,
        'a1'    : 1,
        'a2'    : 1,
        'a3'    : 1,
        'Ts'    : 10,
    },
    'vars':{
        'r' : \
            'mu - kr * (r - 3 * c1 - 6 * c2 - 3 * c3) \
            - b1 * 1 * 3 * kr * c1 \
            - b2 * 2 * 3 * kr * c2 \
            - b3 * 3 * 1 * kr * c3',
        'R' : \
            '1 - kR * (R - 3 * c1 - 3 * c2 - c3) \
            - 3 * kR * c1 * a1 - 3 * kR * c2 * a2 -kR * c3 * a3',
        'c1':\
            'Ts * (K1 * (R - 3 * c1 - 3 * c2 - c3) * (r - 3 * c1 - 6 * c2 - 3 * c3) - c1)',
        'c2':\
            'Ts * (K2 * c1 * (r - 3 * c1 - 6 * c2 - 3 * c3) - c2)',
        'c3':\
            'Ts * (K3 * c2 * (r - 3 * c1 - 6 * c2 - 3 * c3) - c3)',
    },

'fns': {}, 'aux': [], 'name':'mmi3'}

ics_1 = {'r': 0.8, 'R': 0.1, 'c1': 0.0, 'c2': 0.0,  'c3': 0.0}

r = bf.model2te(model_mmi3, ics=ics_1)

K = 1E5
pars={'b1':1, 'b2':1, 'b3':1, 'a1':1, 'a2':1, 'a3':1, 'mu':1, 'sR':1, 'kr':1,'kR':1,
        'K1':K, 'K2': K, 'K3':K}
for k in pars:
    r[k] = pars[k]

m = r.simulate(0, 1000, 1000)


pars_base={'b1':1, 'b2':1, 'b3':1, 'a1':1, 'a2':1, 'a3':1}

np.random.seed(1)
dp = np.random.uniform(0.125, 16, size=(1000, 6))

a, b, c, d, e, f = dp[:,3], dp[:,4], dp[:,5], dp[:,0], dp[:,1], dp[:,2]

err, fn = 0, 0
bs_fcs, bs_ranges, ids_fn, data2plot, bounds2plot, i2plot = [], [], [], [], [], []
data_all = []
inuerr = []

for i in range(dp.shape[0]):
    fcs = dp[i,:]
    for ii, k in enumerate(('b1', 'b2', 'b3', 'a1', 'a2', 'a3')):
        r[k] = pars_base[k] * fcs[ii]
    data, bounds = bf.run_bf(r, auto, dirc="+", par="mu", lims=[0,15], dsmin=1E-5, dsmax=0.02, nmx=30000)
    data_all.append(data)
    if len(data) == 0:
        inuerr.append(i)
        continue
    if len(bounds) < 3: # if it is bistable, then the length of data should be four
        err += 1
        continue
    data2plot.append(data)
    bounds2plot.append(bounds)
    i2plot.append(i)
    bs_range = data[bounds[1]-1, 4] - data[bounds[2]-1, 4]  # The range of bistability
    bs_ranges.append(bs_range)

print('No Bistability: ', err, 'FN: ', fn)

ptt = len(data2plot)
tt = dp.shape[0] - len(inuerr)


bs_ranges = np.array(bs_ranges)
i_sel = i2plot[np.argsort(bs_ranges)[-1]]


cmap = plt.cm.RdYlGn

def plot_bif(var='Free mRNA', ax=None):
    for data, bounds, ip in zip(data2plot, bounds2plot, i2plot):
        j = 0
        for k, n in enumerate(bounds[1:]):
            st, en = bounds[k], n
            if var == 'Free mRNA':
                a5_arr = data[st:en, 6] - data[st:en,7] * 3 - data[st:en,8] * 3 - data[st:en,9]
            elif var == 'Total mRNA':
                a5_arr = data[st:en, 6]
            elif var == 'Free miRNA':
                a5_arr = data[st:en, 10] - data[st:en,7] * 3 - data[st:en,8] * 3 * 2 - data[st:en,9] * 3
            elif var == 'Total miRNA':
                a5_arr = data[st:en, 10]
            elif var == '1:2 complex':
                a5_arr = data[st:en, 7] * 3
            elif var == '1:1 complex':
                a5_arr = data[st:en, 8] * 3
            elif var == '1:3 complex':
                a5_arr = data[st:en, 9] 
            # 11:mir   6:P   7:mRNA  8:mmR1  9:mmR2  10:mmR3
            x_arr = data[st:en, 4]
            al, zorder = 1, 1
            lw = 0.15
            if data[int((st+en)/2), 1] > 0:
                stab = 'u'
                ls = '-'
                co = 'lightblue'
            else:
                stab = 's'
                ls = '-'
                if j == 1:
                    co = '0.1'
                elif j == 0:
                    co = cmap(0)
                j += 1
            if ip != i_sel:
                al = 0.7
                pass
            else:
                al = 1
                lw = lw*40
                zorder = 10
            ax.plot(x_arr, a5_arr, ls=ls, c=co, alpha=al, lw=lw, zorder=zorder)

    ax.set_yscale('symlog', linthreshy=2E-3)
    ax.set_yticks([0.001, 0.01, 0.1, 1])
    ax.set_xlabel('miRNA synthesis rate ' + r'$\it{s}_{r}$', size=12)
    ax.set_ylabel(var +' (A.U.)', size=12)
    ax.set_title(var, size=12, rotation=8)
    ax.set_ylim([-0.00005, 3])
    ax.set_xlim([0, 15])

# Bifurcation diagrams for free mRNA
fig, ax = plt.subplots(figsize=(3,3))
fig.subplots_adjust(bottom=0.17, left=0.28, right=0.97)
plot_bif(var='Free mRNA', ax=ax)
ax.set_title('')
ax.set_ylabel('Free mRNA\nconcentration (A.U.)')
ax.set_yscale('symlog', linthreshy=1E-4)


# Bifurcation diagrams for four other variables
fig, axs = plt.subplots(figsize=(6,3), ncols=5)
fig.subplots_adjust(bottom=0.17, left=0.18, right=0.97, wspace=0.1)
plot_bif(var='Free mRNA', ax=axs[0])
plot_bif(var='Total mRNA', ax=axs[1])
plot_bif(var='1:1 complex', ax=axs[2])
plot_bif(var='1:2 complex', ax=axs[3])
plot_bif(var='1:3 complex', ax=axs[4])
for ax in axs[1:]:
    ax.set_yticklabels([])
    ax.set_ylabel('')
axs[0].set_ylabel('Concentration (A.U.)')
for ax in axs[[0,1,3,4]]:
    ax.set_xlabel('')
for ax in axs:
    ax.set_xlim([-0.2, 15])
    ax.set_ylim([-0.00005, 10])
    ax.set_yscale('symlog', linthreshy=1E-4)
plt.show()



# Scatter plot
bsdata = np.array(bs_ranges)
idsp = np.where(bsdata>0)[0]
ids0 = np.array(list(set(range(len(dp))) - set(i2plot)))
bsdata[idsp] = np.log(bsdata[idsp])

def id2idp(ids):
    return [i2plot.index(x) for x in ids]

fig, axs = plt.subplots(figsize=(8,3), ncols=3)
fig.subplots_adjust(right=0.98, left=0.2, bottom=0.23)

cax = fig.add_axes([0.1, 0.2, 0.01, 0.7])

s2 = 2
axs[0].scatter((a/d)[ids0], (b/e/1)[ids0], c='gray', s=s2, lw=0)
sca = axs[0].scatter((a/d)[i2plot], (b/e/1)[i2plot], c=bsdata[id2idp(i2plot)], alpha=0.8, s=5, lw=0)

axs[1].scatter((a/d)[ids0], (c/f/1)[ids0], c='gray', s=s2, lw=0)
axs[1].scatter((a/d)[i2plot], (c/f/1)[i2plot], c=bsdata[id2idp(i2plot)], alpha=0.8, s=5, lw=0)

axs[2].scatter((b/e)[ids0], (c/f/1)[ids0], c='gray', s=s2, lw=0)
axs[2].scatter((b/e)[i2plot], (c/f/1)[i2plot], c=bsdata[id2idp(i2plot)], alpha=0.8, s=5, lw=0)

fig.colorbar(sca, cax=cax, orientation='vertical')
cax.set_ylabel('Range of bistability\nlog(miRNA synthesis rate ' + r'$\it{s}_{r}$)')
cax.yaxis.set_ticks_position('left')
cax.yaxis.set_label_position('left')

for ax in axs:
    ax.scatter(1, 1, c='orange',s=10, zorder=10)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.01, 100)
    ax.set_ylim(0.01, 100)
    ax.minorticks_off()
    ax.set_yticks([0.01, 1, 100])
for ax in axs[1:]:
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r'$\it{a_3/b_3}$', size=12, labelpad=-4)
axs[0].set_xlabel(r'$\it{a_1/b_1}$', size=12)
axs[0].set_ylabel(r'$\it{a_2/b_2}$', size=12, labelpad=-4)
axs[1].set_xlabel(r'$\it{a_1/b_1}$', size=12)
axs[2].set_xlabel(r'$\it{a_2/b_2}$', size=12)
plt.show()


