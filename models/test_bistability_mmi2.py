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

model_mmi2 = {
    'pars':{
        'mu'    : 1,
        'b1'    : 1,
        'b2'    : 1,
        'sR'    : 1,
        'kR'    : 1,
        'kr'    : 1,
        'K1'    : 100000,
        'K2'    : 100000,
        'a1'    : 1,
        'a2'    : 1,
        'Ts'    : 10,
    },
    'vars':{
        'r' : \
            'mu - kr * (r - 1 * 2 * c1 - 2 * 1 * c2) \
            - b1 * 1 * 2 * kr * c1 \
            - b2 * 2 * 1 * kr * c2',
        'R' : \
            'sR - kR * (R - 2 * c1 - 1 * c2) \
            - 2 * kR * c1 * a1 - 1 * kR * c2 * a2',
        'c1':\
            'Ts * (K1 * (R - 2 * c1 - c2) * (r - 2 * c1 - 2 * c2) - c1) - c1',
        'c2':\
            'Ts * (K2 * c1 * (r - 2 * c1 - 2 * c2) - c2) - c2',
    },

'fns': {}, 'aux': [], 'name':'mmi2'}


ics_1 = {'r': 0.8, 'R': 0.1, 'c1': 0.0, 'c2': 0.0}


r = bf.model2te(model_mmi2, ics=ics_1)

# set all key parameters
# The first four will be adjusted later

pars={'b1':1, 'b2':1, 'a1':1, 'a2':1, 'sR':1,'mu':1, 'kr':1.0, 'kR':1}
for k in pars:
    r[k] = pars[k]

m  = r.simulate(0, 1000, 1000)

pars_base = {'b1':1, 'b2':1, 'a1':1, 'a2':1}

np.random.seed(1)
dp = np.random.uniform(0.125, 16, size=(1000, 4))

a, b, c, d = dp[:,2], dp[:,3], dp[:,0], dp[:,1]

err, fn = 0, 0
bs_fcs, bs_ranges, ids_fn, data2plot, bounds2plot, i2plot = [], [], [], [], [], []
inuerr = []

for i in range(dp.shape[0]):
    fcs = dp[i,:]
    det = (fcs[2]/fcs[0]-fcs[3]/fcs[1]/2)
    for ii, k in enumerate(('b1', 'b2', 'a1', 'a2')):
        r[k] = pars_base[k] * fcs[ii]
    data, bounds, dummy = bf.run_bf(r, auto, dirc="+", par="mu", lims=[0,15], dsmin=1E-4, dsmax=0.1)
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


dets_fp, dets = [], []
tp, fn, fp, tn = 0, 0, 0, 0
hpsn = 0
for i, fcs in enumerate(dp[:]):
    if i in inuerr:
        continue
    det = (fcs[2]/fcs[0]-fcs[3]/fcs[1]/2)
    dets.append(det)
    if det < 0: # Among the positives
        if not i in i2plot: # False positives
            fp += 1
            dets_fp.append(det)
        else: # True positives
            tp += 1
            pt_tys = data2plot[tp-1][:,2]
            if sum(pt_tys==3) > 0:
                hpsn += 1
    else: # Among the negatives
        if not i in i2plot: # True negatives
            tn += 1
        else: # False negatives
            fn += 1


print('TP', 'FP', 'FN', 'TN', 'NumErr')
print(tp, fp, fn, tn, len(inuerr))
print(hpsn)


cmap = plt.cm.RdYlGn

def plot_bif(var='Free mRNA', ax=None):
    for data, bounds, ip in zip(data2plot, bounds2plot, i2plot):
        j = 0
        for k, n in enumerate(bounds[1:]):
            #print(i,n)
            st, en = bounds[k], n
            if var == 'Free mRNA':
                a5_arr = data[st:en, 6] - data[st:en,7] * 2 - data[st:en,8]
            elif var == 'Total mRNA':
                a5_arr = data[st:en, 6]
            elif var == 'Free miRNA':
                a5_arr = data[st:en, 9] - data[st:en,7] * 2 - data[st:en,8] * 2
            elif var == 'Total miRNA':
                a5_arr = data[st:en, 9]
            elif var == '1:2 complex':
                a5_arr = data[st:en, 8]
            elif var == '1:1 complex':
                a5_arr = data[st:en, 7] * 2
            # 10:mir   6:P   6:mRNA  8:mmR1  9:mmR2
            x_arr = data[st:en, 4]
            al, zorder = 1, 1
            al = 0.7
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
            if ip == 986:
                lw = 4
                al = 1
            ax.plot(x_arr, a5_arr, ls=ls, c=co, alpha=al, lw=lw, zorder=zorder)

    ax.set_yscale('symlog', linthreshy=2E-4)
    ax.set_yticks([0.0001, 0.001, 0.01, 0.1, 1])
    ax.set_xlabel('miRNA synthesis rate ' + r'$\it{s}_{r}$', size=12)
    ax.set_ylabel(var +' (A.U.)', size=12)
    ax.set_title(var, size=12)
    ax.set_ylim([-0.00005, 3])
    ax.set_xlim([0, 5])


# Bifurcation diagrams for free mRNA
fig, ax = plt.subplots(figsize=(3,3))
fig.subplots_adjust(bottom=0.17, left=0.30, right=0.97)
plot_bif(var='Free mRNA', ax=ax)
ax.set_title('')
ax.set_ylabel('Free mRNA\nconcentration (A.U.)');


# Bifurcation diagrams for three other variables
fig, axs = plt.subplots(figsize=(4,3), ncols=3)
fig.subplots_adjust(bottom=0.17, left=0.18, right=0.97, wspace=0.1)
plot_bif(var='Total mRNA', ax=axs[0])
plot_bif(var='1:1 complex', ax=axs[1])
plot_bif(var='1:2 complex', ax=axs[2])
for ax in axs[1:]:
    ax.set_yticklabels([])
    ax.set_ylabel('')
axs[0].set_ylabel('Concentration (A.U.)')
for ax in axs[[0,2]]:
    ax.set_xlabel('')
for ax in axs:
    ax.set_xlim([-0.1, 5])
    ax.set_ylim([0, 10])

# Scatter plot
bsdata = np.array(bs_ranges)
ids_tr = [ i for i in range(dp.shape[0]) if not i in ids_fn]
idsp = np.where(bsdata>0)[0]
ids0 = np.array(list(set(range(len(a))) - set(i2plot) - set(inuerr)))
bsdata[idsp] = np.log(bsdata[idsp])

fig, ax = plt.subplots(figsize=(4,3))
fig.subplots_adjust(right=0.97, left=0.5, bottom=0.23)

cax = fig.add_axes([0.2, 0.2, 0.03, 0.7])

sca = ax.scatter((a/c)[i2plot], (b/d/1)[i2plot], c=bsdata, alpha=0.8, s=15, lw=0)
ax.scatter((a/c)[ids0], (b/d/1)[ids0], c='gray', s=15, lw=0)

fig.colorbar(sca, cax=cax, orientation='vertical')
cax.set_ylabel('Range of bistability\nlog(miRNA synthesis rate ' + r'$\it{s}_{r}$)')
cax.yaxis.set_ticks_position('left')
cax.yaxis.set_label_position('left')

ax.plot([0.01, 50], [0.02, 100], c='indianred')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.01, 100)
ax.set_ylim(0.01, 100)
ax.set_xlabel(r'$a_1/b_1$', size=12)
ax.set_ylabel(r'$a_2/b_2$', size=12)
ax.minorticks_off()
ax.set_yticks([0.01, 1, 100])
plt.show()
