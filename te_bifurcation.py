# Bifurcation utility functions with Tellurium
# To be placed under a Python path
import tellurium as te
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
import re
import os
import sympy
import numpy as np

def extract_data():
    with open('fort.7', 'r') as f:
        lines = f.readlines()

    data = []
    for line in lines[12:]:
        if re.match(r'\s+0', line):
            break
        fs = re.split(r' +', line.strip())
        data.append([float(f) for f in fs])
    data = np.array(data)
    if len(data.shape) == 1:
        return [], []
    data = data[data[:,3]>0,:]
    idsn = np.where(data[:,1]<0)[0]
    idsp = np.where(data[:,1]>0)[0]

    bksn = np.where((idsn[1:]-idsn[:-1])>1)[0]
    bksp = np.where((idsp[1:]-idsp[:-1])>1)[0]
    bks = np.where((data[:,2]==2)|(data[:,2]==1))[0]

    bounds = [0]+list(bks)+[len(data)]

    with open('fort.8', 'r') as f:
        f_str = f.read()

    blks = re.split('\n +1 +.*\n', f_str)
    half_blk = int(blks[0].count('\n')/2)
    numlines = [re.split("\n",blk) for blk in blks]
    numlines[0] = numlines[0][1:]
    states = [[float(num)
        for num in re.split(' +', "".join(lines[:half_blk]).strip())[1:]
        ] for lines in numlines]
    data8 = np.array(states)[:,1:]
    data = np.hstack([data, data8])

    return data, bounds

def load_jl(jl_file):
    with open(jl_file) as f:
        lines = f.readlines()
    model = {'vars':{}, 'pars':{}, 'fns':{}, 'aux':[], 'name':'mmiS'}
    seen_par = 0
    for line in lines:
        if line.lstrip().startswith('du'):
            m = re.search("d(\w+) *= *(.*)", line)
            eq = re.sub("μ", "mu", m.group(2))
            model['vars'][m.group(1)] = eq
        if not seen_par:
            if re.match("p *= *Dict\(", line):
                seen_par = 1
            else:
                continue
        m = re.search("\"(\w+)\" *=> *(.*)$", line)
        par = m.group(1)
        pv = sympy.N(re.sub("\)", "", m.group(2)))
        model['pars'][par] = pv
        if re.search("\)", line):
            seen_par = 0
    for line in lines:
        if re.match(r"p\[\"\w+\"\] *= *", line):
            m = re.search(r"p\[\"(\w+)\"\] *= *(.*)$", line)
            par = m.group(1)
            pv = sympy.N(re.sub("\)", "", m.group(2)))
            model['pars'][par] = pv
    return model

def model2te(model, ics={}):
    model_str = '// Reactions\n\t'

    for i, var in enumerate(sorted(model['vars'], reverse=False)):
        de = model['vars'][var]
        model_str += 'J'+ str(i) + ': -> ' + var + '; ' + de + '\n\t'

    model_str += '\n// Species Init\n\t'

    for k, v in ics.items():
        model_str += k + ' = ' + str(round(v,4)) + '; '

    model_str += '\n\n// Parameters\n\t'

    for k, v in model['pars'].items():
        model_str += k + ' = ' + str(v) + '; '

    r = te.loada(model_str)
    return r

def run_bf(r, auto, dirc="Positive", par="", lims=[0, 1],
        ds=0.001, dsmin=1E-5, dsmax=1,
        pre_sim_dur=10, nmx=10000):
    if dirc.lower()[:3] == "pos" or dirc == "+":
        dirc = "Positive"
    elif dirc.lower()[:3] == "neg" or dirc == "-":
        dirc = "Negative"
    # Setup properties
    auto.setProperty("SBML", r.getCurrentSBML())
    auto.setProperty("ScanDirection", dirc)
    auto.setProperty("PrincipalContinuationParameter", par)
    auto.setProperty("PreSimulation", "True")
    auto.setProperty("PreSimulationDuration", pre_sim_dur)
    auto.setProperty("RL0", lims[0])
    auto.setProperty("RL1", lims[1])
    auto.setProperty("NMX", nmx)
    auto.setProperty("NPR", 2)
    auto.setProperty("KeepTempFiles", True)
    auto.setProperty("DS", ds)
    auto.setProperty("DSMIN", dsmin)
    auto.setProperty("DSMAX", dsmax)
    auto.execute()
    pts     = auto.BifurcationPoints
    #print(pts1)
    data, bounds = extract_data()
    if os.path.exists('fort.7'):
        os.remove('fort.7')
        os.remove('fort.8')
    return data, bounds

def plot_bfdata58(data, bounds, color=None, refig=None):
    data = np.vstack(data)
    bounds_all = bounds[0]
    for i in range(1, len(bounds)):
        bounds_all = bounds_all + [b+bounds_all[-1] for b in bounds[i]]

    cmap = plt.cm.RdYlGn
    cs = [cmap(0.999), 'gray', cmap(0)]
    if not refig:
        fig, (ax0, ax1) = plt.subplots(ncols=2, nrows=1, figsize=(3,1.3), sharex=True, sharey=True)
        fig.subplots_adjust(hspace=0.15, left=0.23, bottom=0.27, \
                top=0.8, right=0.75, wspace=0.25)
    else:
        ax0, ax1 = refig.axes[:2]

    for ax in (ax0, ax1):
        ax.set_xlim([10, 80])
        ax.set_yscale('symlog', linthreshy=1E-4)

    j = 0
    for i, n in enumerate(bounds_all[1:]):
        st, en = bounds_all[i], n
        a5_arr = data[st:en, 8]
        c8_arr = data[st:en, 18]
        x_arr = data[st:en, 4]
        zorder = 5
        if data[int((st+en)/2), 1] > 0:
            stab = 'u'
            ls = '--'
            c = 'lightblue'
            al = 0.9
            lw = 2
            zorder = 10
        else:
            stab = 's'
            ls = '-'
            lw = 4
            j += 1
            al = 0.9
            if a5_arr[-1] > 0.01 and c8_arr[-1] < 0.07:
                c = cmap(0)
            elif a5_arr[-1] < 0.01 and c8_arr[-1] > 0.07:
                c = cmap(0.999)
            elif a5_arr[-1] < 0.01 and c8_arr[-1] < 0.07:
                c = 'gray'
                lw = 8
                label = None
                zorder = 0
            elif a5_arr[-1] > 0.01 and c8_arr[-1] > 0.07:
                c = cmap(0.5)
                al=1
                lw = 8
                zorder = 0
            else:
                print('Unknown phenotype.')
                c = 'gray'
        #print(st, en, stab)
        if color:
            c = color
        ax0.plot(data[st:en, 4], data[st:en, 8], ls=ls, c=c, alpha=al, lw=lw, zorder=zorder)
        ax1.plot(data[st:en, 4], data[st:en, 18], ls=ls, c=c, alpha=al, lw=lw, zorder=zorder)
        x_arr[x_arr>80] = np.nan
        x_arr[x_arr<10] = np.nan

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    if not refig:
        cax = fig.add_axes([0.8, 0.1, 0.02, 0.35])
        cbar = fig.colorbar(sm, cax=cax, orientation='vertical', cmap=cmap)
        cbar.ax.set_yticklabels([r'$a5^+c8^-$', 'Hybrid', r'$a5^-c8^+$'])
        cax.tick_params(axis='both', which='major', labelsize=8)

    for ax in (ax0, ax1):
        ax.set_xlabel('x/L', size=8)
        ax.xaxis.set_ticks([10, (80-10)/2+10, 80])
        ax.xaxis.set_ticklabels([0, 0.5, 1], size=8)
        ax.tick_params(axis='both', which='major', labelsize=8)

    ax0.set_yticks([1E-1, 1E-3, 0])
    ax0.set_ylabel('Protein\nsteady state (A.U.)', size=8)
    #ax0.set_ylabel('Steady State\nHoxc8 (a.u.)', size=14)
    ax0.set_title('Hoxa5', size=9)
    ax1.set_title('Hoxc8', size=9)

    if not refig:
        return fig

