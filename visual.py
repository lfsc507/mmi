# This is a visualization utility file that should be under
# a python path that julia can find
import numpy as np
import matplotlib.animation as manimation
import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
import matplotlib as mpl
from sys import exit
import pandas as pd
import seaborn as sns
mpl.rcParams["font.family"] = "sans-serif"

def anim(d, svlist, if_show=True, if_anim=True, anim_file=None,
        fig_name = './figures/pat_all_genes.svg', fmt='svg', fig=None, axs=None):
    if fig:
        axs2 = np.array(axs)
        for ax in axs2.flatten():
            ax.clear()
    else:
        fig, axs2 = plt.subplots(figsize=(3,2), ncols=2, nrows=4)
        fig.subplots_adjust(hspace=0.15, left=0.08, bottom=0.1, right=0.65, wspace=0.35)
    for ax in [axs2[3,0]]:
        ax.axis("off")

    d = np.swapaxes(np.swapaxes(d, 1, 3), 0, 1)

    vars_sorted = svlist
    d5 = d[:, vars_sorted.index('P5'),:,:]
    d8 = d[:, vars_sorted.index('P8'),:,:]

    nrows, ncols = d5[0].shape[0], d5[0].shape[1]

    cax = fig.add_axes([0.7, 0.2, 0.02, 0.3])
    caxg = fig.add_axes([0.7, 0.6, 0.02, 0.3])
    ax2 = axs2[3,1]

    d5ln, d8ln = np.log(d5+0.011), np.log(d8+0.011)
    d5ln = (d5ln - d5ln.min()) / (d5ln.max() - d5ln.min())
    d8ln = (d8ln - d8ln.min()) / (d8ln.max() - d8ln.min())
    sums = np.clip(d8ln+d5ln, 0, 1)

    ratios = (d8ln)/(d5ln+1)

    sums5 = np.clip(d5ln, 0, 1)
    ratios5 = d5ln

    sums8 = np.clip(d8ln, 0, 1)
    ratios8 = d8ln

    cmap = plt.cm.RdYlGn
    ratiosc, ratiosc5, ratiosc8 = cmap(ratios), cmap(ratios5), cmap(ratios8)
    ratiosc[..., -1] = sums
    ratiosc5[..., -1] = sums5
    ratiosc8[..., -1] = sums8

    greys = np.empty(d5.shape + (3,), dtype=np.uint8)
    greys.fill(100)

    dR5 = d[:, vars_sorted.index('R5'),:,:]
    dR8 = d[:, vars_sorted.index('R8'),:,:]
    d2 = d[:, vars_sorted.index('r2'),:,:]
    d9 = d[:, vars_sorted.index('r9'),:,:]

    tp = -1
    cmap = plt.cm.Greys
    for d, ax in zip([dR5, dR8, d2, d9, d5, d8],\
            [axs2[0,0], axs2[0,1], axs2[1,0], axs2[1,1], axs2[2,0], axs2[2,1]]):
        d_lg = np.log(d+0.0001)
        d_no = d_lg
        im = ax.imshow(d_no[tp], cmap=cmap, aspect='auto', vmax=-1.0, vmin=-5)
        ax.set_xticks([])
        ax.set_yticks([])
    axs2[2,0].set_xlabel('x/L', size=8, fontname='sans-serif')
    axs2[0,0].set_ylabel('Total a5\nmRNA', size=7, labelpad=0, fontname="sans-serif")
    axs2[1,0].set_ylabel('Total\nmiR-27', size=7, labelpad=0, fontname="sans-serif")
    axs2[0,1].set_ylabel('Total c8\nmRNA', size=7, labelpad=0, fontname="sans-serif")
    axs2[1,1].set_ylabel('Total\nmiR-196', size=7, labelpad=0, fontname="sans-serif")
    axs2[2,0].set_ylabel('a5\nprotein', size=7, labelpad=0, fontname="sans-serif")
    axs2[2,1].set_ylabel('c8\nprotein', size=7, labelpad=0, fontname="sans-serif")
    cbar = fig.colorbar(im, cax=caxg, orientation='vertical', cmap=cmap)
    cbar.set_ticks([-1, -3, -5])
    cbar.set_ticklabels([r'$10^{-1}$', r'$10^{-3}$', r'$10^{-5}$'])
    cbar.set_label('Concentration\n(A.U.)', size=7, fontname='sans-serif')
    for l in cbar.ax.get_yticklabels():
        l.set_fontname("sans-serif")
    caxg.tick_params(axis='both', which='major', labelsize=8)
    ratios2plot = ratiosc

    ax2.tick_params(labelbottom=False)
    ax2.tick_params(labelleft=False)

    ax2.imshow(greys[0])

    rmin, rmax = ratios2plot.min(), ratios2plot.max()
    mock = ratios2plot[0].copy()
    mock[0,0], mock[-1,-1] = rmin, rmax
    cmap = plt.cm.RdYlGn
    im = ax2.imshow(mock, cmap=cmap, aspect='auto')
    tx = ax2.text(16, 13, 'x/L', size=8, fontname='sans-serif')
    ax2.set_xticks([])
    ax2.set_yticks([])
    cbar = fig.colorbar(im, cax=cax, orientation='vertical', cmap=cmap)
    cbar.ax.set_yticklabels([r'$a5^{\rm on}c8^{\rm off}$', 'Hybrid', r'$a5^{\rm off}c8^{\rm on}$'])
    for l in cbar.ax.get_yticklabels():
        l.set_fontname("monospace")
    cax.tick_params(axis='both', which='major', labelsize=8)
    if if_anim == False:
        ax2.imshow(ratios2plot[tp], aspect='auto')
        fig.savefig(fig_name, format=fmt, dpi=600)
        if if_show == True:
            plt.show()

    if if_anim == True:
        def init():
            im.set_data(ratios2plot[0])
            tx.set_text('x/L' + ' '*2 +  'Time:' +\
                        "{:5}".format((1)*1) + '')
            return im, tx

        def update_im(n):
            im.set_data(ratios2plot[n+1])
            tx.set_text('x/L' + ' '*2 +  'Time:' +\
                    "{:5}".format((n+1)*1+1) + '')
            return im, tx

        ani = manimation.FuncAnimation(fig, update_im, init_func=init,\
                frames=len(ratios)-1,\
                interval=200, blit=False, repeat_delay=500)

        mpl.rcParams['animation.codec'] = 'mpeg4'
        mpl.rcParams['savefig.dpi'] = 900
        if anim_file:
            ani.save(anim_file)
        plt.show()

def plot_xL(d, svlist, vars2plot=[],
        fig_name = './figures/ssLine.svg',
        if_show=True, ax=None):
    d = np.swapaxes(np.swapaxes(d, 1, 3), 0, 1)
    if ax == None:
        fig, ax = plt.subplots(figsize=(3,2))
        fig.subplots_adjust(left=0.32, bottom=0.29)
    for var in vars2plot:
        iv = svlist.index(var)
        dfa = pd.DataFrame(d[-1,iv,:,:]).melt()
        sns.lineplot(x="variable", y="value", data=dfa, ax=ax, ci=95,
            err_style="bars", label=var,
            err_kws={"ecolor":"indianred","capsize":1})
    ax.set_yscale("symlog", linthreshy=1E-2)
    ax.set_ylabel("Expression\nlevel(a.u.)",size=12)
    ax.set_xticks([])
    ax.set_xlabel("x/L",size=12)
    #fig.savefig(fig_name, format="svg")
    if if_show == True:
        plt.show()
    return ax

