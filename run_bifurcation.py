import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
import re
import numpy as np
from rrplugins import Plugin
auto = Plugin("tel_auto2000")
import te_bifurcation as bf

jl_file = "./models/model_mmiS.jl"
model = bf.load_jl(jl_file)
ics_1 = {"R8": 0.00900023576123889, "c15": 1.761337400332505e-06, "P8": 2.9560034347325604e-06, "P5": 3.4020485285993286e-07, "c25": 9.118941517979407e-05, "c35": 0.004721137474792935, "A": 1.5e-323, "c18": 8.844673354103184e-06, "R5": 0.000500023861019239, "c38": 0.001867111128611335, "F": 1.5e-323, "r2": 0.3524443000060298, "r9": 0.33001432105772416, "c28": 0.0002687664938900425, "c48":0.001, 'A':0.001, 'F':0.001}
r = bf.model2te(model=model, ics=ics_1)

r.x = 50


data_all, bounds_all, ras, fgfs = [], [], [], []

r.x = 10
r.maxA = r.maxA *0.65 #/0.65 # 0.1
r.maxF = r.maxF *0.65 #/0.65 # 0.6
m = r.simulate (0, 400, 100)

for k in ics_1.keys():
    if k == 'r2' or k == 'r9':
        ics_1[k] = 1
    else:
        ics_1[k] = 0

if 1:
    r.l05 = 0.01
    r.l15 = 0.5
    r.l25 = 0.1
    d1, b1 = bf.run_bf(r,auto,dirc="-",par="x",lims=[-100, 200])
    for k, v in ics_1.items():
        r[k] = v
    d2, b2 = bf.run_bf(r,auto,dirc="-",par="x",lims=[0, 100])
    d3, b3 = bf.run_bf(r,auto,dirc="+",par="x",lims=[0, 300])
    fig = bf.plot_bfdata58([d1, d2, d3], [b1, b2, b3])
    plt.show()
#fig.savefig('./figures/frames_bf.svg', format='svg')

r.x = 0
r.maxA = 0
r.l05 = 0.03 # 1.8
r.l15 = 0.5 # 0
r.K5A = 7.2 * 0.1
r.n5A = 2.0 * 3
r.K2A = 3.2 * 0.2
d, b = bf.run_bf(r,auto,dirc="+",par="maxA",lims=[0, 10],dsmax=1E-2, dsmin=1E-5, ds=1E-2)
fig2 = bf.plot_bfdata58([d], [b])

r.s2 = 0
dko, bko = bf.run_bf(r,auto,dirc="+",par="maxA",lims=[0, 10],dsmax=1E-2, dsmin=1E-5, ds=1E-2)
bf.plot_bfdata58([dko], [bko], color='purple', refig=fig2)
for ax in fig2.axes:
    ax.set_xlim(0.2, 1.0)
    ax.set_yscale('symlog', linthreshy=1E-2)
    ax.set_xlabel("RA (a.u.)")
    ax.set_xticks([0.2, 0.6, 1.0])
    ax.set_xticklabels([0.2, 0.6, 1.0])
fig2.axes[0].set_ylabel("Hoxa5 steady\nstate (a.u.)")
fig2.savefig('./figures/bf_ra_wtvs27ko.svg', format='svg')
plt.show()
