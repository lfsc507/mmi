# mmi
Code for research article 

<a href="http://dx.doi.org/10.15252/msb.20209945">__"MicroRNAs govern bistable cell differentiation and lineage segregation via noncanonical feedback"<br>
Li C-J, Liau ES, Lee Y-H, Huang Y-Z, Liu Z, Willems A, Garside V, McGlinn E, Chen J-A, Hong T<br>
Mol Syst Biol (2021) 17:e9945__</a>

If you find the code helpful for your research, please cite the paper.
<a href="https://scholar.googleusercontent.com/scholar.bib?q=info:c8UoTMnAt4AJ:scholar.google.com/&output=citation&scisdr=CgVBh267EOOI5vol5d0:AAGBfm0AAAAAYmgj_d0onqU_2D2z-VtAbl7KPdDP3hTm&scisig=AAGBfm0AAAAAYmgj_RLCtsMDvEuBh0gnEcj-ftN0BAbc&scisf=4&ct=citation&cd=-1&hl=en">BibTex</a>
<a href="https://scholar.googleusercontent.com/scholar.enw?q=info:c8UoTMnAt4AJ:scholar.google.com/&output=citation&scisdr=CgWWDS1BEPDG6k3vMWw:AAGBfm0AAAAAYLvqKWwBCv_OcL9ymEQQmpjk8RBxKGU_&scisig=AAGBfm0AAAAAYLvqKVo7MQLsGXzPLzxamZxEv9ragJW7&scisf=3&ct=citation&cd=-1&hl=en&scfhb=1">EndNote</a>
<a href="https://scholar.googleusercontent.com/scholar.ris?q=info:c8UoTMnAt4AJ:scholar.google.com/&output=citation&scisdr=CgWWDS1BEPDG6k3vMWw:AAGBfm0AAAAAYLvqKWwBCv_OcL9ymEQQmpjk8RBxKGU_&scisig=AAGBfm0AAAAAYLvqKVo7MQLsGXzPLzxamZxEv9ragJW7&scisf=2&ct=citation&cd=-1&hl=en&scfhb=1">RefMan</a>
<a href="https://scholar.googleusercontent.com/scholar.rfw?q=info:c8UoTMnAt4AJ:scholar.google.com/&output=citation&scisdr=CgWWDS1BEPDG6k3vMWw:AAGBfm0AAAAAYLvqKWwBCv_OcL9ymEQQmpjk8RBxKGU_&scisig=AAGBfm0AAAAAYLvqKVo7MQLsGXzPLzxamZxEv9ragJW7&scisf=1&ct=citation&cd=-1&hl=en&scfhb=1">RefWorks</a>


__Requirements:__<br>
Code has been tested with:<br>
Python 3.6<br>
Tellurium 2.1.5 (rrplugins 2.1.1)<br>
Numpy 1.18.1<br>
Matplotlib 3.2.0<br>
Seaborn 0.9.0<br><br>
Julia 1.4.1<br>
DifferentialEquations 6.14.0<br>
Distributions 0.23.4<br>
PyCall 1.91.4<br><br>
XPP-AUT 6.11

__Simulation of mmi-S Model, or other spatial models:__<br>
For steady state distributions and movies, run `run_simulation.jl`.<br>


__Bifurcation diagrams:__<br>
For mmi-2 and mmi-3 Models with sampled parameters, run `test_bistability_mmi2.py` and `test_bistability_mmi3.py` under `models`. <br>
For mmi-S Model and its variants, run `run_bifurcation_RA.py`.


__Phaseplanes:__<br>
XPP-AUT code is included in `nullclines.ode`.<br>
Note: parameter values of &#956; (miRNA production rate constant) for 'Cooperative mRNA degradation', 'Cooperative miRNA degradation' and 'TDMD' are 0.3, 0.4 and 1.2 respectively (information missing from figure legend).
