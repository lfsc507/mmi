# mmi
Code for research article 

<a href="http://dx.doi.org/10.15252/msb.20209945">__"MicroRNA governs bistable cell differentiation and lineage segregation via a noncanonical feedback"<br>
Li C-J, Liau ES, Lee Y-H, Huang Y-Z, Liu Z, Willems A, Garside V, McGlinn E, Chen J-A, Hong T<br>
Mol Syst Biol (2021) 17:e9945__</a>


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
PyCall 1.91.4<br>


__Simulation of mmi-S Model:__<br>
For steady state distributions and movies, run `run_simulation.jl`.<br>


__Bifurcation diagrams:__<br>
For mmi-2 and mmi-3 Models with sampled parameters, run `test_bistability_mmi2.py` and `test_bistability_mmi3.py` under `models`. <br>
For mmi-S Model and its variants, run `run_bifurcation_RA.py`.


__Phaseplanes:__<br>
XPPAUT code included in `nullclines.ode`.
Note: parameter values of &mu for 'Cooperative mRNA degradation', 'Cooperative miRNA degradation' and 'TDMD' are 0.3, 0.4 and 1.2 respectively.
