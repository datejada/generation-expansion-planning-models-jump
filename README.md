# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using Julia/JuMP

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datejada/generation-expansion-planning-models-jump/main)

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

- Stochastic-GEP-two-stage-nb.ipynb
- Stochastic-GEP-two-stage-explicit-nb.ipynb
- Stochastic-GEP-two-stage-Benders-nb.ipynb
- Stochastic-GEP-two-stage-Benders-multicut-nb.ipynb
- Stochastic-GEP-two-stage-LR-nb.ipynb
- Stochastic-GEP-multi-stage-nb.ipynb
- Static-robust-optimization-GEP-nb.ipynb
- Adaptive-robust-optimization-GEP-nb.ipynb

These examples show basic concepts for learning optimization under uncertainty in power systems.

The models are developed in [Julia](https://julialang.org/), using the package [JuMP](https://jump.dev/JuMP.jl/stable/), and solved using [HiGHS](https://highs.dev/).

The main references to model the optimization problems are:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)

[2] [A. J. Conejo, L. Baringo, S. J. Kazempour and A. S. Siddiqui, Investment in Electricity Generation and Transmission, Cham, Zug, Switzerland:Springer, 2016.](https://link.springer.com/book/10.1007/978-3-319-29501-5)

[3] [Sun X.A., Conejo A.J. (2021) Static Robust Optimization. In: Robust Optimization in Electric Energy Systems. International Series in Operations Research & Management Science, vol 313. Springer, Cham.]( https://doi.org/10.1007/978-3-030-85128-6_2)
