# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using Julia/JuMP

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datejada/generation-expansion-planning-models-jump/main)

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP-two-stage-nb.ipynb**: Two-Stage Stochastic GEP
* **Stochastic-GEP-two-stage-Benders-nb.ipynb**: Two-Stage Stochastic GEP solved using Benders' Decomposition
* **Stochastic-GEP-two-stage-Benders-multicut-nb.ipynb**: Two-Stage Stochastic GEP solved using Benders' Decomposition with multi-cuts
* **Stochastic-GEP-two-stage-LR-nb.ipynb**: Two-Stage Stochastic GEP solved using Lagrangian Relaxation
* **Stochastic-GEP-multi-stage-nb.ipynb**: Multi-Stage Stochastic GEP
* **Static-robust-optimization-GEP-nb.ipynb**: Static Robust Optimization (Scenario-based) GEP
* **Adaptive-robust-optimization-GEP-nb.ipynb**: Adaptive Robust Optimization (ARO) GEP

These examples show basic concepts for learning optimization under uncertainty in power systems.

The models are developed in [Julia](https://julialang.org/), using the package [JuMP](https://jump.dev/JuMP.jl/stable/), and solved using [HiGHS](https://highs.dev/).

The main references to model the optimization problems are:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)

[2] [A. J. Conejo, L. Baringo, S. J. Kazempour and A. S. Siddiqui, Investment in Electricity Generation and Transmission, Cham, Zug, Switzerland:Springer, 2016.](https://link.springer.com/book/10.1007/978-3-319-29501-5)

[3] [Sun X.A., Conejo A.J. (2021) Static Robust Optimization. In: Robust Optimization in Electric Energy Systems. International Series in Operations Research & Management Science, vol 313. Springer, Cham.]( https://doi.org/10.1007/978-3-030-85128-6_2)
