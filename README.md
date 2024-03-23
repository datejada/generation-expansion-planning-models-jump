# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using Julia/JuMP

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datejada/generation-expansion-planning-models-jump/main)

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP-notebook.ipynb**: Two-Stage Stochastic Generation Expansion Planning

These examples show basic concepts for learning optimization under uncertainty in power systems.

The models are developed in [Julia](https://julialang.org/), using the package [JuMP](https://jump.dev/JuMP.jl/stable/), and solved using [HiGHS](https://highs.dev/).

The primary reference to model the optimization problems is:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)
