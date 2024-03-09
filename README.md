# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using Julia/JuMP

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP.jl**: Two-Stage Stochastic Generation Expansion Planning

The models are developed in [Julia](https://julialang.org/), using the package [JuMP](https://jump.dev/JuMP.jl/stable/), and solved using [HiGHS](https://highs.dev/).

The primary reference to model the optimization problems is:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)
