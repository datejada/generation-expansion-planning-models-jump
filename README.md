# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES) using Julia/JuMP

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/datejada/generation-expansion-planning-models-jump/main)

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP.jl**: Two-Stage Stochastic Generation Expansion Planning

These examples show basic concepts for learning optimization under uncertainty in power systems.

The models are developed in [Julia](https://julialang.org/), using the package [JuMP](https://jump.dev/JuMP.jl/stable/), and solved using [HiGHS](https://highs.dev/).

The primary reference to model the optimization problems is:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)

## GEP Formulation

### Indices

| **Name** | **Description**                   |
|----------|-----------------------------------|
| $p$      | time periods                      |
| $g$      | generation technologies           |
| $r(g)$   | subset of renewable techonologies |
| $sc$     | scenarios                         |

### Parameters

| **Name**   | **Domains** | **Description**                                             |
|------------|-------------|-------------------------------------------------------------|
| $pVOLL   $ |             | Value of Lost Load [\$/MWh]                                 |
| $pWeight $ |             | Representative period weight [hours]                        |
| $pInvCost$ | $g$         | Investment cost [\$/MW]                                     |
| $pVarCost$ | $g$         | Variable production cost [\$/MWh]                           |
| $pUnitCap$ | $g$         | Capacity per each invested unit [MW/unit]                   |
| $pRenProf$ | $r,p,sc$    | Renewable profile (e.g., load factor) [p.u.]                |
| $pDemand $ | $p$         | Demand [MW]                                                 |
| $pScProb $ | $sc$        | Scenario probability [p.u.]                                 |

### Variables

| **Name**    | **Domains** | **Description**              |
|-------------|-------------|------------------------------|
| $vTotCost $ |             | Total system cost [\$]       |
| $vInvCost $ |             | Total investment cost [\$]   |
| $vOpeCost $ |             | Total operating cost [\$]    |
| $vGenInv  $ | $g$         | Generation investment [1..N] |
| $vGenProd $ | $g,p,sc$    | Generation production [MW]   |
| $vLossLoad$ | $p,sc$      | Loss of load [MW]            |

### Equations

| **Name**                                    | **Domains** | **Description**                    |
|---------------------------------------------|-------------|------------------------------------|
| [eObjFun](#eobjfun)                         |             | Total system cost      [\$]        |
| [eInvCost](#einvcost)                       |             | Total investment cost      [\$]    |
| [eOpeCost](#eopecost)                       |             | Total operating cost      [\$]     |
| [eBalance](#ebalance)                       | $p,sc$      | Power system balance   [MWh]       |
| [eMaxProd](#emaxprod)                       | $g,p,sc$    | Maximum generation production [MW] |
| [eRenProd](#erenprod)                       | $r,p,sc$    | Maximum renewable production [MW]  |

#### *eObjFun*

$$
\displaystyle{\min{vTotCost = vInvCost + vOpeCost}}
$$

#### *eInvCost*

$$
vInvCost = \displaystyle \sum_{g}(pInvCost_{g} \cdot pUnitCap_{g} \cdot vGenInv_{g})
$$

#### *eOpeCost*

$$
vOpeCost = pWeight \cdot {\left(\displaystyle \sum_{sc}pScProb_{sc}\cdot{\left(\sum_{g,p}pVarCost_{g} \cdot vGenProd_{g,p,sc} + \sum_{p,sc}pVOLL \cdot vLossLoad_{p,sc}\right)}\right)}
$$

#### *eBalance*

$$
\displaystyle \sum_{g}vGenProd_{g,p,sc} + vLossLoad_{p,sc} = pDemand_{p} \quad \forall{p,sc}
$$

#### *eMaxProd*

$$
vGenProd_{g,p,sc} \leq pUnitCap_{g} \cdot vGenInv_{g} \quad \forall{g,p,sc}
$$

#### *eRenProd*

$$
vGenProd_{r,p,sc} \leq pRenProf_{r,p,sc} \cdot pUnitCap_{r} \cdot vGenInv_{r} \quad \forall{r,p,sc}
$$

#### *Bounds*

$vGenProd_{g,p,sc}\geq 0 ~ \forall g, p, sc $

$vLossLoad_{p,sc}\geq 0 ~ \forall p, sc $

$vGenInv_{g} \in \mathbb{Z}^{+} ~ \forall g $
