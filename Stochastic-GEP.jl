# Two-Stage Stochastic Generation Expansion Planning

# Packages
using JuMP
using HiGHS
using CSV
using DataFrames
using Plots
using StatsPlots

# Folders names
const INPUT_FOLDER = joinpath(@__DIR__, "inputs")
const OUTPUT_FOLDER = joinpath(@__DIR__, "outputs")

# Files names
demand_file = joinpath(INPUT_FOLDER, "iGEP_Data_Demand.csv")
generation_file = joinpath(INPUT_FOLDER, "iGEP_Data_Generation.csv")
availability_file = joinpath(INPUT_FOLDER, "iGEP_Data_Availability.csv")
scenario_file = joinpath(INPUT_FOLDER, "iGEP_Data_Scenario.csv")

# Read data
demand_df = CSV.read(demand_file, DataFrames.DataFrame)
generation_df = CSV.read(generation_file, DataFrames.DataFrame)
availability_df = CSV.read(availability_file, DataFrames.DataFrame)
scenario_df = CSV.read(scenario_file, DataFrames.DataFrame)

# Scalar parameters
p_rp_weight = 365  #weight of representative period [days]
p_ens_cost = 0.180 #energy not supplied cost    [kEUR/MWh]

# Sets
P = demand_df.p     #time periods (e.g., hours)
SC = scenario_df.sc #scenarios
G = generation_df.g #generation units

# Parameters
## availability profile [p.u.]
p_availability = Dict(
    (row.sc, row.g, row.p) => row.pAviProf
    for row in eachrow(availability_df)
)
## demand per time period [MW]
p_demand = Dict(
    (row.p) => row.pDemand
    for row in eachrow(demand_df)
)
## investment cost of generation units [kEUR/MW/year]
p_investment_cost = Dict(
    (row.g) => row.pInvCost
    for row in eachrow(generation_df)
)
## variable cost of generation units [kEUR/MWh]
p_variable_cost = Dict(
    (row.g) => row.pVarCost
    for row in eachrow(generation_df)
)
## capacity of generation units [MW]
p_unit_capacity = Dict(
    (row.g) => row.pUnitCap
    for row in eachrow(generation_df)
)
## probability of scenario [p.u.]
p_sc_prob = Dict(
    (row.sc) => row.pScProb
    for row in eachrow(scenario_df)
)

# Model
model = Model(optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0))

# Variables
@variable(model, 0 ≤ v_production[SC, G, P])          #production [MW] 
@variable(model, 0 ≤ v_ens[SC, p in P] ≤ p_demand[p]) #energy not supplied [MW]
@variable(model, 0 ≤ v_investment[G], Int)            #number of installed generation units [N]

# Expressions
e_investment_cost = @expression(model,
    sum(p_investment_cost[g] * p_unit_capacity[g] * v_investment[g]
        for g in G)
)
e_variable_cost = @expression(model,
    p_rp_weight * sum(p_sc_prob[sc] * p_variable_cost[g] * v_production[sc, g, p]
                      for sc in SC, g in G, p in P)
)
e_ens_cost = @expression(model,
    p_rp_weight * sum(p_sc_prob[sc] * p_ens_cost * v_ens[sc, p]
                      for sc in SC, p in P)
)

# Objective function
@objective(model, Min,
    e_investment_cost + e_variable_cost + e_ens_cost
)

# Constraints
# - balance equation
@constraint(model, c_balance[sc in SC, p in P],
    sum(v_production[sc, g, p] for g in G) + v_ens[sc, p] == p_demand[p]
)

# - maximum generation
@constraint(model, c_max_prod[sc in SC, g in G, p in P],
    v_production[sc, g, p] <= get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] * v_investment[g]
)

# print lp file
write_to_file(model, "model.lp")

# Solve model
optimize!(model)

# Objective function value
@show objective_value(model)

# Writing the investment results to a CSV file
output_file = open(joinpath(OUTPUT_FOLDER, "oGEP_Invest_Result.csv"), "w")
write(output_file, "g,InstalUnits,InstalCap_MW\n")
for g in G
    write(output_file, "$g,$(value(v_investment[g])),$(value(p_unit_capacity[g]) * value(v_investment[g]))\n")
end
close(output_file)

# Plot the investment results
#plotly() # uncomment to use plotly backend
df = CSV.read(joinpath(OUTPUT_FOLDER, "oGEP_Invest_Result.csv"), DataFrames.DataFrame)
@df df bar(:g,
    :InstalCap_MW,
    title="Investment Results",
    xlabel="Generation Units",
    ylabel="Installed Capacity [MW]",
    legend=false,
    dpi=300,
)
#savefig(joinpath(OUTPUT_FOLDER, "investment_results.png"))
