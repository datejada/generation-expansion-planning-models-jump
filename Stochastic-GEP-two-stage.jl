# Two-Stage Stochastic Generation Expansion Planning

# Packages
using JuMP
using HiGHS
using CSV
using DataFrames
using Plots
using StatsPlots

# Functions

"""
    read_data(input_folder)

Reads data from CSV files and returns sets and parameters.

# Arguments
- `input_folder::String`: The path to the folder containing the input CSV files.

# Returns
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets P, SC, and G.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters `p_availability`, `p_demand`, `p_investment_cost`, `p_variable_cost`, `p_unit_capacity`, `p_sc_prob`, `p_rp_weight`, and `p_ens_cost`.
"""
function read_data(input_folder)
    # Files names
    demand_file       = joinpath(input_folder, "iGEP_Data_Demand.csv")
    generation_file   = joinpath(input_folder, "iGEP_Data_Generation.csv")
    availability_file = joinpath(input_folder, "iGEP_Data_Availability.csv")
    scenario_file     = joinpath(input_folder, "iGEP_Data_Scenario.csv")

    # Read data
    demand_df       = CSV.read(demand_file, DataFrames.DataFrame)
    generation_df   = CSV.read(generation_file, DataFrames.DataFrame)
    availability_df = CSV.read(availability_file, DataFrames.DataFrame)
    scenario_df     = CSV.read(scenario_file, DataFrames.DataFrame)

    # Sets
    P  = demand_df.p     #time periods (e.g., hours)
    SC = scenario_df.sc  #scenarios
    G  = generation_df.g #generation units

    sets = Dict(:P => P, :SC => SC, :G => G)

    # Parameters
    p_availability    = Dict((row.sc, row.g, row.p) => row.pAviProf for row in eachrow(availability_df)) #availability profile [p.u.]
    p_demand          = Dict((row.p) => row.pDemand for row in eachrow(demand_df))        #demand per time period [MW]
    p_investment_cost = Dict((row.g) => row.pInvCost for row in eachrow(generation_df))   #investment cost of generation units [kEUR/MW/year]
    p_variable_cost   = Dict((row.g) => row.pVarCost for row in eachrow(generation_df))   #variable   cost of generation units [kEUR/MWh]
    p_unit_capacity   = Dict((row.g) => row.pUnitCap for row in eachrow(generation_df))   #capacity        of generation units [MW]
    p_sc_prob         = Dict((row.sc) => row.pScProb for row in eachrow(scenario_df))     #probability of scenario [p.u.]
    p_rp_weight       = 365   #weight of representative period [days]
    p_ens_cost        = 0.180 #energy not supplied cost    [kEUR/MWh]

    params = Dict(
        :p_availability    => p_availability,
        :p_demand          => p_demand,
        :p_investment_cost => p_investment_cost,
        :p_variable_cost   => p_variable_cost,
        :p_unit_capacity   => p_unit_capacity,
        :p_sc_prob         => p_sc_prob,
        :p_rp_weight       => p_rp_weight,
        :p_ens_cost        => p_ens_cost,
    )
    return sets, params
end

"""
    create_and_solve_model(sets, params)

This function creates and solves a mathematical optimization model for a central planner problem.

# Arguments
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.

# Returns
- `model::Model`: The optimized model if a solution is found.
- `String`: A warning message if no solution is found.
"""
function create_and_solve_model(sets, params)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Extract parameters
    p_availability    = params[:p_availability]
    p_demand          = params[:p_demand]
    p_investment_cost = params[:p_investment_cost]
    p_variable_cost   = params[:p_variable_cost]
    p_unit_capacity   = params[:p_unit_capacity]
    p_sc_prob         = params[:p_sc_prob]
    p_rp_weight       = params[:p_rp_weight]
    p_ens_cost        = params[:p_ens_cost]

    # Model
    model = Model(optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0))

    # Variables
    @variable(model, 0 ≤ v_production[SC, G, P])          #production [MW] 
    @variable(model, 0 ≤ v_ens[SC, p in P] ≤ p_demand[p]) #energy not supplied [MW]
    @variable(model, 0 ≤ v_investment[G], Int)            #number of installed generation units [N]

    # Expressions
    e_investment_cost = @expression(model, sum(p_investment_cost[g] * p_unit_capacity[g] * v_investment[g] for g in G))
    e_variable_cost   = @expression(model, p_rp_weight * sum(p_sc_prob[sc] * p_variable_cost[g] * v_production[sc, g, p] for sc in SC, g in G, p in P))
    e_ens_cost        = @expression(model, p_rp_weight * sum(p_sc_prob[sc] * p_ens_cost * v_ens[sc, p] for sc in SC, p in P))

    # Objective function
    @objective(model, Min, e_investment_cost + e_variable_cost + e_ens_cost)

    # Constraints
    # - balance equation
    @constraint(
        model,
        c_balance[sc in SC, p in P],
        sum(v_production[sc, g, p] for g in G) + v_ens[sc, p] == p_demand[p]
    )

    # - maximum generation
    @constraint(
        model,
        c_max_prod[sc in SC, g in G, p in P],
        v_production[sc, g, p] <=
        get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] * v_investment[g]
    )

    # Solve model
    optimize!(model)

    # Check if the model is optimal
    if termination_status(model) == MOI.OPTIMAL
        @show objective_value(model)
        return model
    else
        return "Warning: No solution was found."
    end
end

"""
    plot_investment(model)

Plot the investment in generation units.

# Arguments
- `model`: The model containing the investment variable.

"""
function plot_investment(model, sets, params)
    # Extract the variable v_investment from the model
    v_investment = value.(model[:v_investment])

    # Calculate the investment capacity in MW
    v_investment_cap = [v_investment[g] * params[:p_unit_capacity][g] for g in sets[:G]]

    # Extract generator names
    generator_names = [string(k[1]) for k in keys(v_investment)]

    # Create a bar chart subplots
    p = plot(; layout = (2, 1))

    # Create a bar chart for the investment in MW
    bar!(
        generator_names,
        v_investment_cap;
        xlabel = "",
        ylabel = "Capacity [MW]",
        title = "Investment Results",
        legend = false,
        subplot = 1,
    )

    # Create a bar chart for the number of installed units
    bar!(
        generator_names,
        Array(v_investment); # Convert the DenseArray from JuMP to an array
        xlabel = "Generation Technology",
        ylabel = "Units [N]",
        title = "",
        legend = false,
        subplot = 2,
    )

    return p
end

"""
    get_production_tmp(scenario::String, model)

This function takes a scenario name and a model as input and returns a temporary production table for that scenario.

# Arguments
- `scenario::String`: The name of the scenario.
- `model`: The model containing the production and energy storage data.

# Returns
- `production_tmp`: The temporary production table for the given scenario.
"""
function get_production_tmp(scenario, model)
    production_table =
        Containers.rowtable(value, model[:v_production]; header = [:sc, :g, :p, :value])
    production_df = DataFrames.DataFrame(production_table)

    ens_table = Containers.rowtable(value, model[:v_ens]; header = [:sc, :p, :value])
    ens_df = DataFrames.DataFrame(ens_table)

    filter!(row -> row.sc == scenario, production_df)
    filter!(row -> row.sc == scenario, ens_df)

    # drop the scenario column 
    production_tmp = production_df[:, 2:end]

    # unstack column g
    production_tmp = unstack(production_tmp, :g, :value)

    # drop column P
    production_tmp = production_tmp[:, 2:end]

    # add column ens to production_tmp
    production_tmp = hcat(production_tmp, ens_df[:, 3])

    return production_tmp
end

"""
    plot_production(model, params)

Plot the production of generators and the demand for each scenario.

# Arguments
- `model`: The optimization model containing the production variables.
- `params`: A dictionary containing the parameters of the model.

# Returns
- `p`: A plot object showing the production and demand for each scenario.

"""
function plot_production(model, params)
    # Extract the variable v_production and v_ens from the model
    v_production = value.(model[:v_production])

    # Extract unique scenarios, generators, and periods
    scenarios  = unique([k[1] for k in keys(v_production)])
    generators = unique([k[2] for k in keys(v_production)])
    generators = push!(generators, "ens")
    periods    = unique([k[3] for k in keys(v_production)])

    # Create a plot for each scenario
    p = plot(; layout = (length(scenarios), 1)) # Create subplots for each scenario
    for (i, sc) in enumerate(scenarios)
        # 
        total_production_df = get_production_tmp(sc, model)
        groupedbar!(
            Matrix(total_production_df);
            bar_position = :stack,
            labels = reshape(generators, 1, length(generators)),
            legend = :topleft,
            title = sc,
            subplot = i,
        )

        # Add demand as a black line
        demand = [params[:p_demand][p] for p in periods]
        plot!(p[i], demand; label = "demand", color = :black, linewidth = 2, linestyle = :dash)
    end

    # Set plot attributes
    plot!(; ylabel = "Production", legend = :outertopright)

    return p
end

"""
    plot_dual_balance(model, params)

Plot the dual variable of the balance equation for different scenarios and periods.

# Arguments
- `model`: The optimization model.
- `params`: A dictionary of parameters.

# Returns
- `p`: The plot object.

"""
function plot_dual_balance(model, params)
    # Extract the dual variable of the balance equation
    dual_balance = dual.(model[:c_balance])

    # Extract unique scenarios and periods
    scenarios = unique([k[1] for k in keys(dual_balance)])
    periods = unique([k[2] for k in keys(dual_balance)])

    # Create a plot
    p = plot() # Create a single plot

    # For each scenario, create a series
    for sc in scenarios
        dual = [dual_balance[sc, p] for p in periods]
        plot!(p, periods, 1000 * dual / params[:p_rp_weight]; label = sc, linewidth = 2)
    end

    # Set plot attributes
    plot!(; ylabel = "EUR/MWh", legend = :topleft)

    return p
end

"""
    save_results_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""

function save_results_to_csv(output_folder, model)
    CSV.write(
        joinpath(output_folder, "oGEP_Investment.csv"),
        Containers.rowtable(value, model[:v_investment]; header = [:g, :instal_units]),
    )
    CSV.write(
        joinpath(output_folder, "oGEP_Production.csv"),
        Containers.rowtable(value, model[:v_production]; header = [:sc, :g, :h, :production]),
    )

    CSV.write(
        joinpath(output_folder, "oGEP_ENS.csv"),
        Containers.rowtable(value, model[:v_ens]; header = [:sc, :h, :ens]),
    )

    return nothing
end
