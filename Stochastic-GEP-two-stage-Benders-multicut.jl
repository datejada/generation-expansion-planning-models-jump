# Two-Stage Stochastic Generation Expansion Planning

# Packages
using JuMP
using HiGHS
using CSV
using DataFrames
using Plots
using StatsPlots
using Printf

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
    availability    = Dict((row.sc, row.g, row.p) => row.pAviProf for row in eachrow(availability_df)) #availability profile [p.u.]
    demand          = Dict((row.p) => row.pDemand for row in eachrow(demand_df))        #demand per time period [MW]
    investment_cost = Dict((row.g) => row.pInvCost for row in eachrow(generation_df))   #investment cost of generation units [kEUR/MW/year]
    variable_cost   = Dict((row.g) => row.pVarCost for row in eachrow(generation_df))   #variable   cost of generation units [kEUR/MWh]
    unit_capacity   = Dict((row.g) => row.pUnitCap for row in eachrow(generation_df))   #capacity        of generation units [MW]
    sc_prob         = Dict((row.sc) => row.pScProb for row in eachrow(scenario_df))     #probability of scenario [p.u.]
    rp_weight       = 365   #weight of representative period [days]
    ens_cost        = 0.180 #energy not supplied cost    [kEUR/MWh]

    params = Dict(
        :availability    => availability,
        :demand          => demand,
        :investment_cost => investment_cost,
        :variable_cost   => variable_cost,
        :unit_capacity   => unit_capacity,
        :sc_prob         => sc_prob,
        :rp_weight       => rp_weight,
        :ens_cost        => ens_cost,
    )
    return sets, params
end

"""
    create_first_stage_model(sets, params)

This function creates the first-stage subproblem in the Benders' Decomposition.

# Arguments
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.

# Returns
- `model::Model`: JuMP model with the first-stage problem.
"""
function create_first_stage_model(sets, params)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]

    # Extract parameters
    p_investment_cost = params[:investment_cost]
    p_unit_capacity   = params[:unit_capacity]
    p_sc_prob         = params[:sc_prob]

    # Scalar values
    M = -1000 # to avoid the problem is unbounded 

    # Model
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0, "output_flag" => false),
    )

    # Variables
    @variable(model, 0 ≤ v_investment[G], Int)  #number of installed generation units [N]
    @variable(model, v_theta[SC] ≥ M)           #Benders' cut per scenario [kEUR]

    # Expressions
    e_investment_cost = @expression(
        model,
        sum(p_investment_cost[g] * p_unit_capacity[g] * v_investment[g] for g in G)
    )

    @expression(model, e_benders_cut, sum(p_sc_prob[sc] * v_theta[sc] for sc in SC))

    # Objective function
    @objective(model, Min, e_investment_cost + e_benders_cut)

    return model
end

"""
    create_and_solve_subproblem(sets, params)

This function creates the first-stage subproblem in the Benders' Decomposition.

# Arguments
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.
- `p_investment::Array{Float64}`: The optimal investment solution from the first-stage problem.

# Returns
- `model::Model`: JuMP model with the subproblem.
"""
function create_and_solve_subproblem(sets, params, p_investment)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Extract parameters
    p_availability  = params[:availability]
    p_demand        = params[:demand]
    p_variable_cost = params[:variable_cost]
    p_unit_capacity = params[:unit_capacity]
    p_sc_prob       = params[:sc_prob]
    p_rp_weight     = params[:rp_weight]
    p_ens_cost      = params[:ens_cost]

    # Model
    model = Model(optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => false))

    # Variables
    @variable(model, 0 ≤ v_production[SC, G, P])          #production [MW] 
    @variable(model, 0 ≤ v_ens[SC, p in P] ≤ p_demand[p]) #energy not supplied [MW]

    # Expressions
    @expression(
        model,
        e_variable_cost[sc in SC],
        p_rp_weight * sum(p_variable_cost[g] * v_production[sc, g, p] for g in G, p in P)
    )

    @expression(
        model,
        e_ens_cost[sc in SC],
        p_rp_weight * sum(p_ens_cost * v_ens[sc, p] for p in P)
    )

    @expression(model, e_cost_per_scenario[sc in SC], e_variable_cost[sc] + e_ens_cost[sc])

    @expression(model, e_expected_cost, sum(p_sc_prob[sc] * e_cost_per_scenario[sc] for sc in SC))

    # Objective function
    @objective(model, Min, sum(e_cost_per_scenario[sc] for sc in SC))

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
        get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] * p_investment[g]
    )

    # Solve model
    optimize!(model)

    # Check if the model is optimal
    @assert is_solved_and_feasible(model; dual = true)

    # Return results
    return model
end

"""
    add_cut(sets, params)

This function creates the cut to add in Benders' Decomposition algorithm.

# Arguments
- `first_stage_model::Model`: The first-stage model.
- `subproblem::NamedTuple`: The subproblem results.
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.
- `iteration::Int`: The current iteration number.

# Returns
- `nothing`
"""
function add_cut(first_stage_model, subproblem, sets, params, p_investment, iteration)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Extract parameters
    p_availability  = params[:availability]
    p_unit_capacity = params[:unit_capacity]
    p_sc_prob       = params[:sc_prob]

    # Get parameters from the subproblem
    p_subproblem_obj_per_sc = value.(subproblem[:e_cost_per_scenario])
    p_dual = dual.(subproblem[:c_max_prod])

    # Add Benders' cut
    for sc in SC
        @constraint(
            first_stage_model,
            base_name = "cut_iter_$(iteration)_$(sc)",
            first_stage_model[:v_theta][sc] >=
            p_subproblem_obj_per_sc[sc] + sum(
                -p_dual[sc, g, p] *
                get(p_availability, (sc, g, p), 1.0) *
                p_unit_capacity[g] *
                (p_investment[g] - first_stage_model[:v_investment][g]) for g in G, p in P
            )
        )
    end

    return nothing
end

"""
    print_iteration(iteration)

Prints the current iteration number.

# Arguments
- `iteration`: An integer representing the current iteration number.
"""
function print_iteration(k, args...)
    f(x) = Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

"""
    plot_benders_gap(df_Benders_interations)

Plot the Benders' gap.

# Arguments
- `df_Benders_interations`: A DataFrame containing the Benders' iterations.

# Returns
- `p`: A plot object showing the Benders' gap.

"""
function plot_benders_gap(df_Benders_interations)
    p = @df df_Benders_interations plot(
        :iteration,
        :gap,
        label = "Gap",
        xlabel = "Iteration",
        ylabel = "Gap",
        title = "Benders Decomposition",
        lw = 2,
        xticks = 0:1:MAXIMUM_ITERATIONS,
        color = :darkblue,
    )
    return p
end

"""
    plot_benders_bounds(df_Benders_interations)

Plot the Benders' bounds.

# Arguments
- `df_Benders_interations`: A DataFrame containing the Benders' bounds.

# Returns
- `p`: A plot object showing the Benders' bounds.

"""
function plot_benders_bounds(df_Benders_interations)
    p = @df df_Benders_interations plot(
        :iteration,
        [:lower_bound, :upper_bound],
        label = ["Lower Bound" "Upper Bound"],
        xlabel = "Iteration",
        ylabel = "kEUR",
        title = "Benders Decomposition",
        lw = 2,
        xticks = 0:1:MAXIMUM_ITERATIONS,
    )
    return p
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
    v_investment_cap = [v_investment[g] * params[:unit_capacity][g] for g in sets[:G]]

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
        demand = [params[:demand][p] for p in periods]
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
        plot!(p, periods, 1000 * dual / params[:rp_weight]; label = sc, linewidth = 2)
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
function save_results_first_stage_to_csv(output_folder, model)
    CSV.write(
        joinpath(output_folder, "oGEP_Investment.csv"),
        Containers.rowtable(value, model[:v_investment]; header = [:g, :instal_units]),
    )

    return nothing
end

"""
    save_results_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""
function save_results_subproblem_to_csv(output_folder, model)
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

