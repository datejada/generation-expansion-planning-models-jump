# Multi-Stage Stochastic Generation Expansion Planning

# Packages
using JuMP
using HiGHS
using CSV
using DataFrames
using Plots
using StatsPlots

# Functions

"""
Function to find the stochastic path from the root to a given node in the scenario tree.

Parameters:
- scenario_tree (DataFrame): The scenario tree data.
- stage (String): The stage of the target node.
- scenario (String): The scenario identifier of the target node.

Returns:
- path (Vector{Tuple{String, String}}): A list of tuples, where each tuple represents a node in the path,
  formatted as (stage, scenario).
"""
function get_stochastic_path_to_node(scenario_tree, stage, scenario)
    # Initialize the path list
    path = []

    # Find the target node
    target_node = filter(row -> row[:st] == stage && row[:sc] == scenario, scenario_tree)

    # Check if the target node exists
    if isempty(target_node)
        return push!(path, ("empty", "empty"))
    end

    # Loop to trace back to the root from the target node
    while target_node[1, :sc_ancestor] != "none"
        # Extract details of the current node
        node_info = target_node[1, :]
        pushfirst!(path, (node_info[:st], node_info[:sc]))

        # Move to the ancestor node
        ancestor_scenario = node_info[:sc_ancestor]
        target_node = filter(row -> row[:sc] == ancestor_scenario, scenario_tree)
    end

    # Add the root node to the path
    root_node = filter(row -> row[:sc] == "root", scenario_tree)[1, :]
    pushfirst!(path, (root_node[:st], root_node[:sc]))

    return path
end

"""
    read_data(input_folder)

Reads data from CSV files and returns sets and parameters.

# Arguments
- `input_folder::String`: The path to the folder containing the input CSV files.

# Returns
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets P, ST, SC, and G.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters `availability`, `demand`, `investment_cost`, `variable_cost`, `unit_capacity`, `tree_node_prob`, `stochastic_paths`, `rp_weight`, `ens_cost`, and `discount_rate`.
"""
function read_data(input_folder)
    # Files names
    demand_file        = joinpath(input_folder, "iGEP_Data_Demand.csv")
    generation_file    = joinpath(input_folder, "iGEP_Data_Generation.csv")
    availability_file  = joinpath(input_folder, "iGEP_Data_Availability.csv")
    scenario_tree_file = joinpath(input_folder, "iGEP_Data_ScenarioTree.csv")

    # Read data
    demand_df        = CSV.read(demand_file, DataFrames.DataFrame)
    generation_df    = CSV.read(generation_file, DataFrames.DataFrame)
    availability_df  = CSV.read(availability_file, DataFrames.DataFrame)
    scenario_tree_df = CSV.read(scenario_tree_file, DataFrames.DataFrame)

    # Sets
    P  = unique(demand_df.p)         #time periods (e.g., hours)
    ST = unique(scenario_tree_df.st) #stages
    SC = unique(scenario_tree_df.sc) #scenarios
    G  = generation_df.g             #generation units

    sets = Dict(:P => P, :ST => ST, :SC => SC, :G => G)

    # Parameters
    availability    = Dict((row.st, row.sc, row.g, row.p) => row.pAviProf for row in eachrow(availability_df)) #availability profile [p.u.]
    demand          = Dict((row.st, row.p) => row.pDemand for row in eachrow(demand_df))                       #demand per time period [MW]
    investment_cost = Dict((row.g) => row.pInvCost for row in eachrow(generation_df))                          #investment cost of generation units [kEUR/MW/year]
    variable_cost   = Dict((row.g) => row.pVarCost for row in eachrow(generation_df))                          #variable   cost of generation units [kEUR/MWh]
    unit_capacity   = Dict((row.g) => row.pUnitCap for row in eachrow(generation_df))                          #capacity        of generation units [MW]
    tree_node_prob  = Dict((row.st, row.sc) => row.pNodeProb for row in eachrow(scenario_tree_df))             #probability of scenario [p.u.]
    st_order        = Dict(value => index for (index, value) in enumerate(ST))                                 #stage order
    rp_weight       = 365   #weight of representative period [days]
    ens_cost        = 0.18  #energy not supplied cost    [kEUR/MWh]
    discount_rate   = 0.05  #discount rate

    stochastic_paths = Dict(
        (st, sc) => get_stochastic_path_to_node(scenario_tree_df, st, sc) for st in ST for sc in SC
    )

    params = Dict(
        :availability     => availability,
        :demand           => demand,
        :investment_cost  => investment_cost,
        :variable_cost    => variable_cost,
        :unit_capacity    => unit_capacity,
        :tree_node_prob   => tree_node_prob,
        :stochastic_paths => stochastic_paths,
        :st_order         => st_order,
        :rp_weight        => rp_weight,
        :ens_cost         => ens_cost,
        :discount_rate    => discount_rate,
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
    ST = sets[:ST]
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Extract parameters
    p_availability     = params[:availability]
    p_demand           = params[:demand]
    p_investment_cost  = params[:investment_cost]
    p_variable_cost    = params[:variable_cost]
    p_unit_capacity    = params[:unit_capacity]
    p_tree_node_prob   = params[:tree_node_prob]
    p_stochastic_paths = params[:stochastic_paths]
    p_st_order         = params[:st_order]
    p_rp_weight        = params[:rp_weight]
    p_ens_cost         = params[:ens_cost]
    p_discount_rate    = params[:discount_rate]

    # Subsets that depens on parameters values
    tree_node_operational = Dict(
        (st, sc) => get(p_tree_node_prob, (st, sc), 0) > 0 && sc != "root" for st in ST, sc in SC
    )

    tree_node_investment = Dict(
        (st, sc) => get(p_tree_node_prob, (st, sc), 0) > 0 && st != ST[end] for st in ST, sc in SC
    )

    # Model
    model = Model(optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0))

    # Variables
    @variable(model, 0 ≤ v_production[st in ST, sc in SC, G, P; tree_node_operational[(st, sc)]]) #production [MW] 

    @variable(
        model,
        0 ≤ v_ens[st in ST, sc in SC, p in P; tree_node_operational[(st, sc)]] ≤ p_demand[st, p]
    ) #energy not supplied [MW]

    @variable(model, 0 ≤ v_investment[st in ST, sc in SC, G; tree_node_investment[(st, sc)]], Int)                  #number of installed generation units [N]

    # Expressions
    e_investment_cost = @expression(
        model,
        sum(
            (1 / (1 + p_discount_rate)^(p_st_order[st] - 1)) *
            p_tree_node_prob[st, sc] *
            p_investment_cost[g] *
            p_unit_capacity[g] *
            v_investment[st, sc, g] for
            st in ST, sc in SC, g in G if tree_node_investment[(st, sc)]
        )
    )

    e_variable_cost = @expression(
        model,
        p_rp_weight * sum(
            (1 / (1 + p_discount_rate)^(p_st_order[st] - 1)) *
            p_tree_node_prob[st, sc] *
            p_variable_cost[g] *
            v_production[st, sc, g, p] for
            st in ST, sc in SC, g in G, p in P if tree_node_operational[(st, sc)]
        )
    )

    e_ens_cost = @expression(
        model,
        p_rp_weight * sum(
            (1 / (1 + p_discount_rate)^(p_st_order[st] - 1)) *
            p_tree_node_prob[st, sc] *
            p_ens_cost *
            v_ens[st, sc, p] for st in ST, sc in SC, p in P if tree_node_operational[(st, sc)]
        )
    )

    @expression(
        model,
        e_accumulated_investment[st in ST, sc in SC, g in G],
        sum(
            v_investment[stst, scsc, g] for
            stst in ST, scsc in SC if tree_node_investment[(stst, scsc)] &&
            p_st_order[stst] < p_st_order[st] &&
            (stst, scsc) in p_stochastic_paths[(st, sc)]
        )
    )

    # Objective function
    @objective(model, Min, e_investment_cost + e_variable_cost + e_ens_cost)

    # Constraints
    # - balance equation
    @constraint(
        model,
        c_balance[st in ST, sc in SC, p in P; tree_node_operational[(st, sc)]],
        sum(v_production[st, sc, g, p] for g in G) + v_ens[st, sc, p] == p_demand[st, p]
    )

    # - maximum generation
    @constraint(
        model,
        c_max_prod[st in ST, sc in SC, g in G, p in P; tree_node_operational[(st, sc)]],
        v_production[st, sc, g, p] <=
        get(p_availability, (st, sc, g, p), 1.0) *
        p_unit_capacity[g] *
        e_accumulated_investment[st, sc, g]
    )

    # Write the model to a file
    #write_to_file(model, "model.lp") # Uncomment to write the model to a file

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
function plot_investment(model, params)
    investment_table =
        Containers.rowtable(value, model[:v_investment]; header = [:st, :sc, :g, :units])
    investment_df = DataFrames.DataFrame(investment_table)

    # Create a new column with the MW capacity
    investment_df[!, :capacity] =
        [row.units * params[:unit_capacity][row.g] for row in eachrow(investment_df)]

    # Create a bar chart subplots
    p = plot(; layout = (2, 1))

    # Create a bar chart for the investment in MW
    @df investment_df groupedbar!(
        :g,
        :capacity,
        group = (:st, :sc),
        xlabel = "",
        ylabel = "Capacity [MW]",
        title = "Investment Results",
        legend = true,
        #bar_position = :stack,
        subplot = 1,
    )

    # Create a bar chart for the number of installed units
    @df investment_df groupedbar!(
        :g,
        :units,
        group = (:st, :sc),
        xlabel = "Generation Technology",
        ylabel = "Units [N]",
        title = "",
        legend = false,
        #bar_position = :stack,
        subplot = 2,
    )

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
        Containers.rowtable(value, model[:v_investment]; header = [:st, :sc, :g, :instal_units]),
    )
    CSV.write(
        joinpath(output_folder, "oGEP_Production.csv"),
        Containers.rowtable(value, model[:v_production]; header = [:st, :sc, :g, :h, :production]),
    )

    CSV.write(
        joinpath(output_folder, "oGEP_ENS.csv"),
        Containers.rowtable(value, model[:v_ens]; header = [:st, :sc, :h, :ens]),
    )

    return nothing
end
