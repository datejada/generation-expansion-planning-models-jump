"""
    read_two_stage_data(input_folder)

Reads data from CSV files and returns sets and parameters.

# Arguments
- `input_folder::String`: The path to the folder containing the input CSV files.

# Returns
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets P, SC, and G.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters `p_availability`, `p_demand`, `p_investment_cost`, `p_variable_cost`, `p_unit_capacity`, `p_sc_prob`, `p_rp_weight`, and `p_ens_cost`.
"""
function read_two_stage_data(input_folder)
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
    read_multi_stage_data(input_folder)

Reads data from CSV files and returns sets and parameters.

# Arguments
- `input_folder::String`: The path to the folder containing the input CSV files.

# Returns
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets P, ST, SC, and G.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters `availability`, `demand`, `investment_cost`, `variable_cost`, `unit_capacity`, `tree_node_prob`, `stochastic_paths`, `rp_weight`, `ens_cost`, and `discount_rate`.
"""
function read_multi_stage_data(input_folder)
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
    read_aro_data(input_folder)

Reads data from CSV files and returns sets and parameters.

# Arguments
- `input_folder::String`: The path to the folder containing the input CSV files.

# Returns
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets P, and G.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters `p_availability`, `p_demand`, `p_investment_cost`, `p_variable_cost`, `p_unit_capacity`, `p_rp_weight`, and `p_ens_cost`.
"""
function read_aro_data(input_folder)
    # Files names
    demand_file       = joinpath(input_folder, "iGEP_Data_Demand.csv")
    generation_file   = joinpath(input_folder, "iGEP_Data_Generation.csv")
    availability_file = joinpath(input_folder, "iGEP_Data_Availability.csv")

    # Read data
    demand_df       = CSV.read(demand_file, DataFrames.DataFrame)
    generation_df   = CSV.read(generation_file, DataFrames.DataFrame)
    availability_df = CSV.read(availability_file, DataFrames.DataFrame)

    # Sets
    P = demand_df.p     #time periods (e.g., hours)
    G = generation_df.g #generation units

    sets = Dict(:P => P, :G => G)

    # Parameters
    max_availability    = Dict((row.g, row.p) => row.pMaxAviProf for row in eachrow(availability_df)) #maximum availability profile [p.u.]
    min_availability    = Dict((row.g, row.p) => row.pMinAviProf for row in eachrow(availability_df)) #minimum availability profile [p.u.]
    demand              = Dict((row.p) => row.pDemand for row in eachrow(demand_df))         #demand per time period [MW]
    investment_cost     = Dict((row.g) => row.pInvCost for row in eachrow(generation_df))    #investment cost of generation units [kEUR/MW/year]
    variable_cost       = Dict((row.g) => row.pVarCost for row in eachrow(generation_df))    #variable   cost of generation units [kEUR/MWh]
    unit_capacity       = Dict((row.g) => row.pUnitCap for row in eachrow(generation_df))    #capacity        of generation units [MW]
    availability_factor = Dict((row.g) => row.pAvaiFactor for row in eachrow(generation_df)) #generation availability factor [p.u.]
    rp_weight           = 365   #weight of representative period [days]
    ens_cost            = 0.180 #energy not supplied cost   [kEUR/MWh]
    exc_cost            = 0.180 #excess cost                [kEUR/MWh]
    uncertainty_budget  = 0.2   #uncertainty budget

    params = Dict(
        :max_availability    => max_availability,
        :min_availability    => min_availability,
        :demand              => demand,
        :investment_cost     => investment_cost,
        :variable_cost       => variable_cost,
        :unit_capacity       => unit_capacity,
        :availability_factor => availability_factor,
        :rp_weight           => rp_weight,
        :ens_cost            => ens_cost,
        :exc_cost            => exc_cost,
        :uncertainty_budget  => uncertainty_budget,
    )
    return sets, params
end