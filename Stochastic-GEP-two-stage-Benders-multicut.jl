# Two-Stage Stochastic Generation Expansion Planning

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
