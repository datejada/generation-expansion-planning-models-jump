# Two-Stage Stochastic Generation Expansion Planning

"""
    create_first_stage_model(sets, params)

This function creates the first-stage problem in the Lagrangian Relaxation Decomposition.

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

    # Scalar values
    M = 1e9 # to avoid the problem is unbounded

    # Model
    model = Model(
        optimizer_with_attributes(
            HiGHS.Optimizer,
            "output_flag" => false,
            "mip_rel_gap" => 0.0,
            "dual_feasibility_tolerance" => 1e-9,
            "mip_feasibility_tolerance" => 1e-9,
            "primal_feasibility_tolerance" => 1e-9,
        ),
    )

    # Variables
    @variable(model, v_theta)                            # recourse function
    @variable(model, 0 ≤ v_lambda[SC, G] ≤ M, start = 0) # lagrangian multipliers

    # Objective function
    @objective(model, Max, v_theta)

    return model
end

"""
    create_and_solve_subproblem(sets, params)

This function creates the subproblem in the Lagrangian Relaxation Decomposition.

# Arguments
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.
- `p_lambda::Array{Float64}`: The optimal lagrangian multiplier solution from the first-stage problem.

# Returns
- `model::Model`: JuMP model with the subproblem.
"""
function create_and_solve_subproblem(sets, params, p_lambda)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Create set scenario lag
    SC_lag = circshift(SC, 1)

    # Extract parameters
    p_availability    = params[:availability]
    p_demand          = params[:demand]
    p_investment_cost = params[:investment_cost]
    p_variable_cost   = params[:variable_cost]
    p_unit_capacity   = params[:unit_capacity]
    p_sc_prob         = params[:sc_prob]
    p_rp_weight       = params[:rp_weight]
    p_ens_cost        = params[:ens_cost]

    # get the maximum units needed to cover the demand each hour per generator
    p_max_req_units_per_period = Dict(
        (sc, g, p) =>
            ceil(p_demand[p] / (get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] + 1))
        for sc in SC, g in G, p in P
    )

    # It is important to put bounds on the subproblem's variables
    # to avoid an unbounded subproblem. If the variables
    # do not have initial bounds, one can impose maximum bounds
    # that makes sense to the problem.
    p_max_investment = Dict(
        g => maximum([v for ((_, type, _), v) in p_max_req_units_per_period if type == g]) for
        g in G
    )

    # Model
    model = Model(
        optimizer_with_attributes(
            HiGHS.Optimizer,
            "output_flag" => false,
            "mip_rel_gap" => 0.0,
            "dual_feasibility_tolerance" => 1e-9,
            "mip_feasibility_tolerance" => 1e-9,
            "primal_feasibility_tolerance" => 1e-9,
        ),
    )

    # Variables
    @variable(model, 0 ≤ v_production[SC, G, P])                              #production [MW] 
    @variable(model, 0 ≤ v_ens[SC, p in P] ≤ p_demand[p])                     #energy not supplied [MW]
    @variable(model, 0 ≤ v_investment[SC, g in G] ≤ p_max_investment[g], Int) #number of installed generation units [N]

    # Expressions
    e_investment_cost = @expression(
        model,
        sum(
            p_sc_prob[sc] * p_investment_cost[g] * p_unit_capacity[g] * v_investment[sc, g] for
            sc in SC, g in G
        )
    )

    e_variable_cost = @expression(
        model,
        p_rp_weight * sum(
            p_sc_prob[sc] * p_variable_cost[g] * v_production[sc, g, p] for sc in SC, g in G,
            p in P
        )
    )

    e_ens_cost = @expression(
        model,
        p_rp_weight * sum(p_sc_prob[sc] * p_ens_cost * v_ens[sc, p] for sc in SC, p in P)
    )

    e_LR_cost = @expression(
        model,
        sum(
            p_lambda[SC[index], g] * (v_investment[SC[index], g] - v_investment[SC_lag[index], g]) for
            index in 1:length(SC), g in G
        )
    )

    # Objective function
    @objective(model, Min, e_investment_cost + e_variable_cost + e_ens_cost + e_LR_cost)

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
        get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] * v_investment[sc, g]
    )

    # Solve model
    optimize!(model)

    # Check if the model is optimal
    @assert is_solved_and_feasible(model) # Check if the model is optimal

    # Return results
    return model
end

"""
    add_cut(sets, params)

This function creates the cut to add in LR Decomposition algorithm.

# Arguments
- `first_stage_model::Model`: The first-stage model.
- `sets::Dict{Symbol, Any}`: A dictionary containing the sets of the problem.
- `params::Dict{Symbol, Any}`: A dictionary containing the parameters of the problem.
- `iteration::Int`: The current iteration number.

# Returns
- `nothing`
"""
function add_cut(first_stage_model, sets, params, p_investment, p_production, p_ens, iteration)
    # Extract sets
    SC = sets[:SC]
    G  = sets[:G]
    P  = sets[:P]

    # Create set scenario lag
    SC_lag = circshift(SC, 1)

    # Extract parameters
    p_investment_cost = params[:investment_cost]
    p_variable_cost   = params[:variable_cost]
    p_unit_capacity   = params[:unit_capacity]
    p_sc_prob         = params[:sc_prob]
    p_rp_weight       = params[:rp_weight]
    p_ens_cost        = params[:ens_cost]

    # Expressions
    e_investment_cost = @expression(
        first_stage_model,
        sum(
            p_sc_prob[sc] * p_investment_cost[g] * p_unit_capacity[g] * p_investment[sc, g] for
            sc in SC, g in G
        )
    )

    e_variable_cost = @expression(
        first_stage_model,
        p_rp_weight * sum(
            p_sc_prob[sc] * p_variable_cost[g] * p_production[sc, g, p] for sc in SC, g in G,
            p in P
        )
    )

    e_ens_cost = @expression(
        first_stage_model,
        p_rp_weight * sum(p_sc_prob[sc] * p_ens_cost * p_ens[sc, p] for sc in SC, p in P)
    )

    e_LR_cost = @expression(
        first_stage_model,
        sum(
            first_stage_model[:v_lambda][SC[index], g] *
            (p_investment[SC[index], g] - p_investment[SC_lag[index], g]) for
            index in 1:length(SC), g in G
        )
    )

    # Add LR cut
    first_stage_model[Symbol("cut_iter_$(iteration)")] = @constraint(
        first_stage_model,
        base_name = "cut_iter_$(iteration)",
        first_stage_model[:v_theta] ≤ e_investment_cost + e_variable_cost + e_ens_cost + e_LR_cost
    )

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
