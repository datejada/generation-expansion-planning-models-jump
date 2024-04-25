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
    G = sets[:G]

    # Extract parameters
    p_investment_cost = params[:investment_cost]
    p_unit_capacity   = params[:unit_capacity]

    # Scalar values
    M = -1000 # to avoid the problem is unbounded 

    # Model
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0, "output_flag" => false),
    )

    # Variables
    @variable(model, 0 ≤ v_investment[G], Int) #number of installed generation units [N]
    @variable(model, v_theta ≥ M)                    #Benders' cut

    # Expressions
    e_investment_cost = @expression(
        model,
        sum(p_investment_cost[g] * p_unit_capacity[g] * v_investment[g] for g in G)
    )

    # Objective function
    @objective(model, Min, e_investment_cost + v_theta)

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
    G = sets[:G]
    P = sets[:P]

    # Extract parameters
    p_max_avai      = params[:max_availability]
    p_min_avai      = params[:min_availability]
    p_demand        = params[:demand]
    p_variable_cost = params[:variable_cost]
    p_unit_capacity = params[:unit_capacity]
    p_is_renewable  = params[:is_renewable]
    p_avai_factor   = params[:availability_factor]
    p_rp_weight     = params[:rp_weight]
    p_ens_cost      = params[:ens_cost]
    p_exc_cost      = params[:exc_cost]
    p_uncert_budget = params[:uncertainty_budget]

    # Create subsets
    R = [r for r in G if p_is_renewable[r]]

    # Scalar values
    p_BigM = 5000 # Big M value

    # Model
    model = Model(
        optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0, "output_flag" => false),
    )

    # Positive Variables
    @variable(model, 0 ≤ v_production[G, P])                #production [MW] 
    @variable(model, 0 ≤ v_ens[p in P] ≤ p_demand[p])       #energy not supplied [MW]
    @variable(model, 0 ≤ v_exc[p in P] ≤ p_demand[p])       #excess of energy [MW]
    @variable(model, 0 ≤ v_avai[g in G] ≤ p_avai_factor[g]) #uncertain generation availability[p.u.]
    @variable(model, 0 ≤ v_dual_max_prod[G, P]) #dual variable of max production constraint   [kEUR/MWh]
    @variable(model, 0 ≤ v_dual_min_prod[G, P]) #dual variable of min production constraint   [kEUR/MWh]  
    @variable(model, 0 ≤ v_dual_max_ens[P])     #dual variable of max ENS        constraint   [kEUR/MWh]
    @variable(model, 0 ≤ v_dual_min_ens[P])     #dual variable of min ENS        constraint   [kEUR/MWh]
    @variable(model, 0 ≤ v_dual_max_exc[P])     #dual variable of max EXC        constraint   [kEUR/MWh]
    @variable(model, 0 ≤ v_dual_min_exc[P])     #dual variable of min EXC        constraint   [kEUR/MWh]

    # Free Variables
    @variable(model, v_dual_balance[P])         #dual variable of balance        constraint   [kEUR/MWh]

    # Binary Variables
    @variable(model, v_csc_max_prod[G, P], Bin) #aux binary for Complementary Slackness Condition of max production
    @variable(model, v_csc_min_prod[G, P], Bin) #aux binary for Complementary Slackness Condition of min production
    @variable(model, v_csc_max_ens[P], Bin)     #aux binary for Complementary Slackness Condition of max ENS       
    @variable(model, v_csc_min_ens[P], Bin)     #aux binary for Complementary Slackness Condition of min ENS       
    @variable(model, v_csc_max_exc[P], Bin)     #aux binary for Complementary Slackness Condition of max EXC       
    @variable(model, v_csc_min_exc[P], Bin)     #aux binary for Complementary Slackness Condition of min EXC       

    # Expressions
    e_variable_cost = @expression(
        model,
        p_rp_weight * sum(p_variable_cost[g] * v_production[g, p] for g in G, p in P)
    )

    e_ens_cost = @expression(model, p_rp_weight * sum(p_ens_cost * v_ens[p] for p in P))

    e_exc_cost = @expression(model, p_rp_weight * sum(p_exc_cost * v_exc[p] for p in P))

    # Objective function
    @objective(model, Max, e_variable_cost + e_ens_cost + e_exc_cost)

    # Constraints
    # - uncertainty budget (only for renewable sources)
    @constraint(
        model,
        c_uncertainty_budget,
        sum(p_avai_factor[r] - v_avai[r] for r in R) ≤
        sum(p_avai_factor[r] for r in R) * p_uncert_budget
    )

    # - fix the v_avai for non-renewable sources
    for g in G
        if !p_is_renewable[g]
            fix(v_avai[g], p_avai_factor[g]; force = true)
        end
    end

    # - balance equation
    @constraint(
        model,
        c_balance[p in P],
        sum(v_production[g, p] for g in G) + v_ens[p] == p_demand[p] + v_exc[p]
    )

    # - maximum generation
    @constraint(
        model,
        c_max_prod[g in G, p in P],
        v_production[g, p] <=
        (
            get(p_min_avai, (g, p), 0.0) +
            (get(p_max_avai, (g, p), 1.0) - get(p_min_avai, (g, p), 0.0)) * v_avai[g]
        ) *
        p_unit_capacity[g] *
        p_investment[g]
    )

    # - minimum generation
    @constraint(model, c_min_prod[g in G, p in P], -v_production[g, p] ≤ 0)

    # - maximum ENS
    @constraint(model, c_max_ens[p in P], v_ens[p] ≤ p_demand[p])

    # - minimum ENS
    @constraint(model, c_min_ens[p in P], -v_ens[p] ≤ 0)

    # - maximum EXC
    @constraint(model, c_max_exc[p in P], v_exc[p] ≤ p_demand[p])

    # - minimum EXC
    @constraint(model, c_min_exc[p in P], -v_exc[p] ≤ 0)

    # - dual generation constraints
    @constraint(
        model,
        c_dual_prod[g in G, p in P],
        p_rp_weight * p_variable_cost[g] - v_dual_balance[p] + v_dual_max_prod[g, p] -
        v_dual_min_prod[g, p] == 0
    )

    # - dual ENS constraints
    @constraint(
        model,
        c_dual_ens[p in P],
        p_rp_weight * p_ens_cost - v_dual_balance[p] + v_dual_max_ens[p] - v_dual_min_ens[p] == 0
    )

    # - dual EXC constraints
    @constraint(
        model,
        c_dual_exc[p in P],
        p_rp_weight * p_exc_cost + v_dual_balance[p] + v_dual_max_exc[p] - v_dual_min_exc[p] == 0
    )

    # - complementary slackness conditions max production part a
    @constraint(
        model,
        c_csc_max_prod_a[g in G, p in P],
        v_dual_max_prod[g, p] ≤ p_BigM * v_csc_max_prod[g, p]
    )

    # - complementary slackness conditions max production part b
    @constraint(
        model,
        c_csc_max_prod_b[g in G, p in P],
        (
            get(p_min_avai, (g, p), 0.0) +
            (get(p_max_avai, (g, p), 1.0) - get(p_min_avai, (g, p), 0.0)) * v_avai[g]
        ) *
        p_unit_capacity[g] *
        p_investment[g] - v_production[g, p] ≤ p_BigM * (1 - v_csc_max_prod[g, p])
    )

    # - complementary slackness conditions min production part a
    @constraint(
        model,
        c_csc_min_prod_a[g in G, p in P],
        v_dual_min_prod[g, p] ≤ p_BigM * v_csc_min_prod[g, p]
    )

    # - complementary slackness conditions min production part b
    @constraint(
        model,
        c_csc_min_prod_b[g in G, p in P],
        v_production[g, p] ≤ p_BigM * (1 - v_csc_min_prod[g, p])
    )

    # - complementary slackness conditions max ENS part a
    @constraint(model, c_csc_max_ens_a[p in P], v_dual_max_ens[p] ≤ p_BigM * v_csc_max_ens[p])

    # - complementary slackness conditions max ENS part b
    @constraint(
        model,
        c_csc_max_ens_b[p in P],
        p_demand[p] - v_ens[p] ≤ p_BigM * (1 - v_csc_max_ens[p])
    )

    # - complementary slackness conditions min ENS part a
    @constraint(model, c_csc_min_ens_a[p in P], v_dual_min_ens[p] ≤ p_BigM * v_csc_min_ens[p])

    # - complementary slackness conditions min ENS part b
    @constraint(model, c_csc_min_ens_b[p in P], v_ens[p] ≤ p_BigM * (1 - v_csc_min_ens[p]))

    # - complementary slackness conditions max EXC part a
    @constraint(model, c_csc_max_exc_a[p in P], v_dual_max_exc[p] ≤ p_BigM * v_csc_max_exc[p])

    # - complementary slackness conditions max EXC part b
    @constraint(
        model,
        c_csc_max_exc_b[p in P],
        p_demand[p] - v_exc[p] ≤ p_BigM * (1 - v_csc_max_exc[p])
    )

    # - complementary slackness conditions min EXC part a
    @constraint(model, c_csc_min_exc_a[p in P], v_dual_min_exc[p] ≤ p_BigM * v_csc_min_exc[p])

    # - complementary slackness conditions min EXC part b
    @constraint(model, c_csc_min_exc_b[p in P], v_exc[p] ≤ p_BigM * (1 - v_csc_min_exc[p]))

    # Solve model
    optimize!(model)

    # Check if the model is optimal
    @assert is_solved_and_feasible(model) # We don't check the dual solution since the cut is added in the primal information

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
function add_cut(first_stage_model, subproblem, sets, params, iteration)
    # Here the cut is added in the primal information since the subproblem has binary variables
    # and the dual information is not available for MIPs

    # Extract sets
    G = sets[:G]
    P = sets[:P]

    # Extract parameters
    p_max_avai      = params[:max_availability]
    p_min_avai      = params[:min_availability]
    p_variable_cost = params[:variable_cost]
    p_unit_capacity = params[:unit_capacity]
    p_rp_weight     = params[:rp_weight]
    p_demand        = params[:demand]
    p_ens_cost      = params[:ens_cost]
    p_exc_cost      = params[:exc_cost]

    # Get parameters from the subproblem
    p_avai = value.(subproblem[:v_avai])

    # Add variables
    v_production_iter = @variable(
        first_stage_model,
        base_name = "v_production_iter_$(iteration)",
        [G, P],
        lower_bound = 0
    ) #production [MW] 
    v_ens_iter = @variable(
        first_stage_model,
        base_name = "v_ens_iter_$(iteration)",
        [p in P],
        lower_bound = 0,
        upper_bound = p_demand[p]
    ) #energy not supplied [MW]
    v_exc_iter = @variable(
        first_stage_model,
        base_name = "v_exc_iter_$(iteration)",
        [p in P],
        lower_bound = 0,
        upper_bound = p_demand[p]
    ) #excess of energy [MW]

    # Add Benders' cut
    @constraint(
        first_stage_model,
        base_name = "cut_iter_$(iteration)",
        first_stage_model[:v_theta] >=
        p_rp_weight * (
            sum(p_variable_cost[g] * v_production_iter[g, p] for g in G, p in P) +
            sum(p_ens_cost * v_ens_iter[p] + p_exc_cost * v_exc_iter[p] for p in P)
        )
    )

    # balance equation
    @constraint(
        first_stage_model,
        base_name = "balance_iter_$(iteration)",
        [p in P],
        sum(v_production_iter[g, p] for g in G) + v_ens_iter[p] == p_demand[p] + v_exc_iter[p]
    )

    # maximum generation
    @constraint(
        first_stage_model,
        base_name = "max_prod_iter_$(iteration)",
        [g in G, p in P],
        v_production_iter[g, p] <=
        (
            get(p_min_avai, (g, p), 0.0) +
            (get(p_max_avai, (g, p), 1.0) - get(p_min_avai, (g, p), 0.0)) * p_avai[g]
        ) *
        p_unit_capacity[g] *
        first_stage_model[:v_investment][g]
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
