# Static Robust Optimization (Scenario-based) - Generation Expansion Planning

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
    p_availability    = params[:availability]
    p_demand          = params[:demand]
    p_investment_cost = params[:investment_cost]
    p_variable_cost   = params[:variable_cost]
    p_unit_capacity   = params[:unit_capacity]
    p_rp_weight       = params[:rp_weight]
    p_ens_cost        = params[:ens_cost]

    # Model
    model = Model(optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.0))

    # Variables
    @variable(model, 0 ≤ v_production[SC, G, P])          #production [MW] 
    @variable(model, 0 ≤ v_ens[SC, p in P] ≤ p_demand[p]) #energy not supplied [MW]
    @variable(model, 0 ≤ v_investment[G], Int)            #number of installed generation units [N]
    @variable(model, β)                                   #worst operation cost [kEUR]

    # Expressions
    e_investment_cost = @expression(
        model,
        sum(p_investment_cost[g] * p_unit_capacity[g] * v_investment[g] for g in G)
    )

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

    # Objective function
    @objective(model, Min, e_investment_cost + β)

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
        v_production[sc, g, p] ≤
        get(p_availability, (sc, g, p), 1.0) * p_unit_capacity[g] * v_investment[g]
    )

    # - operation cost per scenario
    @constraint(model, c_oper_cost[sc in SC], β ≥ e_variable_cost[sc] + e_ens_cost[sc])

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
