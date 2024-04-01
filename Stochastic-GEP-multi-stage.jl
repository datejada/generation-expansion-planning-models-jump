# Multi-Stage Stochastic Generation Expansion Planning

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
