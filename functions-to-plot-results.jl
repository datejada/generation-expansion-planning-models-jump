"""
    plot_two_stage_investment(model)

Plot the investment in generation units for the two-stage problem.

# Arguments
- `model`: The model containing the investment variable.
- `sets`: A dictionary containing the sets of the problem.
- `params`: A dictionary containing the parameters of the problem.

"""
function plot_two_stage_investment(model, sets, params)
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
    plot_multi_stage_investment(model)

Plot the investment in generation units for the multi-stage problem.

# Arguments
- `model`: The model containing the investment variable.
- `params`: A dictionary containing the parameters of the model.

"""
function plot_multi_stage_investment(model, params)
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
    plot_two_stage_production(model, params)

Plot the production of generators and the demand for each scenario.

# Arguments
- `model`: The optimization model containing the production variables.
- `params`: A dictionary containing the parameters of the model.

# Returns
- `p`: A plot object showing the production and demand for each scenario.

"""
function plot_two_stage_production(model, params)
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
    plot_investment_per_iteration(df_investment_per_iteration)

Plot the investment per iteration.

# Arguments
- `df_investment_per_iteration`: A DataFrame containing the investment per iteration.

# Returns
- `p`: A plot object showing the investment per iteration.

"""
function plot_investment_per_iteration(df_investment_per_iteration)
    p = @df df_investment_per_iteration plot(
        :iteration,
        :investment,
        group = :generator,
        xlabel = "Iteration",
        ylabel = "Investment [units]",
        title = "Investment per Iteration",
        lw = 2,
        xticks = 0:1:MAXIMUM_ITERATIONS,
    )
    return p
end