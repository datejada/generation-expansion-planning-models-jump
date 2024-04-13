"""
save_two_stage_results_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""

function save_two_stage_results_to_csv(output_folder, model)
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

"""
save_two_stage_explicit_formulation_results_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""

function save_two_stage_explicit_formulation_results_to_csv(output_folder, model)
    CSV.write(
        joinpath(output_folder, "oGEP_Investment_Explicit.csv"),
        Containers.rowtable(value, model[:v_investment]; header = [:sc, :g, :instal_units]),
    )
    CSV.write(
        joinpath(output_folder, "oGEP_Production_Explicit.csv"),
        Containers.rowtable(value, model[:v_production]; header = [:sc, :g, :h, :production]),
    )

    CSV.write(
        joinpath(output_folder, "oGEP_ENS_Explicit.csv"),
        Containers.rowtable(value, model[:v_ens]; header = [:sc, :h, :ens]),
    )

    return nothing
end

"""
    save_multi_stage_results_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""

function save_multi_stage_results_to_csv(output_folder, model)
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

"""
    save_results_benders_first_stage_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""
function save_results_benders_first_stage_to_csv(output_folder, model)
    CSV.write(
        joinpath(output_folder, "oGEP_Investment.csv"),
        Containers.rowtable(value, model[:v_investment]; header = [:g, :instal_units]),
    )

    return nothing
end

"""
    save_results_benders_subproblem_to_csv(output_folder, model)

Save the results of the model to CSV files.

# Arguments
- `output_folder`: The folder path where the CSV files will be saved.
- `model`: The model containing the results to be saved.
"""
function save_results_benders_subproblem_to_csv(output_folder, model)
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
