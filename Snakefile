configfile: "config/default.yaml"

# ---------------------- Model Runs ----------------------

wildcard_constraints:
    penalty_scenario = "low|medium|high"

rule run_eu_model:
    message: "Run the full Calliope base model."
    input:
        model_yaml = "models/ehighways/model.yaml"
    output:
        nc_file = "results/model_eu.nc"
    conda: "envs/calliope.yaml"
    shell:
        "calliope run {input.model_yaml} --save_netcdf {output.nc_file}"

rule run_gbr_model:
    message: "Run Calliope model using GBR-only scenario."
    input:
        model_yaml = "models/ehighways/model.yaml"
    output:
        nc_file = "results/model_gbr_only.nc"
    conda: "envs/calliope.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario only_gbr --save_netcdf {output.nc_file}"

rule run_gbr_model_with_penalty:
    message: "Run GBR model with penalty factors with {wildcards.penalty_scenario} weighting (uses the penalty_factors override)."
    input:
        model_yaml = "models/ehighways/model.yaml",
        penalty_csv = "data/outputs/penalty_factors_computed.csv"
    output:
        nc_file = "results/model_gbr_with_penalty_{penalty_scenario}.nc"
    conda: "envs/calliope.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario {wildcards.penalty_scenario} --save_netcdf {output.nc_file}"

rule run_gbr_model_with_imports:
    message: "Run GBR model with imports (no penalties)."
    input:
        model_yaml = "models/ehighways/model.yaml"
    output:
        nc_file = "results/model_gbr_imports.nc"
    conda: "envs/calliope.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario only_gbr,add_uk_import_export --save_netcdf {output.nc_file}"


rule run_gbr_model_with_penalty_and_imports:
    message: "Run GBR model with penalty factors with {wildcards.penalty_scenario} weighting and imports."
    input:
        model_yaml = "models/ehighways/model.yaml",
        penalty_csv = "data/outputs/penalty_factors_computed.csv"
    output:
        nc_file = "results/model_gbr_penalty_imports_{penalty_scenario}.nc"
    conda: "envs/calliope.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario {wildcards.penalty_scenario},add_uk_import_export --save_netcdf {output.nc_file}"

# ---------------------- Supporting Rules ----------------------

rule visualise_gbr_model:
    message: "Launch Calliope's visualisation tool for GBR-only model."
    input:
        nc_file = rules.run_gbr_model.output.nc_file
    conda: "envs/calliope.yaml"
    shell:
        "calligraph {input.nc_file}"

rule visualise_gbr_model_with_penalty:
    message: "Launch Calliope's visualisation tool for GBR-only model with penalty factors."
    input:
        nc_file = rules.run_gbr_model_with_penalty.output.nc_file
    conda: "envs/calliope.yaml"
    shell:
        "calligraph {input.nc_file}"

rule visualise_gbr_model_with_imports:
    message: "Launch Calliope's visualisation tool for GBR model with imports (no penalties)."
    input:
        nc_file = rules.run_gbr_model_with_imports.output.nc_file
    conda: "envs/calliope.yaml"
    shell:
        "calligraph {input.nc_file}"

rule plot_imports_exports:
    message: "Plotting imports and exports from model results."
    input:
        nc_file="results/model_gbr_imports.nc"
    output:
        png="results/imports_exports_plot.png"
    shell:
        "python scripts/plot_imports_exports_v1.py {input.nc_file} {output.png}"

rule explore_results_manually:
    message: "Plotting imports and exports from model results."
    input:
        nc_file="results/model_gbr_imports.nc"
    output:
        plot="results/imports_exports_plot.png"
    shell:
        "python scripts/explore_results_manually.py {input.nc_file} {output.plot}"


# ---------------------- Data Preprocessing ----------------------

rule compute_penalty_factors:
    message: "Compute penalty factors for each scenario"
    input:
        "data/inputs/penalty_factors_input.csv"
    params:
        penalty_params = config["penalty_parameters"]
    output:
        "data/outputs/penalty_factors_computed.csv"
    script:
        "scripts/penalty_factors_preprocess.py"

rule convert_interconnector_prices:
    message: "Format interconnector prices CSV"
    input: "data/inputs/interconnector_prices_avg.csv"
    output: "data/outputs/interconnector_prices_formatted.csv"
    script: "scripts/convert_interconnector_prices.py"

# ---------------------- Build All ----------------------

rule all:
    input:
        "data/outputs/penalty_factors_computed.csv",
        expand(
            "results/model_gbr_with_penalty_{penalty_scenario}.nc",
            penalty_scenario=["low", "medium", "high"]
        ),
        expand(
            "results/model_gbr_penalty_imports_{penalty_scenario}.nc",
            penalty_scenario=["low", "medium", "high"]
        )
    default_target: True
