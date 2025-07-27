configfile: "config/default.yaml"

# ---------------------- Model Runs ----------------------

rule run_eu_model:
    message: "Run the full Calliope base model."
    input:
        model_yaml = "ehighways/model.yaml"
    output:
        nc_file = "results/model_eu.nc"
    conda: "environment.yaml"
    shell:
        "calliope run {input.model_yaml} --save_netcdf {output.nc_file}"

rule run_gbr_model:
    message: "Run Calliope model using GBR-only scenario."
    input:
        model_yaml = "../../ehighways/model.yaml"
    output:
        nc_file = "results/model_gbr_only.nc"
    conda: "environment.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario only_gbr --save_netcdf {output.nc_file}"

rule run_gbr_model_with_penalty:
    message: "Run GBR model with penalty factors (uses the penalty_factors override)."
    input:
        model_yaml = "../../ehighways/model.yaml",
        penalty_csv = "../../data/outputs/penalty_factors_computed.csv"
    output:
        nc_file = "results/model_gbr_with_penalty.nc"
    conda: "environment.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario only_gbr,penalty_factors --save_netcdf {output.nc_file}"

rule run_gbr_model_with_penalty_and_imports:
    message: "Run GBR model with penalty factors and imports."
    input:
        model_yaml = "models/ehighways/model.yaml",
        penalty_csv = "data/outputs/penalty_factors_computed.csv"
    output:
        nc_file = "results/model_gbr_penalty_imports.nc"
    conda: "environment.yaml"
    shell:
        "calliope run {input.model_yaml} --scenario only_gbr,penalty_factors,add_uk_import_export --save_netcdf {output.nc_file}"

# ---------------------- Supporting Rules ----------------------

rule visualise_gbr_model:
    message: "Launch Calliope's visualisation tool for GBR-only model."
    input:
        nc_file = rules.run_gbr_model.output.nc_file
    conda: "environment.yaml"
    shell:
        "calligraph {input.nc_file}"

rule visualise_gbr_model_with_penalty:
    message: "Launch Calliope's visualisation tool for GBR-only model with penalty factors."
    input:
        nc_file = rules.run_gbr_model_with_penalty.output.nc_file
    conda: "environment.yaml"
    shell:
        "calligraph {input.nc_file}"

# ---------------------- Data Preprocessing ----------------------

rule compute_penalty_factors:
    message: "Compute penalty factors (scenario chosen via --config scenario=low|medium|high)"
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
        "results/model_gbr_with_penalty.nc",
        "results/model_gbr_penalty_imports.nc"
