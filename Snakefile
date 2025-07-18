configfile: "config/default.yaml"

rule run_eu_model:
    message: "Run the full Calliope base model."
    input:
        model_yaml = "ehighways/model.yaml"
    output:
        nc_file = "results/model_eu.nc"
    conda:
        "environment.yaml"
    shell:
         "calliope run ehighways/model.yaml --save_netcdf results/eu-model.nc"

rule run_gbr_model:
    message: "Run Calliope model using GBR-only scenario."
    input:
        model_yaml = "ehighways/model.yaml"
    output:
        nc_file = "results/model_gbr_only.nc"
    conda:
        "environment.yaml"
    shell:
        "calliope run {input} --scenario only_gbr --save_netcdf {output}" #use this for model name results 

rule run_gbr_model_with_penalty:
    message: "Run Calliope model using GBR-only scenario with penalty factors."
    input:
        model_yaml = "ehighways/model.yaml"
    output:
        nc_file = "results/model_gbr_with_penalty.nc"
    conda:
        "environment.yaml"
    shell:
        "calliope run {input} --scenario only_gbr penalty_factors --save_netcdf {output}"  # Note: Two scenarios chained

rule visualise_gbr_model:
    message: "Launch Calliope's visualisation tool for GBR-only model."
    input:
        nc_file = rules.run_gbr_model.output.nc_file
    conda:
        "environment.yaml"
    shell:
        "calligraph {input.nc_file}"

rule visualise_gbr_model_with_penalty:
    message: "Launch Calliope's visualisation tool for GBR-only model with penalty factors."
    input:
        nc_file = rules.run_gbr_model_with_penalty.output.nc_file
    conda:
        "environment.yaml"
    shell:
        "calligraph {input.nc_file}"

rule compute_penalty_factors:
    message: "compute penalty factors for the scenario"
    input: "data/inputs/penalty_factors_input.csv"
    params:
        penalty_params = config["penalty_parameters"]
    output: "data/outputs/penalty_factors_computed.csv"
    script: "scripts/penalty_factors_preprocess.py"


