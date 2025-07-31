import numpy as np
import pandas as pd

DEFAULT_PENALTY = 1.0
penalty_params = snakemake.params.penalty_params

# Load opposition data (must have: nodes, techs, O_1, O_2)
df = pd.read_csv(snakemake.input[0])


def compute_omega(row, k, alpha, beta):
    O_1, O_2 = row["O_1"], row["O_2"]
    if not pd.isna(O_1) and not pd.isna(O_2):
        return k * (
            1 - np.exp(-alpha * (1 - O_1 / 100)) * np.exp(-beta * (1 - O_2 / 100))
        )
    elif not pd.isna(O_1):
        return k * (1 - np.exp(-alpha * (1 - O_1 / 100)))
    elif not pd.isna(O_2):
        return k * (1 - np.exp(-beta * (1 - O_2 / 100)))
    else:
        return DEFAULT_PENALTY


# Calculate omega for each scenario and store in new_df
penalty_values = {}
for scenario_name, params in penalty_params.items():
    penalty_values[scenario_name] = df.apply(
        compute_omega, axis=1, k=params["k"], alpha=params["alpha"], beta=params["beta"]
    )
new_df = pd.concat([df[["techs", "nodes"]], pd.DataFrame(penalty_values)], axis=1)

# Write output
new_df.to_csv(snakemake.output[0], index=False)
print(f"✅ Penalty factors written: {snakemake.output[0]}")
