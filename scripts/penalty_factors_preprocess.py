import pandas as pd
import numpy as np

# User-defined parameters for this scenario
DEFAULT_PENALTY = 1.0

penalty_params = snakemake.params.penalty_params
# Load opposition data
df = pd.read_csv(snakemake.input[0], index_col=["nodes", "techs"])  # must contain techs, nodes, O_1, O_2

# Function to compute omega
def compute_omega(row, k, alpha, beta):
    O_1 = row['O_1']
    O_2 = row['O_2']
    if not pd.isna(O_1) and not pd.isna(O_2):
        return k * (1 - np.exp(-alpha * (1 - O_1/100)) * np.exp(-beta * (1 - O_2/100)))
    elif not pd.isna(O_1):
        return k * (1 - np.exp(-alpha * (1 - O_1/100)))
    elif not pd.isna(O_2):
        return k * (1 - np.exp(-beta * (1 - O_2/100)))
    else:
        return DEFAULT_PENALTY

# Apply the function
new_df = pd.DataFrame(DEFAULT_PENALTY, index=df.index, columns=penalty_params.keys())
for penalty_scenario, params in penalty_params.items():
    new_df[penalty_scenario] = df.apply(compute_omega, k=params["k"], alpha=params["alpha"], beta=params["beta"], axis=1)

# Save the result
new_df.to_csv(snakemake.output[0])
