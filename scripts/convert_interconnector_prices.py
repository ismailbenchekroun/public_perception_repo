import pandas as pd

# Load input file from snakemake
df = pd.read_csv(snakemake.input[0])  # expects columns: timesteps, Import price, Export Price

# Build long-format dataframe for Calliope
records = []
for _, row in df.iterrows():
    records.append([row["timesteps"], "import_electricity", "cost_flow_out", row["Import price"]])
    records.append([row["timesteps"], "export_electricity", "cost_flow_in", row["Export price"]])

long_df = pd.DataFrame(records, columns=["timesteps", "techs", "parameters", "value"])

# ✅ Save the file in Calliope's expected format
long_df.to_csv(snakemake.output[0], index=False)

print(f"✅ Formatted interconnector prices saved to {snakemake.output[0]}")