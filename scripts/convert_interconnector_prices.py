import pandas as pd

# Load input file (with renamed columns: timesteps, import_electricity, export_electricity)
df = pd.read_csv(snakemake.input[0])

# Convert timestamps to ISO8601 format
df["timesteps"] = pd.to_datetime(df["timesteps"], dayfirst=True).dt.strftime("%Y-%m-%d %H:%M:%S")

# Build long-format dataframe
records = []
for _, row in df.iterrows():
    records.append([row["timesteps"], "import_electricity", "cost_flow_out", row["import_electricity"]])
    records.append([row["timesteps"], "export_electricity", "cost_flow_in", row["export_electricity"]])

# Create output dataframe
long_df = pd.DataFrame(records, columns=["timesteps", "techs", "parameters", "value"])

# Save to output
long_df.to_csv(snakemake.output[0], index=False)
print(f"✅ Correctly formatted interconnector prices saved to {snakemake.output[0]}")
