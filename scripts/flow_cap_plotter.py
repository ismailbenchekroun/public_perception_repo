import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Scenario Files ===
scenarios = {
    "Cost Optimal": "results/model_gbr_imports.nc",
    "Low Penalty": "results/model_gbr_penalty_imports_low.nc",
    "Medium Penalty": "results/model_gbr_penalty_imports_medium.nc",
    "High Penalty": "results/model_gbr_penalty_imports_high.nc",
}

# === Setup ===
P_FACTOR = 1e2
supply_techs = [
    "electrified_biofuel", "historic_electrified_heat", "hydro_reservoir",
    "hydro_run_of_river", "nuclear", "open_field_pv", "roof_mounted_pv",
    "wind_offshore", "wind_onshore_competing", "wind_onshore_monopoly", "import_electricity", "battery", "pumped_hydro"
]

# === Load Installed Capacity Data ===
cap_dfs = []
for scenario_name, filepath in scenarios.items():
    ds = xr.open_dataset(filepath)
    cap = ds["flow_cap"].sel(techs=supply_techs).sum(dim="carriers", min_count=1)
    df = cap.to_series().reset_index()
    df.columns = ["nodes", "techs", "value"]
    df["scenario"] = scenario_name
    cap_dfs.append(df)

# === Combine All Scenarios ===
df_all = pd.concat(cap_dfs)

# === Pivot Table: tech + node ===
pivot = df_all.pivot_table(index=["techs", "nodes"], columns="scenario", values="value", fill_value=0).reset_index()

# === Compute Differences vs Cost Optimal ===
for scen in scenarios:
    if scen != "Cost Optimal":
        pivot[f"Δ {scen}"] = pivot[scen] - pivot["Cost Optimal"]

# === Make Output Folder ===
output_dir = "results/plots/capacity_delta_by_tech"
os.makedirs(output_dir, exist_ok=True)

# === Plot per Tech ===
for tech in pivot["techs"].unique():
    df_tech = pivot[pivot["techs"] == tech].copy()
    df_melt = df_tech.melt(
        id_vars=["nodes"],
        value_vars=[f"Δ {s}" for s in scenarios if s != "Cost Optimal"],
        var_name="Scenario",
        value_name="Δ Capacity (GW)"
    )

    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=df_melt, x="nodes", y="Δ Capacity (GW)", hue="Scenario", palette="Set2")
    plt.axhline(0, color="black", linestyle="--", linewidth=1)
    plt.title(f"Δ Installed Capacity vs Cost Optimal – Tech: {tech}")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Annotate bars
    for bar in ax.patches:
        if abs(bar.get_height()) > 0.01:
            ax.annotate(f'{bar.get_height():.2f}', 
                        (bar.get_x() + bar.get_width() / 2., bar.get_height()), 
                        ha='center', 
                        va='bottom' if bar.get_height() > 0 else 'top',
                        fontsize=8, rotation=90)

    # Save plot
    out_file = os.path.join(output_dir, f"delta_cap_{tech}.png")
    plt.savefig(out_file, dpi=300)
    plt.close()

print(f"✅ All plots grouped by tech saved in: {output_dir}")
