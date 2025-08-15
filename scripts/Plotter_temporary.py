# === INPUT ===
import sys
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt



# === 🔹 1. MANUAL FILE PATH ===
file_path = "results/model_gbr_penalty_imports_medium.nc"
ds = xr.open_dataset(file_path)

# === CONVERSION FACTORS ===
# Euro-Calliope default units:
# power: 100 GW → convert to GW multiply by 1e2
# energy: 100 GWh → convert to GWh multiply by 1e2
P_FACTOR = 1e2   # to GW
E_FACTOR = 1e2  # to GWh (per timestep)

# === SUPPLY & DEMAND TECH LISTS ===
supply_techs = [
    "electrified_biofuel", "historic_electrified_heat", "hydro_reservoir",
    "hydro_run_of_river", "nuclear", "open_field_pv", "roof_mounted_pv",
    "wind_offshore", "wind_onshore_competing", "wind_onshore_monopoly", "import_electricity"
]
demand_techs = ["demand_elec", "demand_heat_electrified", "export_electricity"]

# === TIME CONVERSION ===
time_index = pd.to_datetime(ds["timesteps"].values)
weekly = pd.Grouper(freq="W", sort=True)

# === HELPER: AGGREGATE FLOWS ===
def aggregate_flow(variable, techs, scale):
    """Returns a weekly aggregated dataframe for given variable and tech list."""
    da = ds[variable].sel(techs=techs).sum(dim=["nodes", "carriers"]) * scale
    df = da.to_pandas()
    df.index = time_index
    return df.groupby(weekly).sum()

# === 1. FLOW_OUT per supply tech ===
flow_out_weekly = {}
for tech in supply_techs:
    if tech in ds.techs:
        da = ds["flow_out"].sel(techs=tech).sum(dim=["nodes", "carriers"]) * P_FACTOR
        df = da.to_pandas()
        df.index = time_index
        flow_out_weekly[tech] = df.groupby(weekly).mean()

df_flow_out = pd.DataFrame(flow_out_weekly)

# === 2. FLOW_IN per demand tech ===
flow_in_weekly = {}
for tech in demand_techs:
    if tech in ds.techs:
        da = ds["flow_in"].sel(techs=tech).sum(dim=["nodes", "carriers"]) * P_FACTOR
        df = da.to_pandas()
        df.index = time_index
        flow_in_weekly[tech] = df.groupby(weekly).mean()

df_flow_in = pd.DataFrame(flow_in_weekly)

# === 3. SOURCE_USE for supply techs ===
source_use_weekly = {}
if "source_use" in ds:
    for tech in supply_techs:
        if tech in ds.techs:
            da = ds["source_use"].sel(techs=tech).sum(dim=["nodes"]) * E_FACTOR
            df = da.to_pandas()
            df.index = time_index
            source_use_weekly[tech] = df.groupby(weekly).sum()

df_source_use = pd.DataFrame(source_use_weekly)

# === 🔹 INSTALLED CAPACITY TABLE (flow_cap) ===
if "flow_cap" in ds:
    cap = ds["flow_cap"].sel(techs=supply_techs).sum(dim="carriers", min_count=1)
    df_cap = cap.to_series().reset_index()
    df_cap.columns = ["nodes", "techs", "value"]

    # Pivot into a nice table: techs as rows, nodes as columns
    cap_table = df_cap.pivot(index="techs", columns="nodes", values="value").fillna(0)
    
    print("\n📊 === Installed Capacities [100 GW] ===")
    print(cap_table.round(3))

    # (Optional) save to CSV for later analysis
    cap_table.to_csv("results/installed_capacities.csv")
    print("✅ Installed capacities saved to results/installed_capacities.csv")
else:
    print("⚠️ flow_cap not found in dataset.")



# === 4. PLOT RESULTS ===
for df, title, ylabel in [
    (df_flow_out, "Weekly Average Power Output per Supply Tech", "Power [MW]"),
    (df_flow_in, "Weekly Average Power Input per Demand Tech", "Power [MW]"),
    (df_source_use, "Weekly Energy Supplied (Source Use)", "Energy [MWh]"),
]:
    plt.figure(figsize=(14, 6))
    df.plot(ax=plt.gca())
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("Week")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()

# === 5. CHECK CONSTRAINTS ===
imports = (ds["flow_out"].sel(techs="import_electricity").sum() * E_FACTOR).item()
exports = (ds["flow_in"].sel(techs="export_electricity").sum() * E_FACTOR).item()
gen_total = (ds["flow_out"].sel(techs=supply_techs).sum() * E_FACTOR).item()

net_import = imports - exports

# Conditions:
cond_net_exporter = imports <= 0.95 * exports
cond_max_share = net_import <= 0.05 * gen_total

print("\n=== Constraint Check ===")
print(f"Total Imports: {imports/1e3:.2f} GWh")
print(f"Total Exports: {exports/1e3:.2f} GWh")
print(f"Total Generation: {gen_total/1e3:.2f} GWh")
print(f"Net Imports: {net_import/1e3:.2f} GWh")
print(f"Net Exporter Constraint Satisfied? {'✅' if cond_net_exporter else '❌'}")
print(f"Max Import Share Constraint Satisfied? {'✅' if cond_max_share else '❌'}")
