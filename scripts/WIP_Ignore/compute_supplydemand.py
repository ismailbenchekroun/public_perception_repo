import xarray as xr

# === Load your results ===
file_path = "results/model_gbr_imports.nc"
ds = xr.open_dataset(file_path)

# ✅ Check variable existence
if "sink_use_equals" not in ds:
    raise KeyError("`sink_use_equals` not found in dataset.")
if "source_use" not in ds:
    raise KeyError("`source_use` not found in dataset.")

# === Total demand energy ===
total_demand = ds["sink_use_equals"].sum(dim=["nodes", "techs", "timesteps"]).item()

# === Total supply energy ===
total_supply = ds["source_use"].sum(dim=["nodes", "techs", "timesteps"]).item()

# === Imports energy ===
imports = ds["source_use"].sel(techs="import_electricity").sum(dim=["nodes", "timesteps"]).item()

# === Print results ===
print(f"🔹 Total demand energy: {total_demand:.2f} MWh")
print(f"🔹 Total supply energy: {total_supply:.2f} MWh")
print(f"🔹 Imports energy: {imports:.2f} MWh")
print(f"🔹 Imports share: {(imports / total_supply * 100 if total_supply > 0 else 0):.2f}%")

# === Sanity check: demand ≈ supply? ===
if abs(total_demand - total_supply) / total_demand < 0.05:  # within 5%
    print("✅ Supply and demand are balanced (as expected).")
else:
    print("⚠️ Supply-demand mismatch! Check constraints or units.")
