import sys
import xarray as xr
import matplotlib.pyplot as plt

# ✅ Get the NetCDF file path from command line argument
file_path = sys.argv[1]  

ds = xr.open_dataset(file_path)

# ✅ Known variables in Calliope output
tech_import = "import_electricity"   # exists in your YAML
tech_export = "export_electricity"   # exists in your YAML

# ✅ Extract relevant data
imports_total = ds['carrier_prod'].sel(techs=tech_import).sum(dim="nodes")
exports_total = ds['carrier_con'].sel(techs=tech_export).sum(dim="nodes")

# ✅ Plot
plt.figure(figsize=(12, 5))
imports_total.plot(label="Imports")
exports_total.plot(label="Exports")
plt.title("Electricity Imports and Exports")
plt.ylabel("Power [GW]")
plt.legend()
plt.savefig(sys.argv[2], dpi=300)  # ✅ save output file
