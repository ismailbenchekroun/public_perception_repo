import sys
import calliope
import matplotlib.pyplot as plt

file_path = sys.argv[1]
output_plot = sys.argv[2]

model = calliope.read_netcdf(file_path)

print("\n=== Available results variables ===")
print(model.results.data_vars)

tech_import = "import_electricity"
tech_export = "export_electricity"

# ✅ Sum over nodes and carriers to get (timesteps,)
imports = model.results.flow_out.sel(techs=tech_import).sum(dim=["nodes", "carriers"])
exports = model.results.flow_in.sel(techs=tech_export).sum(dim=["nodes", "carriers"])

# ✅ Ensure proper alignment
time = imports['timesteps']

plt.figure(figsize=(12, 5))
plt.plot(time, imports, label="Imports")
plt.plot(time, exports, label="Exports")
plt.title("Electricity Imports and Exports")
plt.ylabel("Power [GW]")
plt.legend()
plt.tight_layout()
plt.savefig(output_plot)
print(f"✅ Plot saved at {output_plot}")
