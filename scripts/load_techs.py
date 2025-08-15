#import calliope

#model = calliope.Model("models/ehighways/model.yaml", scenario="only_gbr")
#techs = model.inputs.name.to_series()
#print(techs)

import calliope

# === Load your model with the correct scenario ===
model = calliope.Model("models/ehighways/model.yaml", scenario="only_gbr")

# === Extract all technologies actually instantiated in the model ===
all_techs = model.inputs.name.to_series()
print("✅ All instantiated technologies:\n", all_techs)

# === Separate by type if desired (e.g., supply, demand, transmission) ===
base_techs = model.inputs.base_tech.to_series()

gen_techs = base_techs[base_techs == "supply"].index.tolist()
demand_techs = base_techs[base_techs == "demand"].index.tolist()
transmission_techs = base_techs[base_techs == "transmission"].index.tolist()

print("\n🔹 Generation technologies:\n", gen_techs)
print("\n🔹 Demand technologies:\n", demand_techs)
print("\n🔹 Transmission technologies:\n", transmission_techs)

# === Optional: Write to a CSV to use directly for penalties ===
import pandas as pd
pd.DataFrame({
    "techs": all_techs.index,
    "base_tech": base_techs.values,
    "name": all_techs.values
}).to_csv("tech_list_for_penalties.csv", index=False)

print("\n✅ Technology list exported to tech_list_for_penalties.csv")
