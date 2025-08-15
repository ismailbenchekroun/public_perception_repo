import pandas as pd

# === CONFIGURE THESE FILE PATHS ===
penalty_factors_techs_path = "data/outputs/penalty_factors_manual_techs.csv"
penalty_factors_links_path = "data/outputs/penalty_factors_manual_links.csv"

output_techs_path = "data/outputs/penalty_costs_techs.csv"
output_links_path = "data/outputs/penalty_costs_links.csv"

# === STEP 1: Provide base cost_flow_cap ((10,000 EUR2015/MW) for each technology

base_capex = {
    'open_field_pv': 30.125,
    'roof_mounted_pv': 88,
    'wind_onshore_competing': 96.3,
    'wind_onshore_monopoly': 96.3,
    'wind_offshore': 130,
    'nuclear': 367.2,
    'electrified_biofuel': 230,
    'battery': 8.558536,
    'pumped_hydro': 103.833246,
    'ac_ohl': 17.5,  # Transmission lines
}

# === STEP 2: Load penalty factor files ===
df_techs = pd.read_csv(penalty_factors_techs_path)
df_links = pd.read_csv(penalty_factors_links_path)

# === STEP 3: Compute penalty cost from omega using: CAPEX * (1/ω - 1)
def compute_penalty_costs(df, capex_lookup):
    scenarios = ['low', 'medium', 'high']
    result = df.copy()
    for scenario in scenarios:
        result[scenario] = df.apply(
            lambda row: capex_lookup.get(row['techs'], 1000) * (1 / row[scenario] - 1)
            if row[scenario] != 0 else 0,
            axis=1
        )
    return result[['techs', 'nodes'] + scenarios] if 'nodes' in df.columns else result[['techs'] + scenarios]

# === STEP 4: Apply to techs and links
penalty_costs_techs = compute_penalty_costs(df_techs, base_capex)
penalty_costs_links = compute_penalty_costs(df_links, base_capex)

# === STEP 5: Save results
penalty_costs_techs.to_csv(output_techs_path, index=False)
penalty_costs_links.to_csv(output_links_path, index=False)

print(f"✅ Penalty costs saved to:\n- {output_techs_path}\n- {output_links_path}")
