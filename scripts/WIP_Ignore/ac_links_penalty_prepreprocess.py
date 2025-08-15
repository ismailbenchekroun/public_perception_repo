import pandas as pd

# === Input files ===
input_file = "data/inputs/penalty_factors_preinput_ac_only.csv"
output_file = "data/inputs/penalty_factors_nodes_to_links.csv"

# === Define AC links and their associated regions ===
links = {
    "link_GBR_1_GBR_2_ac_ohl": ["GBR_1", "GBR_2"],
    "link_GBR_1_GBR_3_ac_ohl": ["GBR_1", "GBR_3"],
    "link_GBR_2_GBR_3_ac_ohl": ["GBR_2", "GBR_3"],
    "link_GBR_3_GBR_4_ac_ohl": ["GBR_3", "GBR_4"],
    "link_GBR_4_GBR_5_ac_ohl": ["GBR_4", "GBR_5"]
}

# === Load regional penalties ===
df = pd.read_csv(input_file)

# Filter only rows for ac_ohl
df_ac = df[df["techs"] == "ac_ohl"].set_index("nodes")["O_2"]

# === Build penalty factors for each link ===
records = []
for link, (n1, n2) in links.items():
    if n1 in df_ac and n2 in df_ac:
        penalty = (df_ac[n1] + df_ac[n2]) / 2
        records.append([link, "", penalty])
    else:
        print(f"⚠️ Missing data for {link}, skipping...")

# === Save to CSV ===
out_df = pd.DataFrame(records, columns=["techs", "nodes", "O_1","O_2"])
out_df.to_csv(output_file, index=False)

print(f"✅ Link penalty factors written to {output_file}")
