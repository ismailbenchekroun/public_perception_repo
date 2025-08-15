import pandas as pd
import numpy as np

# === INPUTS & OUTPUTS ===
input_regional = str(snakemake.input.regional)      # Original input penalties (techs + ac_ohl)
output_techs   = str(snakemake.output.techs)        # Output: processed penalties for techs
output_links   = str(snakemake.output.links)        # Output: processed penalties for links
scenarios      = snakemake.params.penalty_params    # Scenario parameters dict

DEFAULT_MULTIPLIER = 1.0  # if no info at all

# ──────────────────────────────────────────────────────────────────────────────
# 0) LOAD + SPLIT (unchanged logic)
# ──────────────────────────────────────────────────────────────────────────────
df = pd.read_csv(input_regional)

df_techs    = df[df["techs"] != "ac_ohl"].copy()
df_ac_nodes = df[df["techs"] == "ac_ohl"].set_index("nodes")["O_2"]

# Links dictionary (kept as you wrote it)
links = {
    "link_GBR_1_GBR_2_ac_ohl": ["GBR_1", "GBR_2"],
    "link_GBR_1_GBR_3_ac_ohl": ["GBR_1", "GBR_3"],
    "link_GBR_2_GBR_3_ac_ohl": ["GBR_2", "GBR_3"],
    "link_GBR_3_GBR_4_ac_ohl": ["GBR_3", "GBR_4"],
    "link_GBR_4_GBR_5_ac_ohl": ["GBR_4", "GBR_5"],
}

# Build link table: O2 = average of node O2; O1 unknown here
link_records = []
for link, (n1, n2) in links.items():
    if n1 in df_ac_nodes and n2 in df_ac_nodes:
        avg_O2 = (df_ac_nodes[n1] + df_ac_nodes[n2]) / 2  # still in [0..100]%
        link_records.append([link, np.nan, avg_O2])
    else:
        print(f"⚠️ Missing AC OHL data for {link}, skipping...")

df_links = pd.DataFrame(link_records, columns=["techs", "O_1", "O_2"])

# ──────────────────────────────────────────────────────────────────────────────
# 1) IMPUTE MISSING O1/O2 (kept, but simplified)
#    - Inputs are % in [0..100]. We'll convert to [0..1] inside the functions.
# ──────────────────────────────────────────────────────────────────────────────
df_techs["O_1_filled"] = df_techs["O_1"]
df_techs["O_2_filled"] = df_techs["O_2"]

for i, row in df_techs.iterrows():
    if pd.isna(row["O_2"]) and pd.notna(row["O_1"]):
        df_techs.at[i, "O_2_filled"] = min(100.0, 3 * row["O_1"])
    elif pd.isna(row["O_1"]) and pd.notna(row["O_2"]):
        df_techs.at[i, "O_1_filled"] = row["O_2"] / 3.0

# ──────────────────────────────────────────────────────────────────────────────
# 2) CALIBRATION UTILITIES (NEW)
#    We support:
#      - 'early' curve:   G = exp(-α (O1^p + λ O2^q)), 0<p,q<1
#      - 'late'  curve:   G = [1 - exp(- (α W1 + β W2))] / [1 - exp(- (α+β))], W=1-O
#    Then wrap M = M_min + (1 - M_min) * G
#    Floors: techs -> M_min = 0.5 ; links -> M_min = 0.2
#    β is always 1.5 * α (enforced below).
# ──────────────────────────────────────────────────────────────────────────────

def _early_G(O1, O2, alpha, p=0.6, q=0.6, lam=1.5):
    """Early-bite: O in [0,1]; missing terms are treated as 0."""
    o1 = np.clip(0.0 if np.isnan(O1) else O1 / 100.0, 0.0, 1.0)
    o2 = np.clip(0.0 if np.isnan(O2) else O2 / 100.0, 0.0, 1.0)
    val = (o1 ** p) + lam * (o2 ** q)
    return np.exp(-alpha * val)

def _late_G(O1, O2, alpha, beta=None):
    """Late-bite (Koecklin-like): W=1-O; normalised so G(1,1)=1, G(0,0)=0."""
    if beta is None:
        beta = 1.5 * alpha  # enforce β = 1.5 α

    o1 = np.clip(0.0 if np.isnan(O1) else O1 / 100.0, 0.0, 1.0)
    o2 = np.clip(0.0 if np.isnan(O2) else O2 / 100.0, 0.0, 1.0)
    w1 = 1.0 - o1
    w2 = 1.0 - o2

    num   = 1.0 - np.exp(- (alpha * w1 + beta * w2))
    denom = 1.0 - np.exp(- (alpha + beta))
    # Guard against denom ~ 0
    if denom <= 1e-12:
        return num  # effectively unnormalised if α+β tiny
    return num / denom

def _wrap_multiplier(G, M_min):
    """Apply floor wrapper, returns M in [M_min, 1]."""
    G = np.clip(G, 0.0, 1.0)
    return M_min + (1.0 - M_min) * G

def compute_multiplier(O1, O2, asset_class, curve, params):
    """
    Compute the cost multiplier M given:
      - O1, O2 in % (0..100) or NaN
      - asset_class: 'tech' or 'link' (decides the floor)
      - curve: 'early' or 'late'
      - params: dict with curve params; β is forced to 1.5 α inside 'late'
    """
    # Floors by asset class
    M_min = 0.5 if asset_class == "tech" else 0.2  # techs vs transmission links

    # Optional one-point calibration:
    # If params has a 'calibrate' dict with (O1*, O2*, M*), solve α so G hits it.
    # Otherwise, use provided α directly.
    alpha = float(params.get("alpha", 1.0))
    if "calibrate" in params and params["calibrate"] is not None:
        cal = params["calibrate"]
        O1c = cal.get("O1", np.nan)
        O2c = cal.get("O2", np.nan)
        Mc  = cal.get("M", None)
        if Mc is not None:
            # target G* from M*
            G_target = (Mc - M_min) / (1.0 - M_min)
            G_target = float(np.clip(G_target, 1e-6, 1 - 1e-6))

            # Solve alpha in closed-form (early) or 1-D find (late)
            if curve == "early":
                p   = float(params.get("p", 0.6))
                q   = float(params.get("q", 0.6))
                lam = float(params.get("lambda", params.get("lam", 1.5)))
                # f* = O1^p + lam O2^q
                o1 = 0.0 if np.isnan(O1c) else np.clip(O1c / 100.0, 0.0, 1.0) ** p
                o2 = 0.0 if np.isnan(O2c) else np.clip(O2c / 100.0, 0.0, 1.0) ** q
                fstar = o1 + lam * o2
                fstar = max(fstar, 1e-9)
                alpha = -np.log(G_target) / fstar
            elif curve == "late":
                r = float(params.get("beta_over_alpha", 1.5))
                # 1-D solve for alpha: G_target = (1 - exp(-(α W1 + β W2)))/(1 - exp(-(α+β)))
                w1 = 1.0 - (0.0 if np.isnan(O1c) else np.clip(O1c / 100.0, 0.0, 1.0))
                w2 = 1.0 - (0.0 if np.isnan(O2c) else np.clip(O2c / 100.0, 0.0, 1.0))

                def G_of_alpha(a):
                    b = r * a
                    num   = 1.0 - np.exp(-(a * w1 + b * w2))
                    denom = 1.0 - np.exp(-(a + b))
                    if denom <= 1e-12:
                        return num
                    return num / denom

                # Simple bracketed search (monotone in α)
                lo, hi = 1e-6, 50.0
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if G_of_alpha(mid) < G_target:
                        lo = mid
                    else:
                        hi = mid
                alpha = 0.5 * (lo + hi)

    # Compute G and wrap to M
    if curve == "early":
        p   = float(params.get("p", 0.6))
        q   = float(params.get("q", 0.6))
        lam = float(params.get("lambda", params.get("lam", 1.5)))
        G   = _early_G(O1, O2, alpha=alpha, p=p, q=q, lam=lam)
        return _wrap_multiplier(G, M_min)

    elif curve == "late":
        r   = float(params.get("beta_over_alpha", 1.5))
        G   = _late_G(O1, O2, alpha=alpha, beta=r * alpha)
        return _wrap_multiplier(G, M_min)

    else:
        raise ValueError(f"Unknown curve type: {curve!r}")

# ──────────────────────────────────────────────────────────────────────────────
# 3) APPLY FOR TECHS
#    scenarios is a dict: {scenario_name: {curve: 'early'|'late', alpha:..., p:..., ...}}
# ──────────────────────────────────────────────────────────────────────────────
penalty_values_techs = {scen: [] for scen in scenarios}
for _, row in df_techs.iterrows():
    for scen, params in scenarios.items():
        M = compute_multiplier(
            O1=row["O_1_filled"],
            O2=row["O_2_filled"],
            asset_class="tech",
            curve=params.get("curve", "early"),
            params=params
        )
        penalty_values_techs[scen].append(M)

df_final_techs = pd.concat(
    [df_techs[["techs", "nodes"]], pd.DataFrame(penalty_values_techs)],
    axis=1
)
df_final_techs.to_csv(output_techs, index=False)
print(f"✅ Tech multipliers saved to {output_techs}")

# ──────────────────────────────────────────────────────────────────────────────
# 4) APPLY FOR LINKS (AC OHL)
#    For links, we typically only have O2; O1 is NaN → handled inside function.
# ──────────────────────────────────────────────────────────────────────────────
penalty_values_links = {scen: [] for scen in scenarios}
for _, row in df_links.iterrows():
    for scen, params in scenarios.items():
        M = compute_multiplier(
            O1=row["O_1"],              # likely NaN for links
            O2=row["O_2"],
            asset_class="link",
            curve=params.get("curve", "late"),  # default to 'late'
            params=params
        )
        penalty_values_links[scen].append(M)

df_final_links = pd.concat(
    [df_links[["techs"]], pd.DataFrame(penalty_values_links)],
    axis=1
)
df_final_links.to_csv(output_links, index=False)
print(f"✅ Link multipliers saved to {output_links}")
