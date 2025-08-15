import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
import yaml
from mpl_toolkits.mplot3d import Axes3D

# === 1. Omega model ===
# The cost multiplier ω is computed as a function of general opposition (O1) and local opposition (O2).
# The function starts at ω = 0 when both opposition terms are at their worst (1.0).
# It rises towards ω = k (usually 1.0) as opposition decreases.
def omega_model(O1, O2, params):
    k, alpha, beta = params  # k: max ω value; alpha: weight of O1; beta: weight of O2
    # ADDITIVE exponential decay (saturates as O1 and O2 increase)
    return k * (1 - np.exp(-(alpha * (1 - O1) + beta * (1 - O2))))
    # ω = k * [1 - e^(-α(1 - O1) - β(1 - O2))]

# === 2. Input data ===
# Your sample data of opposition (O1, O2) and what ω you'd like to get from the function
# O1 and O2 must be paired element-wise (zip) — so you are fitting on the diagonal only
O1_levels = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
O2_levels = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])

# These are the target ω values you'd like your model to replicate
expected_behaviour = {
    "low":    np.array([1.00, 0.94, 0.88, 0.83, 0.78, 0.74, 0.71, 0.70]),
    "medium": np.array([1.00, 0.88, 0.78, 0.70, 0.65, 0.62, 0.60, 0.60]),
    "high":   np.array([1.00, 0.85, 0.75, 0.65, 0.58, 0.54, 0.52, 0.50])
}

# === 3. Fitting function ===
# The optimisation objective: fit model parameters (k, alpha, beta)
# so that predicted ω matches the expected target values
def objective(params, O1, O2, target):
    preds = [omega_model(o1, o2, params) for o1, o2 in zip(O1, O2)]  # zip assumes element-wise matching!
    weights = np.linspace(2, 0.5, len(target))  # weight early points more
    return np.sum(weights * (np.array(preds) - target)**2)  # weighted squared error

fitted_params = {}
for scenario, targets in expected_behaviour.items():
    res = minimize(
        objective,
        x0=[1.0, 0.3, 0.3],  # initial guess
        args=(O1_levels, O2_levels, targets),
        bounds=[(0.1, 4), (0.01, 2), (0.01, 2)]  # k fixed at 1.0, tuning α and β
    )
    fitted_params[scenario] = res.x

# === 4. Plot 3D omega surface ===
def plot_surface(params, title):
    k, alpha, beta = params
    O1_grid, O2_grid = np.meshgrid(np.linspace(0, 0.6, 50), np.linspace(0, 0.6, 50))
    omega_vals = omega_model(O1_grid, O2_grid, params)

    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(O1_grid*100, O2_grid*100, omega_vals, cmap='viridis')
    ax.set_xlabel("O₁: General Opposition (%)")
    ax.set_ylabel("O₂: Local Opposition (%)")
    ax.set_zlabel("ω (Cost Multiplier)")
    ax.set_title(title)
    plt.tight_layout()
    plt.show()

for scenario, params in fitted_params.items():
    plot_surface(params, f"{scenario.capitalize()} Impact")

# === 5. Save the fitted parameters ===
with open("fitted_omega_params.yaml", "w") as f:
    yaml.dump({
        s: {"k": float(p[0]), "alpha": float(p[1]), "beta": float(p[2])}
        for s, p in fitted_params.items()
    }, f)

# Save to CSV
pd.DataFrame([
    {"scenario": s, "k": p[0], "alpha": p[1], "beta": p[2]}
    for s, p in fitted_params.items()
]).to_csv("fitted_omega_params.csv", index=False)
