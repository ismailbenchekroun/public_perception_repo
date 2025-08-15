import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# === Step 1: Define Opposition Levels (in fraction, not %) ===
O_levels = np.array([0.0, 0.05, 0.2, 0.4, 0.6])

# === Step 2: Define Your Expected ω for Each Scenario ===
# You can change these numbers to describe what you expect
expected_behaviour = {
    "low":    np.array([1.0, 0.9, 0.8, 0.75, 0.7]),
    "medium": np.array([1.0, 0.8, 0.6, 0.5, 0.4]),
    "high":   np.array([1.0, 0.7, 0.5, 0.4, 0.3])
}

# === Step 3: Define the Nested Exponential Function ===
def omega_model(O, params):
    k, alpha, beta = params
    return k * (1 - np.exp(-alpha * (1 - O)) * np.exp(-beta * (1 - O)))

# === Step 4: Calibration Function ===
def objective(params, O, target):
    pred = np.array([omega_model(o, params) for o in O])
    return np.sum((pred - target)**2)

# === Step 5: Fit Parameters for Each Scenario ===
fitted = {}
for scenario, targets in expected_behaviour.items():
    init_guess = [1.0, 0.3, 0.3]
    bounds = [(1, 1.0), (0.05, 2), (0.05, 2)]  # k fixed ≤1 to ensure penalty
    res = minimize(objective, init_guess, args=(O_levels, targets), bounds=bounds)
    fitted[scenario] = res.x

# === Step 6: Plot Curves ===
O_plot = np.linspace(0, 0.6, 200)
plt.figure(figsize=(8,6))

for scenario, params in fitted.items():
    plt.plot(O_plot*100, [omega_model(o, params) for o in O_plot], label=f"{scenario} (k={params[0]:.2f}, α={params[1]:.2f}, β={params[2]:.2f})")
    plt.scatter(O_levels*100, expected_behaviour[scenario], s=50)

plt.title("Calibrated Cost Multiplier ω vs Public Opposition")
plt.xlabel("Public Opposition O (%)")
plt.ylabel("Cost Multiplier ω")
plt.legend()
plt.grid(True)
plt.show()

print("=== Fitted Parameters ===")
for s, p in fitted.items():
    print(f"{s.capitalize()} Impact: k={p[0]:.3f}, α={p[1]:.3f}, β={p[2]:.3f}")
