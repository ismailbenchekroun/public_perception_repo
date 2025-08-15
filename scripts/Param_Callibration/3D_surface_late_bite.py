import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- Late-bite function (Koecklin-like, but adjustable) ---
def late_G(O1, O2, alpha, beta_over_alpha=1.5):
    """Late bite penalty factor surface generator."""
    beta = beta_over_alpha * alpha
    o1 = np.clip(O1, 0.0, 1.0)
    o2 = np.clip(O2, 0.0, 1.0)
    w1 = 1.0 - o1
    w2 = 1.0 - o2
    num   = 1.0 - np.exp(-(alpha * w1 + beta * w2))
    denom = 1.0 - np.exp(-(alpha + beta))
    return num / denom if denom > 1e-12 else num

# --- Initial grid ---
res = 50
O1_vals = np.linspace(0, 1, res)
O2_vals = np.linspace(0, 1, res)
O1_grid, O2_grid = np.meshgrid(O1_vals, O2_vals)

# --- Initial parameters ---
alpha0 = 0.7
beta_over_alpha0 = 1.5

# --- Set up figure ---
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

def plot_surface(alpha, beta_over_alpha):
    ax.clear()
    Z = late_G(O1_grid, O2_grid, alpha=alpha, beta_over_alpha=beta_over_alpha)
    ax.plot_surface(O1_grid*100, O2_grid*100, Z, cmap='viridis')
    ax.set_xlabel("O1 (%)")
    ax.set_ylabel("O2 (%)")
    ax.set_zlabel("Penalty Factor")
    ax.set_title(f"Late-Bite Surface (α={alpha:.3f}, β/α={beta_over_alpha:.2f})")
    ax.set_zlim(0, 1)
    plt.draw()

plot_surface(alpha0, beta_over_alpha0)

# --- Slider axes ---
ax_alpha = plt.axes([0.25, 0.05, 0.50, 0.02])
ax_beta_ratio = plt.axes([0.25, 0.01, 0.50, 0.02])

slider_alpha = Slider(ax_alpha, 'Alpha', 0.01, 5.0, valinit=alpha0, valstep=0.01)
slider_beta_ratio = Slider(ax_beta_ratio, 'Beta/Alpha', 0.5, 3.0, valinit=beta_over_alpha0, valstep=0.01)

def update(val):
    plot_surface(slider_alpha.val, slider_beta_ratio.val)

slider_alpha.on_changed(update)
slider_beta_ratio.on_changed(update)

plt.show()
