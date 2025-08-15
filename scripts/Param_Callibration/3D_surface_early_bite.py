import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- Early-bite function ---
def early_G(O1, O2, alpha, p=0.6, q=0.6, lam=1.5):
    """Early-bite penalty factor surface generator."""
    o1 = np.clip(O1, 0.0, 1.0)
    o2 = np.clip(O2, 0.0, 1.0)
    val = (o1 ** p) + lam * (o2 ** q)
    return np.exp(-alpha * val)

# --- Initial grid ---
res = 50
O1_vals = np.linspace(0, 1, res)
O2_vals = np.linspace(0, 1, res)
O1_grid, O2_grid = np.meshgrid(O1_vals, O2_vals)

# --- Initial parameters ---
alpha0 = 0.7
p0 = 0.6
q0 = 0.6
lam0 = 1.5

# --- Set up figure ---
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

def plot_surface(alpha, p, q, lam):
    ax.clear()
    Z = early_G(O1_grid, O2_grid, alpha=alpha, p=p, q=q, lam=lam)
    ax.plot_surface(O1_grid*100, O2_grid*100, Z, cmap='viridis')
    ax.set_xlabel("O1 (%)")
    ax.set_ylabel("O2 (%)")
    ax.set_zlabel("Penalty Factor")
    ax.set_title(f"Early-Bite Surface (α={alpha:.3f}, p={p:.2f}, q={q:.2f}, λ={lam:.2f})")
    ax.set_zlim(0, 1)
    plt.draw()

plot_surface(alpha0, p0, q0, lam0)

# --- Slider axes ---
ax_alpha = plt.axes([0.25, 0.10, 0.50, 0.02])
ax_p     = plt.axes([0.25, 0.07, 0.50, 0.02])
ax_q     = plt.axes([0.25, 0.04, 0.50, 0.02])
ax_lam   = plt.axes([0.25, 0.01, 0.50, 0.02])

slider_alpha = Slider(ax_alpha, 'Alpha', 0.01, 5.0, valinit=alpha0, valstep=0.01)
slider_p     = Slider(ax_p,     'p (O1 exp)', 0.1, 2.0, valinit=p0, valstep=0.01)
slider_q     = Slider(ax_q,     'q (O2 exp)', 0.1, 2.0, valinit=q0, valstep=0.01)
slider_lam   = Slider(ax_lam,   'λ (O2 weight)', 0.5, 3.0, valinit=lam0, valstep=0.01)

def update(val):
    plot_surface(slider_alpha.val, slider_p.val, slider_q.val, slider_lam.val)

slider_alpha.on_changed(update)
slider_p.on_changed(update)
slider_q.on_changed(update)
slider_lam.on_changed(update)

plt.show()

