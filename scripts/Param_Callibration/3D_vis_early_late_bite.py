import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

# ─────────────────────────────
# Early-bite function
# ─────────────────────────────
def early_G(O1, O2, alpha, p=0.6, q=0.6, lam=1.5):
    o1 = np.clip(O1, 0.0, 1.0)
    o2 = np.clip(O2, 0.0, 1.0)
    val = (o1 ** p) + lam * (o2 ** q)
    return np.exp(-alpha * val)

# ─────────────────────────────
# Late-bite (Koecklin-like) function
# ─────────────────────────────
def late_G(O1, O2, alpha, beta_over_alpha=1.5):
    beta = beta_over_alpha * alpha
    o1 = np.clip(O1, 0.0, 1.0)
    o2 = np.clip(O2, 0.0, 1.0)
    w1 = 1.0 - o1
    w2 = 1.0 - o2
    num   = 1.0 - np.exp(- (alpha * w1 + beta * w2))
    denom = 1.0 - np.exp(- (alpha + beta))
    if denom <= 1e-12:
        return num
    return num / denom

# ─────────────────────────────
# Initial grid
# ─────────────────────────────
res = 50
O1_vals = np.linspace(0, 1, res)
O2_vals = np.linspace(0, 1, res)
O1_grid, O2_grid = np.meshgrid(O1_vals, O2_vals)

# ─────────────────────────────
# Initial parameters
# ─────────────────────────────
alpha0 = 0.7
p0 = 0.6
q0 = 0.6
lam0 = 1.5
beta_ratio0 = 1.5

mode = "Early"  # starting mode

# ─────────────────────────────
# Plot setup
# ─────────────────────────────
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')

def plot_surface():
    ax.clear()
    if mode == "Early":
        Z = early_G(O1_grid, O2_grid,
                    alpha=slider_alpha.val,
                    p=slider_p.val,
                    q=slider_q.val,
                    lam=slider_lam.val)
        title = f"Early-Bite Surface (α={slider_alpha.val:.3f}, p={slider_p.val:.2f}, q={slider_q.val:.2f}, λ={slider_lam.val:.2f})"
    else:
        Z = late_G(O1_grid, O2_grid,
                   alpha=slider_alpha.val,
                   beta_over_alpha=slider_beta_ratio.val)
        title = f"Late-Bite Surface (α={slider_alpha.val:.3f}, β/α={slider_beta_ratio.val:.2f})"
    ax.plot_surface(O1_grid*100, O2_grid*100, Z, cmap='viridis')
    ax.set_xlabel("O1 (%)")
    ax.set_ylabel("O2 (%)")
    ax.set_zlabel("Penalty Factor")
    ax.set_title(title)
    ax.set_zlim(0, 1)
    plt.draw()

# ─────────────────────────────
# Sliders & Controls
# ─────────────────────────────
ax_mode = plt.axes([0.02, 0.7, 0.15, 0.15])
radio_mode = RadioButtons(ax_mode, ("Early", "Late"), active=0)

ax_alpha = plt.axes([0.25, 0.10, 0.50, 0.02])
slider_alpha = Slider(ax_alpha, 'Alpha', 0.01, 5.0, valinit=alpha0, valstep=0.01)

# Early-specific sliders
ax_p   = plt.axes([0.25, 0.07, 0.50, 0.02])
ax_q   = plt.axes([0.25, 0.04, 0.50, 0.02])
ax_lam = plt.axes([0.25, 0.01, 0.50, 0.02])
slider_p   = Slider(ax_p,   'p (O1 exp)', 0.1, 2.0, valinit=p0, valstep=0.01)
slider_q   = Slider(ax_q,   'q (O2 exp)', 0.1, 2.0, valinit=q0, valstep=0.01)
slider_lam = Slider(ax_lam, 'λ (O2 weight)', 0.5, 3.0, valinit=lam0, valstep=0.01)

# Late-specific slider
ax_beta_ratio = plt.axes([0.25, 0.07, 0.50, 0.02])
slider_beta_ratio = Slider(ax_beta_ratio, 'β/α', 0.5, 3.0, valinit=beta_ratio0, valstep=0.01)
slider_beta_ratio.ax.set_visible(False)  # hidden initially for early mode

# ─────────────────────────────
# Update functions
# ─────────────────────────────
def update(val):
    plot_surface()

def switch_mode(label):
    global mode
    mode = label
    if mode == "Early":
        slider_p.ax.set_visible(True)
        slider_q.ax.set_visible(True)
        slider_lam.ax.set_visible(True)
        slider_beta_ratio.ax.set_visible(False)
    else:
        slider_p.ax.set_visible(False)
        slider_q.ax.set_visible(False)
        slider_lam.ax.set_visible(False)
        slider_beta_ratio.ax.set_visible(True)
    plot_surface()

slider_alpha.on_changed(update)
slider_p.on_changed(update)
slider_q.on_changed(update)
slider_lam.on_changed(update)
slider_beta_ratio.on_changed(update)
radio_mode.on_clicked(switch_mode)

# ─────────────────────────────
# Initial plot
# ─────────────────────────────
plot_surface()
plt.show()
