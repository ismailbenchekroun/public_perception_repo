import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title("ω Calibration: Public Acceptance Cost Multiplier")
st.markdown("""
This tool lets you interactively shape your public acceptance penalty scalar ω
based on opposition levels O₁ (general) and O₂ (local). 
""")

# Sliders for parameters
scenario = st.selectbox("Scenario", ["low", "medium", "high"])
k = st.slider("k", 1.0, 1.2, 1.0, 0.01)
alpha = st.slider("α (O₁ sensitivity)", 0.01, 2.0, 0.3, 0.01)
beta = st.slider("β (O₂ sensitivity)", 0.01, 2.0, 0.3, 0.01)

# Meshgrid to simulate ω over full range of O₁ and O₂
O1_grid, O2_grid = np.meshgrid(np.linspace(0, 0.6, 50), np.linspace(0, 0.6, 50))
omega_vals = k * (1 - np.exp(-alpha * (1 - O1_grid)) * np.exp(-beta * (1 - O2_grid)))

# 3D plot using matplotlib
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(O1_grid*100, O2_grid*100, omega_vals, cmap='plasma')
ax.set_xlabel("O₁: General Opposition (%)")
ax.set_ylabel("O₂: Local Opposition (%)")
ax.set_zlabel("ω")
ax.set_title(f"ω Surface - {scenario.capitalize()} Scenario")
st.pyplot(fig)
