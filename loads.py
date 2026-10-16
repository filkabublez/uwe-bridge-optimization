"""
loads.py
Derives dead and pedestrian live actions with explicit Eurocode clause references:
  - BS EN 1990:2002 + UK NA (Table NA.A2.4(B)): gamma_G = 1.35, gamma_Q = 1.50
  - BS EN 1991-2:2003 Cl. 5.3.2.1 (Load Model 4): q_fk = 5.0 kN/m^2 for short-span footbridges
  - FRP Grating Deck Datasheet: 19.5 kg/m^2 over clear width b_deck = 1.30 m
  - Dynamic structural steel self-weight recalculated on every sizing iteration (Fixes B3).
"""
import numpy as np

def assemble_load_vector(nodes, members, catalog, config, load_case="ULS_Full_Span"):
    """
    Returns global 6-DOF nodal load vector F (size 6 * n_nodes) for:
      - 'ULS_Full_Span' : 1.35 * Dead + 1.50 * Live (full span)
      - 'ULS_Half_Span' : 1.35 * Dead + 1.50 * Live (x <= L/2 asymmetric shear check)
      - 'SLS_Char'      : 1.00 * Dead + 1.00 * Live (characteristic deflection check)
    """
    n_nodes = len(nodes)
    F = np.zeros(6 * n_nodes, dtype=float)
    g = 9.81
    span = config["geometry"]["span_m"]
    b_deck = config["geometry"]["clear_deck_width_m"]

    if load_case.startswith("ULS"):
        gamma_G = config["loads"]["gamma_G_uls"]
        gamma_Q = config["loads"]["gamma_Q_uls"]
    else:
        gamma_G = 1.0
        gamma_Q = 1.0

    # 1. Exact member steel self-weight (updated dynamically when sections change)
    for m in members:
        sec = catalog[m["section"]]
        L = m["length_m"]
        w_member_N = gamma_G * (sec["mass_per_m"] * L * g)
        F[6 * m["n1"] + 2] -= 0.5 * w_member_N
        F[6 * m["n2"] + 2] -= 0.5 * w_member_N

    # 2. FRP Decking Dead Load + Pedestrian Live Load (BS EN 1991-2 Cl. 5.3.2.1)
    q_deck_N_m2 = config["loads"]["frp_deck_mass_kg_m2"] * g  # 191.3 N/m^2
    q_live_N_m2 = config["loads"]["pedestrian_live_kPa"] * 1e3 # 5000.0 N/m^2

    # Distribute deck tributary lengths to bottom chord nodes (Left: 0..4, Right: 5..9)
    for side_offset in [0, 5]:
        x_coords = [nodes[side_offset + i][0] for i in range(5)]
        for i in range(5):
            n_id = side_offset + i
            L_trib_left = (x_coords[i] - x_coords[i-1]) / 2.0 if i > 0 else 0.0
            L_trib_right = (x_coords[i+1] - x_coords[i]) / 2.0 if i < 4 else 0.0
            L_trib = L_trib_left + L_trib_right
            A_trib = L_trib * (b_deck / 2.0)

            # Dead load of FRP deck
            F[6 * n_id + 2] -= gamma_G * q_deck_N_m2 * A_trib

            # Pedestrian live load (Full span vs Half-span asymmetric)
            if load_case == "ULS_Half_Span":
                live_factor = 1.0 if x_coords[i] < (span / 2.0 - 1e-6) else (0.5 if abs(x_coords[i] - span/2.0) < 1e-6 else 0.0)
            else:
                live_factor = 1.0
            F[6 * n_id + 2] -= gamma_Q * q_live_N_m2 * A_trib * live_factor

    return F