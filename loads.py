"""
loads.py
Derives dead and live actions with Eurocode-compatible load factors and support for
legacy and modern config keys.
"""
import numpy as np


def _cfg_float(dct, *names, default=0.0):
    for key in names:
        if key in dct and dct[key] is not None:
            return float(dct[key])
    return float(default)


def assemble_load_vector(nodes, members, catalog, config, load_case="ULS_Full_Span"):
    """
    Returns global 6-DOF nodal load vector F (size 6 * n_nodes) for:
      - 'ULS_Full_Span' : 1.35 * Dead + 1.35 * Live (traffic) full span
      - 'ULS_Half_Span' : asymmetrical half-span traffic check
      - 'SLS_Char'      : characteristic deflection check
    """
    n_nodes = len(nodes)
    F = np.zeros(6 * n_nodes, dtype=float)
    g = 9.81
    span = config["geometry"]["span_m"]
    b_deck = config["geometry"]["clear_deck_width_m"]
    loads_cfg = config.get("loads", {})

    if load_case.startswith("ULS"):
        gamma_G = float(loads_cfg.get("gamma_G_uls", 1.35))
        gamma_Q = float(loads_cfg.get("gamma_Q_traffic_uls", loads_cfg.get("gamma_Q_uls", 1.35)))
        gamma_Q_wind = float(loads_cfg.get("gamma_Q_wind_uls", 1.50))
    else:
        gamma_G = 1.0
        gamma_Q = 1.0
        gamma_Q_wind = 1.0

    # 1. Member self-weight (legacy behaviour retained for compatibility)
    for m in members:
        sec = catalog[m["section"]]
        L = m["length_m"]
        w_member_N = gamma_G * (sec["mass_per_m"] * L * g)
        F[6 * m["n1"] + 2] -= 0.5 * w_member_N
        F[6 * m["n2"] + 2] -= 0.5 * w_member_N

    # 2. Deck load and crowd load at bottom nodes
    deck_mass_kg_m2 = _cfg_float(loads_cfg, "frp_deck_mass_kg_m2")
    q_deck_N_m2 = deck_mass_kg_m2 * g
    q_live_kPa = _cfg_float(loads_cfg, "pedestrian_udl_kPa", "pedestrian_live_kPa", default=5.0)
    q_live_N_m2 = q_live_kPa * 1e3
    q_crowd_kg_m2 = 50.0
    q_crowd_N_m2 = q_crowd_kg_m2 * g

    for side_offset in [0, 5]:
        x_coords = [nodes[side_offset + i][0] for i in range(5)]
        for i in range(5):
            n_id = side_offset + i
            L_trib_left = (x_coords[i] - x_coords[i-1]) / 2.0 if i > 0 else 0.0
            L_trib_right = (x_coords[i+1] - x_coords[i]) / 2.0 if i < 4 else 0.0
            L_trib = L_trib_left + L_trib_right
            A_trib = L_trib * (b_deck / 2.0)
            F[6 * n_id + 2] -= gamma_G * q_deck_N_m2 * A_trib
            if load_case == "ULS_Half_Span":
                live_factor = 1.0 if x_coords[i] < (span / 2.0 - 1e-6) else (0.5 if abs(x_coords[i] - span / 2.0) < 1e-6 else 0.0)
            else:
                live_factor = 1.0
            F[6 * n_id + 2] -= gamma_Q * q_live_N_m2 * A_trib * live_factor
            if load_case.startswith("ULS"):
                F[6 * n_id + 2] -= gamma_Q * q_crowd_N_m2 * A_trib * 0.10 * live_factor

    # 3. Additional governing Eurocode actions are activated only in their specific cases.
    if load_case in {"ULS_Concentrated_Load", "ULS_Transom_Concentrated"}:
        Qfwk_kN = _cfg_float(loads_cfg, "concentrated_load_Qfwk_kN", default=0.0)
        if Qfwk_kN > 0.0:
            qfwk = gamma_Q * Qfwk_kN * 1e3
            for n_id in (2, 7):
                F[6 * n_id + 2] -= 0.5 * qfwk

    if load_case in {"ULS_Parapet", "ULS_Parapet_Load"}:
        parapet_kN_m = _cfg_float(loads_cfg, "parapet_horizontal_kN_m", default=0.0)
        if parapet_kN_m > 0.0:
            parapet = gamma_Q_wind * parapet_kN_m * 1e3
            for n_id in (10, 13):
                F[6 * n_id + 1] += parapet * 0.5

    if load_case in {"ULS_Wind", "ULS_Transverse_Wind"}:
        wind_kPa = _cfg_float(loads_cfg, "wind_pressure_kPa", default=0.0)
        wind_coeff = _cfg_float(loads_cfg, "wind_force_coeff_Cf", default=1.0)
        if wind_kPa > 0.0:
            wind_force = gamma_Q_wind * wind_kPa * 1e3 * wind_coeff * b_deck
            for side_offset in [0, 5]:
                n_id = side_offset + 2
                F[6 * n_id + 1] += 0.5 * wind_force

    return F