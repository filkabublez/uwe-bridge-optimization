"""
loads.py
Derives dead, pedestrian live, service concentrated, parapet, and wind actions with
Eurocode clause references:
  - FRP grating deck spans longitudinally between floor transoms (UWE Report Sec 4.4,
    spacing s = 1.0 m for interior transoms; s_end = 0.5 m for support transoms).
  - Transom UDL w_t (N/m) is applied over clear deck width b_deck = 1.30 m via equivalent
    nodal forces (w_t * b_deck / 2) and fixed-end moments (+/- w_t * b_deck^2 / 12), and
    stored in member_udl_z for exact mid-span parabolic bending moment recovery in solver.py.
  - Crowd mass (50 kg/m^2 * 10%) is excluded from static ULS/SLS loads (only used in dynamics.py).
"""
import numpy as np


def _cfg_float(dct, *names, default=0.0):
    for key in names:
        if key in dct and dct[key] is not None:
            return float(dct[key])
    return float(default)


def assemble_load_vector(nodes, members, catalog, config, load_case="ULS_Full_Span", return_member_udls=False):
    """
    Returns global 6-DOF nodal load vector F (and optionally member_udl_z dict) for:
      - 'ULS_Full_Span'          : 1.35 * Dead + 1.35 * Live (pedestrian UDL over full span)
      - 'ULS_Half_Span'          : 1.35 * Dead + 1.35 * Live (asymmetric half-span x <= L/2)
      - 'ULS_Concentrated_Load'  : 1.35 * Dead + 1.35 * Q_fwk (10 kN concentrated service load)
      - 'ULS_Parapet'            : 1.35 * Dead + 1.35 * Live + 1.35 * Parapet horizontal load (0.8 kN/m)
      - 'ULS_Wind'               : 1.35 * Dead + 1.50 * Wind transverse action on exposed bridge depth
      - 'SLS_Char'               : 1.00 * Dead + 1.00 * Live (characteristic deflection check)
    """
    n_nodes = len(nodes)
    F = np.zeros(6 * n_nodes, dtype=float)
    member_udl_z = {}
    g = 9.81

    geom_cfg = config["geometry"]
    loads_cfg = config.get("loads", {})
    span = float(geom_cfg["span_m"])
    b_deck = float(geom_cfg["clear_deck_width_m"])
    s_transom_nom = float(geom_cfg.get("transom_tributary_spacing_m", 1.0))
    d_exposed = float(geom_cfg.get("exposed_wind_depth_m", 0.30))

    if load_case.startswith("ULS"):
        gamma_G = float(loads_cfg.get("gamma_G_uls", 1.35))
        gamma_Q_traffic = float(loads_cfg.get("gamma_Q_traffic_uls", loads_cfg.get("gamma_Q_uls", 1.35)))
        gamma_Q_wind = float(loads_cfg.get("gamma_Q_wind_uls", 1.50))
    else:
        gamma_G = 1.0
        gamma_Q_traffic = 1.0
        gamma_Q_wind = 1.0

    # 1. Structural steel member self-weight (lumped at member ends + local transom self-weight UDL)
    for m in members:
        sec = catalog[m["section"]]
        L = m["length_m"]
        w_sw_N_m = gamma_G * sec["mass_per_m"] * g
        w_member_N = w_sw_N_m * L
        F[6 * m["n1"] + 2] -= 0.5 * w_member_N
        F[6 * m["n2"] + 2] -= 0.5 * w_member_N
        if m["group"] == "Transom":
            member_udl_z[m["id"]] = w_sw_N_m

    # 2. Transom Load Path: FRP Deck Dead Load + Pedestrian Crowd Live Load (BS EN 1991-2 LM4)
    # Deck spans longitudinally between transoms (Report Sec 4.4: s = 1.0 m for interior transoms)
    deck_mass_kg_m2 = _cfg_float(loads_cfg, "frp_deck_mass_kg_m2", default=19.5)
    q_deck_N_m2 = deck_mass_kg_m2 * g  # 191.3 N/m^2
    q_live_kPa = _cfg_float(loads_cfg, "pedestrian_udl_kPa", "pedestrian_live_kPa", default=5.0)
    q_live_N_m2 = q_live_kPa * 1e3      # 5000.0 N/m^2

    include_crowd_udl = load_case in {"ULS_Full_Span", "ULS_Half_Span", "ULS_Parapet", "SLS_Char"}

    for m in members:
        if m["group"] != "Transom":
            continue
        x_t = nodes[m["n1"]][0]
        is_end_transom = (abs(x_t) < 1e-6) or (abs(x_t - span) < 1e-6)
        # Interior transoms carry s = 1.0 m tributary spacing (Report Sec 4.4); end transoms carry 0.5 * s
        s_trib = 0.5 * s_transom_nom if is_end_transom else s_transom_nom

        if load_case == "ULS_Half_Span":
            live_factor = 1.0 if x_t < (span / 2.0 - 1e-6) else (0.5 if abs(x_t - span / 2.0) < 1e-6 else 0.0)
        elif include_crowd_udl:
            live_factor = 1.0
        else:
            live_factor = 0.0

        # Distributed load w_t (N/m) along the transom over clear deck width b_deck = 1.30 m
        w_deck_and_live = (gamma_G * q_deck_N_m2 + gamma_Q_traffic * q_live_N_m2 * live_factor) * s_trib
        member_udl_z[m["id"]] = member_udl_z.get(m["id"], 0.0) + w_deck_and_live

        # Equivalent nodal vertical forces (w * b_deck / 2) & fixed-end moments (+/- w * b_deck^2 / 12)
        V_end = 0.5 * w_deck_and_live * b_deck
        M_fem = (w_deck_and_live * (b_deck ** 2)) / 12.0
        F[6 * m["n1"] + 2] -= V_end
        F[6 * m["n2"] + 2] -= V_end
        F[6 * m["n1"] + 3] -= M_fem
        F[6 * m["n2"] + 3] += M_fem

    # 3. Concentrated Load Case Q_fwk = 10 kN (BS EN 1991-2 Cl. 5.3.2.2)
    if load_case in {"ULS_Concentrated_Load", "ULS_Transom_Concentrated"}:
        Qfwk_kN = _cfg_float(loads_cfg, "concentrated_load_Qfwk_kN", default=10.0)
        qfwk_N = gamma_Q_traffic * Qfwk_kN * 1e3
        for n_id in (2, 7):
            F[6 * n_id + 2] -= 0.5 * qfwk_N

    # 4. Horizontal Parapet Load Case (BS EN 1991-2 Cl. 4.8: variable traffic action -> gamma_Q = 1.35)
    # Applied along the entire top-chord line (nodes 10..12 on Left, 13..15 on Right) by tributary length
    if load_case in {"ULS_Parapet", "ULS_Parapet_Load"}:
        q_parapet_N_m = gamma_Q_traffic * _cfg_float(loads_cfg, "parapet_horizontal_kN_m", default=0.8) * 1e3
        for top_offset in [10, 13]:
            x_top = [nodes[top_offset + i][0] for i in range(3)]
            for i in range(3):
                n_id = top_offset + i
                L_left = (x_top[i] - x_top[i - 1]) / 2.0 if i > 0 else 0.0
                L_right = (x_top[i + 1] - x_top[i]) / 2.0 if i < 2 else 0.0
                F[6 * n_id + 1] += 0.5 * q_parapet_N_m * (L_left + L_right)

    # 5. Transverse Wind Load Case (BS EN 1991-1-4: pressure * Cf * exposed_depth * tributary_length)
    # Spread along all bottom-chord nodes (0..4 and 5..9)
    if load_case in {"ULS_Wind", "ULS_Transverse_Wind"}:
        wind_kPa = _cfg_float(loads_cfg, "wind_pressure_kPa", default=0.65)
        wind_cf = _cfg_float(loads_cfg, "wind_force_coeff_Cf", default=1.6)
        q_wind_line_N_m = gamma_Q_wind * (wind_kPa * 1e3) * wind_cf * d_exposed
        for side_offset in [0, 5]:
            x_bot = [nodes[side_offset + i][0] for i in range(5)]
            for i in range(5):
                n_id = side_offset + i
                L_left = (x_bot[i] - x_bot[i - 1]) / 2.0 if i > 0 else 0.0
                L_right = (x_bot[i + 1] - x_bot[i]) / 2.0 if i < 4 else 0.0
                F[6 * n_id + 1] += 0.5 * q_wind_line_N_m * (L_left + L_right)

    if return_member_udls:
        return F, member_udl_z
    return F
