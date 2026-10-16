"""
dynamics.py
Modal eigenvalue analysis solving [K]{phi} = omega^2 [M]{phi} on the 3D frame stiffness
and lumped mass matrices, plus a parametric span sweep evaluating BS EN 1990 Annex A2
footbridge vibration comfort triggers (5.0 Hz vertical, 2.5 Hz lateral).
"""
import numpy as np
from scipy.linalg import eigh
from geometry import build_footbridge_geometry
from solver import solve_3d_frame
from loads import assemble_load_vector

def compute_modal_frequencies(nodes, members, catalog, config, K, free_dofs):
    """Assembles 3D lumped mass matrix [M] and solves generalized eigenvalue problem."""
    n_nodes = len(nodes)
    ndof = 6 * n_nodes
    M_diag = np.zeros(ndof, dtype=float)

    # 1. Structural steel member mass lumped at end nodes
    for m in members:
        sec = catalog[m["section"]]
        m_half = 0.5 * sec["mass_per_m"] * m["length_m"]
        rot_inertia = m_half * (sec["h"]**2) / 12.0
        for n_id in (m["n1"], m["n2"]):
            M_diag[6*n_id + 0] += m_half
            M_diag[6*n_id + 1] += m_half
            M_diag[6*n_id + 2] += m_half
            M_diag[6*n_id + 3] += rot_inertia
            M_diag[6*n_id + 4] += rot_inertia
            M_diag[6*n_id + 5] += rot_inertia

    # 2. FRP Decking mass + 10% pedestrian mass participation on bottom chord nodes
    span = config["geometry"]["span_m"]
    b_deck = config["geometry"]["clear_deck_width_m"]
    deck_mass_total = config["loads"]["frp_deck_mass_kg_m2"] * span * b_deck
    mass_per_bot_node = deck_mass_total / 10.0
    for n_id in range(10):
        for j in range(3):
            M_diag[6*n_id + j] += mass_per_bot_node

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    M_ff = np.diag(M_diag[free_dofs])

    eigvals, eigvecs = eigh(K_ff, M_ff)
    pos_mask = eigvals > 1e-2
    freqs_hz = np.sqrt(eigvals[pos_mask]) / (2.0 * np.pi)
    vecs = eigvecs[:, pos_mask]

    # Classify first vertical (Z-dominant) and first lateral (Y-dominant) modes
    f1_vert_hz, f1_lat_hz = None, None
    for idx, f_hz in enumerate(freqs_hz):
        mode_full = np.zeros(ndof)
        mode_full[free_dofs] = vecs[:, idx]
        uy_energy = np.sum(mode_full[1::6]**2)
        uz_energy = np.sum(mode_full[2::6]**2)
        if f1_lat_hz is None and uy_energy > uz_energy:
            f1_lat_hz = float(f_hz)
        if f1_vert_hz is None and uz_energy >= uy_energy:
            f1_vert_hz = float(f_hz)
        if f1_lat_hz is not None and f1_vert_hz is not None:
            break

    return {
        "f1_vertical_hz": f1_vert_hz or float(freqs_hz[0]),
        "f1_lateral_hz": f1_lat_hz or float(freqs_hz[0]),
        "requires_vertical_dynamic_check": (f1_vert_hz or freqs_hz[0]) < 5.0,
        "requires_lateral_dynamic_check": (f1_lat_hz or freqs_hz[0]) < 2.5
    }

def run_span_frequency_sweep(catalog, config, section_map, spans_m=None):
    """Sweeps span L from 3m to 24m to identify critical vibration threshold spans."""
    if spans_m is None:
        spans_m = np.linspace(3.0, 24.0, 15)
    sweep_results = []
    for L in spans_m:
        span = float(L)
        top_chord_length = 0.517 * span
        x_sh = (span - top_chord_length) / 2.0
        x_mid = span / 2.0
        brace_run = max(x_sh, x_mid - x_sh)
        min_height = (
            np.tan(np.deg2rad(config["geometry"]["min_brace_angle_deg"]))
            * brace_run
            * (1.0 + 1e-9)
        )
        height = max(config["geometry"]["height_m"], min_height)
        cfg = {
            **config,
            "geometry": {
                **config["geometry"],
                "span_m": span,
                "top_chord_length_m": top_chord_length,
                "height_m": height
            }
        }
        nodes, members = build_footbridge_geometry(
            span_m=span, width_cc_m=cfg["geometry"]["width_cc_m"],
            height_m=height, top_chord_length_m=top_chord_length,
            section_map=section_map
        )
        F_sls = assemble_load_vector(nodes, members, catalog, cfg, "SLS_Char")
        sol = solve_3d_frame(nodes, members, catalog, F_sls)
        dyn = compute_modal_frequencies(nodes, members, catalog, cfg, sol["K"], sol["free_dofs"])
        sweep_results.append({
            "span_m": float(L),
            "f1_vertical_hz": dyn["f1_vertical_hz"],
            "f1_lateral_hz": dyn["f1_lateral_hz"]
        })
    return sweep_results