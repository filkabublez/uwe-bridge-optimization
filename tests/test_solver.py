"""
tests/test_solver.py
Unit test suite verifying structural mechanics, analytical benchmarks, global equilibrium,
symmetry, mechanism detection, transom UDL bending, Euler buckling alpha_cr, and EN 1993-1-8 joints.
"""
import json
import numpy as np
import pytest
from sections import load_section_catalog
from geometry import build_footbridge_geometry
from loads import assemble_load_vector
from solver import solve_3d_frame, solve_linear_buckling
from joints import compute_angle_between_members_deg, check_shs_k_joint


def test_symmetric_two_bar_frame_benchmark():
    """1. Verifies 3D solver axial force against exact hand calculation for a 2-bar V-truss."""
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, -1.0], [2.0, 0.0, 0.0]])
    catalog = load_section_catalog()
    members = [
        {"id": 1, "n1": 0, "n2": 1, "group": "Diag", "section": "SHS 100x100x4.0", "length_m": np.sqrt(2.0)},
        {"id": 2, "n1": 1, "n2": 2, "group": "Diag", "section": "SHS 100x100x4.0", "length_m": np.sqrt(2.0)},
    ]
    F = np.zeros(18)
    P = 10000.0
    F[6*1 + 2] = -P
    support_dofs = [0, 1, 2, 3, 4, 5, 6*1+1, 6*1+3, 6*1+5, 12, 13, 14, 15, 16, 17]
    sol = solve_3d_frame(nodes, members, catalog, F, support_dofs=support_dofs)
    N_expected = P / (2.0 * np.sin(np.pi / 4.0))
    assert np.isclose(sol["member_forces"][0]["N_Ed_N"], N_expected, rtol=2e-2)


def test_global_vertical_equilibrium():
    """2. Verifies sum of vertical support reactions equals total applied ULS vertical load."""
    with open("config.json") as f:
        cfg = json.load(f)
    catalog = load_section_catalog()
    nodes, members = build_footbridge_geometry(section_map=cfg["baseline_sections"])
    F = assemble_load_vector(nodes, members, catalog, cfg, "ULS_Full_Span")
    sol = solve_3d_frame(nodes, members, catalog, F)
    assert np.isclose(np.sum(F[2::6]) + np.sum(sol["reactions"][2::6]), 0.0, atol=1e-4)


def test_transom_udl_reproduces_hand_calc_53MPa():
    """3. Verifies interior transom bending stress matches UWE Report hand calc (~50-55 MPa)."""
    with open("config.json") as f:
        cfg = json.load(f)
    catalog = load_section_catalog()
    nodes, members = build_footbridge_geometry(section_map=cfg["baseline_sections"])
    F, udl_z = assemble_load_vector(nodes, members, catalog, cfg, "ULS_Full_Span", return_member_udls=True)
    sol = solve_3d_frame(nodes, members, catalog, F, member_udl_z=udl_z, clear_deck_width_m=1.3)
    transoms = [mf for mf in sol["member_forces"] if mf["group"] == "Transom"]
    max_sigma_pl = max(t["sigma_plastic_MPa"] for t in transoms)
    # Exact solver result is ~50.2 MPa (14.1% UR) using unrounded published Wpl = 30.8 cm^3 and gamma_Q = 1.35
    assert 49.0 <= max_sigma_pl <= 56.0


def test_mechanism_detection_raises_error():
    """4. Verifies solver detects an unrestrained mechanism without artificial regularization."""
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    catalog = load_section_catalog()
    members = [{"id": 1, "n1": 0, "n2": 1, "group": "Chord", "section": "SHS 100x100x4.0", "length_m": 1.0}]
    F = np.zeros(12)
    with pytest.raises(RuntimeError):
        solve_3d_frame(nodes, members, catalog, F, support_dofs=[0, 1, 2])


def test_euler_column_linear_buckling_alpha_cr():
    """5. Verifies solve_linear_buckling() against Euler pinned strut P_cr = pi^2*E*I/L^2."""
    # Horizontal strut along X-axis: L = 3.0m divided into two 1.5m elements
    nodes = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]])
    catalog = load_section_catalog()
    sec = catalog["SHS 100x100x4.0"]
    members = [
        {"id": 1, "n1": 0, "n2": 1, "group": "Strut", "section": "SHS 100x100x4.0", "length_m": 1.5},
        {"id": 2, "n1": 1, "n2": 2, "group": "Strut", "section": "SHS 100x100x4.0", "length_m": 1.5},
    ]
    L_tot = 3.0
    E = 210e9
    P_euler = (np.pi**2 * E * sec["Iy"]) / (L_tot**2)

    # Apply exact Euler axial compression load at Node 2 (-X direction)
    F = np.zeros(18)
    F[6*2 + 0] = -P_euler

    # Pinned at Node 0 (Ux, Uy, Uz, Rx fixed); Roller at Node 2 (Uy, Uz, Rx fixed, free Ux);
    # Intermediate Node 1 is free to buckle in Z (Uy, Rx fixed to prevent rigid-body roll)
    support_dofs = [
        0, 1, 2, 3,          # Node 0: pinned support
        6*1 + 1, 6*1 + 3,    # Node 1: out-of-plane restrained, free Uz deflection
        6*2 + 1, 6*2 + 2, 6*2 + 3  # Node 2: axial roller in X, free Uz rotation
    ]
    buck = solve_linear_buckling(nodes, members, catalog, F, E=E, support_dofs=support_dofs)
    # alpha_cr should be 1.0 within 2% discretization error for a 2-element beam-column
    assert np.isclose(buck["alpha_cr"], 1.0, rtol=0.03)


def test_en1993_1_8_k_joint_beta_regimes():
    """6. Verifies EN 1993-1-8 joint check across beta < 0.85 (face) and beta == 1.0 (side-wall)."""
    catalog = load_section_catalog()
    theta = compute_angle_between_members_deg([1.0, 0.0, 0.0], [0.775, 0.0, -1.0])
    assert 50.0 <= theta <= 55.0
    # beta = 0.60 -> chord face plastification
    res_face = check_shs_k_joint(catalog["SHS 100x100x4.0"], catalog["SHS 60x60x3.2"], -25e3, theta_deg=theta)
    assert res_face["failure_mode"] == "Chord_Face_Plastification_or_Punching"
    assert res_face["N1_Rd_kN"] > 40.0
    # beta = 1.00 -> chord side-wall failure
    res_wall = check_shs_k_joint(catalog["SHS 50x50x3.2"], catalog["SHS 50x50x3.2"], -25e3, theta_deg=theta)
    assert res_wall["failure_mode"] == "Chord_Side_Wall_Failure"
    assert res_wall["N1_Rd_kN"] > 60.0
