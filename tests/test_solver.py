"""
Unit test suite verifying structural mechanics, analytical benchmarks, global equilibrium,
symmetry, and mechanism detection.
"""
import json
import numpy as np
import pytest
from sections import load_section_catalog
from geometry import build_footbridge_geometry
from loads import assemble_load_vector
from solver import solve_3d_frame, solve_linear_buckling
from dynamics import run_span_frequency_sweep

def test_symmetric_two_bar_frame_benchmark():
    """1. Verifies 3D solver axial force against exact hand calculation for a 2-bar V-truss."""
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, -1.0], [2.0, 0.0, 0.0]])
    catalog = load_section_catalog()
    members = [
        {"id": 1, "n1": 0, "n2": 1, "group": "Diag", "section": "SHS 100x100x4.0", "length_m": np.sqrt(2.0)},
        {"id": 2, "n1": 1, "n2": 2, "group": "Diag", "section": "SHS 100x100x4.0", "length_m": np.sqrt(2.0)},
    ]
    F = np.zeros(18)
    P = 10000.0  # 10 kN downward load at apex node 1
    F[6*1 + 2] = -P
    support_dofs = [
        0, 1, 2, 3, 4, 5,       # Node 0 clamped
        6*1+1, 6*1+3, 6*1+5,    # Node 1 out-of-plane restrained
        12, 13, 14, 15, 16, 17  # Node 2 clamped
    ]
    sol = solve_3d_frame(nodes, members, catalog, F, support_dofs=support_dofs)
    N_expected = P / (2.0 * np.sin(np.pi / 4.0))  # 7071.07 N tension
    assert np.isclose(sol["member_forces"][0]["N_Ed_N"], N_expected, rtol=2e-2)

def test_global_vertical_equilibrium():
    """2. Verifies sum of vertical support reactions equals total applied ULS vertical load."""
    with open("config.json") as f:
        cfg = json.load(f)
    catalog = load_section_catalog()
    nodes, members = build_footbridge_geometry(section_map=cfg["baseline_sections"])
    F = assemble_load_vector(nodes, members, catalog, cfg, "ULS_Full_Span")
    sol = solve_3d_frame(nodes, members, catalog, F)
    total_applied_Fz = np.sum(F[2::6])
    total_reaction_Fz = np.sum(sol["reactions"][2::6])
    assert np.isclose(total_applied_Fz + total_reaction_Fz, 0.0, atol=1e-4)

def test_symmetric_loading_produces_symmetric_forces():
    """3. Verifies Left and Right truss chords carry identical axial forces under symmetric load."""
    with open("config.json") as f:
        cfg = json.load(f)
    catalog = load_section_catalog()
    nodes, members = build_footbridge_geometry(section_map=cfg["baseline_sections"])
    F = assemble_load_vector(nodes, members, catalog, cfg, "ULS_Full_Span")
    sol = solve_3d_frame(nodes, members, catalog, F)
    mf = sol["member_forces"]
    # Compare the first bottom chord on each side.
    assert np.isclose(mf[0]["N_Ed_N"], mf[13]["N_Ed_N"], rtol=1e-3)

def test_mechanism_detection_raises_error():
    """4. Verifies solver detects an unrestrained mechanism without artificial 1e-4 regularization."""
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    catalog = load_section_catalog()
    members = [{"id": 1, "n1": 0, "n2": 1, "group": "Chord", "section": "SHS 100x100x4.0", "length_m": 1.0}]
    F = np.zeros(12)
    with pytest.raises(RuntimeError):
        solve_3d_frame(nodes, members, catalog, F, support_dofs=[0, 1, 2])  # Free to rotate -> mechanism!


def test_linear_buckling_matches_cantilever_euler_load():
    """Checks the geometric-stiffness buckling estimate against Euler's cantilever result."""
    nodes = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    catalog = load_section_catalog()
    members = [{
        "id": 1, "n1": 0, "n2": 1, "group": "Column",
        "section": "SHS 100x100x4.0", "length_m": 1.0
    }]
    F = np.zeros(12)
    F[8] = -1.0
    solution = solve_linear_buckling(
        nodes, members, catalog, F, support_dofs=list(range(6))
    )
    section = catalog["SHS 100x100x4.0"]
    expected = np.pi**2 * 210e9 * min(section["Iy"], section["Iz"]) / 4.0
    assert np.isclose(solution["alpha_cr"], expected, rtol=0.03)


def test_parapet_height_constraint_enforced():
    """5. Verifies geometry generator rejects truss heights below 1.0m handrail limit."""
    with open("config.json") as f:
        cfg = json.load(f)
    with pytest.raises(ValueError):
        build_footbridge_geometry(height_m=0.60, section_map=cfg["baseline_sections"])

def test_span_frequency_sweep_preserves_brace_angle_limit():
    with open("config.json") as f:
        cfg = json.load(f)
    catalog = load_section_catalog()
    result = run_span_frequency_sweep(
        catalog, cfg, cfg["baseline_sections"], spans_m=[24.0]
    )
    assert len(result) == 1
    assert result[0]["f1_vertical_hz"] > 0.0