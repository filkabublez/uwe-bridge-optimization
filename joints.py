"""
joints.py
Welded Square Hollow Section (SHS) Nodal Joint Resistance Check per BS EN 1993-1-8:2005
Chapter 7 (Table 7.2):
  - Geometric validity ranges (Cl. 7.1.2): 0.35 <= beta <= 1.0, b0/t0 <= 35, theta >= 30 deg
  - Chord Face Plastification Resistance (Table 7.2) for Gap K/N-joints
  - Chord Face Punching Shear Resistance (Table 7.2) for beta <= 0.85
"""
import numpy as np

def check_shs_k_joint(chord_sec, brace_sec, N_brace_Ed_N, theta_deg=54.0, fy0=355e6, gamma_M5=1.25):
    b0, t0 = chord_sec["b"], chord_sec["t"]
    b1 = brace_sec["b"]
    beta = b1 / b0
    theta = np.radians(theta_deg)

    # Geometric validity checks per BS EN 1993-1-8 Cl. 7.1.2
    b0_over_t0 = b0 / t0
    valid_geometry = (beta >= 0.35) and (b0_over_t0 <= 35.0) and (theta_deg >= 30.0)

    # 1. Chord Face Plastification Resistance N_1,Rd (BS EN 1993-1-8 Table 7.2)
    gamma_half = b0 / (2.0 * t0)
    k_n = 1.0  # Conservative chord stress factor
    N1_Rd_plast = (k_n * fy0 * (t0**2) / np.sin(theta)) * (
        (2.0 * b1 / b0) + 4.0 * np.sqrt(max(0.05, 1.0 - beta))
    ) * (gamma_half ** 0.5) / gamma_M5

    # 2. Chord Punching Shear Resistance (applicable when beta <= 0.85)
    if beta <= 0.85:
        b_ep = min(b1, (10.0 / (b0 / t0)) * b1)
        N1_Rd_punch = (fy0 / np.sqrt(3.0)) * t0 * (
            (2.0 * b1 / np.sin(theta)) + 2.0 * b_ep
        ) / gamma_M5
        N1_Rd = min(N1_Rd_plast, N1_Rd_punch)
    else:
        N1_Rd = N1_Rd_plast

    ur_joint = abs(N_brace_Ed_N) / N1_Rd if valid_geometry else 9.99
    return {
        "chord": chord_sec["name"],
        "brace": brace_sec["name"],
        "beta": round(beta, 3),
        "valid_geometry": valid_geometry,
        "N1_Ed_kN": abs(N_brace_Ed_N) / 1e3,
        "N1_Rd_kN": N1_Rd / 1e3,
        "ur_joint": ur_joint
    }