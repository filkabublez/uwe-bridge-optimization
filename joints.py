"""
joints.py
Simplified, indicative hollow-section nodal joint check based on BS EN 1993-1-8
geometric validity limits (0.35 <= beta <= 1.0, b0/t0 <= 35, theta >= 30 deg)
and approximate chord face plastification / side-wall bearing expressions.
"""
import numpy as np


def compute_angle_between_members_deg(v1, v2):
    """Computes acute angle in degrees between two 3D member direction vectors."""
    u1 = np.asarray(v1, dtype=float) / np.linalg.norm(v1)
    u2 = np.asarray(v2, dtype=float) / np.linalg.norm(v2)
    cos_theta = np.clip(abs(np.dot(u1, u2)), 0.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def check_shs_k_joint(chord_sec, brace_sec, N_brace_Ed_N, theta_deg=None, fy0=355e6, gamma_M5=1.25, E=210e9):
    if theta_deg is None:
        raise ValueError("theta_deg must be computed from actual member geometry, not defaulted.")

    b0, h0, t0 = chord_sec["b"], chord_sec["h"], chord_sec["t"]
    b1, h1 = brace_sec["b"], brace_sec["h"]
    beta = b1 / b0
    theta = np.radians(theta_deg)
    sin_t = np.sin(theta)

    b0_over_t0 = b0 / t0
    valid_geometry = (0.35 <= beta <= 1.0 + 1e-6) and (b0_over_t0 <= 35.0) and (theta_deg >= 30.0)

    gamma_half = b0 / (2.0 * t0)
    k_n = 1.0

    # Simplified indicative face plastification and punching shear check (beta <= 0.85)
    def _face_and_punching_rd(beta_val, b1_val):
        N_plast = (k_n * fy0 * (t0**2) / sin_t) * (
            (2.0 * b1_val / b0) + 4.0 * np.sqrt(max(1e-6, 1.0 - beta_val))
        ) * (gamma_half ** 0.5) / gamma_M5
        b_ep = min(b1_val, (10.0 / b0_over_t0) * b1_val)
        N_punch = (fy0 / np.sqrt(3.0)) * t0 * (
            (2.0 * b1_val / sin_t) + 2.0 * b_ep
        ) / gamma_M5
        return min(N_plast, N_punch)

    # Chord side-wall bearing/buckling resistance at beta = 1.0
    def _sidewall_rd(h1_val):
        lambda_w = ((h0 / t0) - 2.0) / sin_t
        lambda_1 = np.pi * np.sqrt(E / fy0)
        lambda_bar_w = 0.546 * (lambda_w / lambda_1)
        phi_w = 0.5 * (1.0 + 0.21 * (lambda_bar_w - 0.2) + lambda_bar_w**2)
        chi_w = min(1.0, 1.0 / (phi_w + np.sqrt(max(1e-9, phi_w**2 - lambda_bar_w**2))))
        f_b = chi_w * fy0 if N_brace_Ed_N < 0 else fy0
        return (2.0 * f_b * t0 / sin_t) * ((h1_val / sin_t) + 10.0 * t0) / gamma_M5

    if beta <= 0.85:
        N1_Rd = _face_and_punching_rd(beta, b1)
        failure_mode = "Chord_Face_Plastification_or_Punching"
    elif beta >= 0.999:
        N1_Rd = _sidewall_rd(h1)
        failure_mode = "Chord_Side_Wall_Failure"
    else:
        rd_085 = _face_and_punching_rd(0.85, 0.85 * b0)
        rd_100 = _sidewall_rd(h1)
        w = (beta - 0.85) / 0.15
        N1_Rd = (1.0 - w) * rd_085 + w * rd_100
        failure_mode = "Interpolated_Face_to_SideWall"

    ur_joint = abs(N_brace_Ed_N) / N1_Rd if valid_geometry else 9.99
    return {
        "chord": chord_sec["name"],
        "brace": brace_sec["name"],
        "beta": round(beta, 3),
        "theta_deg": round(theta_deg, 1),
        "valid_geometry": valid_geometry,
        "failure_mode": failure_mode,
        "N1_Ed_kN": abs(N_brace_Ed_N) / 1e3,
        "N1_Rd_kN": N1_Rd / 1e3,
        "ur_joint": ur_joint
    }


def check_shs_joint_en1993_1_8(chord_sec, brace_sec, N_brace_Ed_N, theta_deg, fy0=355e6, gamma_M5=1.25):
    return check_shs_k_joint(chord_sec, brace_sec, N_brace_Ed_N, theta_deg=theta_deg, fy0=fy0, gamma_M5=gamma_M5)
