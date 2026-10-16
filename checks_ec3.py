"""
checks_ec3.py
BS EN 1993-1-1:2005 Member Resistance Checks for Hot-Finished S355J2H Hollow Sections:
  - Cross-section Classification (Table 5.2): c/t <= 33*epsilon for Class 1
  - Tension Resistance (Cl. 6.2.3): Nt,Rd = A * fy / gamma_M0
  - Compression Resistance (Cl. 6.2.4): Nc,Rd = A * fy / gamma_M0
  - Bending Resistance (Cl. 6.2.5): Mc,Rd = Wpl * fy / gamma_M0 (Class 1/2)
  - Shear Resistance (Cl. 6.2.6): Vpl,Rd = Av * (fy / sqrt(3)) / gamma_M0
  - Flexural Buckling (Cl. 6.3.1 & Table 6.2): Curve 'a' (alpha = 0.21) for EN 10210 hot-finished SHS/RHS
  - System Buckling Lengths (Annex BB.1.3): 0.90*L for chords, 0.75*L for web braces
  - Combined Bending & Axial Compression (Cl. 6.2.9 / Cl. 6.3.3)
"""
import numpy as np

def classify_shs_section(sec, fy=355e6):
    """Classifies compression flange/web slenderness c/t per BS EN 1993-1-1 Table 5.2."""
    epsilon = np.sqrt(235e6 / fy)  # 0.814 for S355
    c_mm = max(sec["h_mm"], sec["b_mm"]) - 3.0 * sec["t_mm"]
    c_over_t = c_mm / sec["t_mm"]
    if c_over_t <= 33.0 * epsilon:
        return 1, c_over_t
    elif c_over_t <= 38.0 * epsilon:
        return 2, c_over_t
    elif c_over_t <= 42.0 * epsilon:
        return 3, c_over_t
    return 4, c_over_t

def check_member_ec3(mf, sec, config):
    E = config["materials"]["E_Pa"]
    fy = config["materials"]["fy_Pa"]
    gamma_M0 = config["materials"]["gamma_M0"]
    gamma_M1 = config["materials"]["gamma_M1"]

    sec_class, c_t = classify_shs_section(sec, fy)
    W_y = sec["Wpl_y"] if sec_class in (1, 2) else sec["Wel_y"]
    W_z = sec["Wpl_z"] if sec_class in (1, 2) else sec["Wel_z"]

    N_Ed = mf["N_Ed_N"]
    V_Ed = mf["V_Ed_N"]
    My_Ed = mf["My_Ed_Nm"]
    Mz_Ed = mf["Mz_Ed_Nm"]
    L = mf["length_m"]
    grp = mf["group"]

    # 1. Cross-section axial, bending & shear resistances (Cl. 6.2)
    N_Rd = (sec["A"] * fy) / gamma_M0
    My_Rd = (W_y * fy) / gamma_M0
    Mz_Rd = (W_z * fy) / gamma_M0
    V_Rd = (sec["Av"] * (fy / np.sqrt(3.0))) / gamma_M0

    ur_axial = abs(N_Ed) / N_Rd
    ur_bending = (My_Ed / My_Rd) + (Mz_Ed / Mz_Rd)
    ur_shear = V_Ed / V_Rd
    ur_cross_section = ur_axial + ur_bending

    # 2. Flexural Buckling Check (Cl. 6.3.1, Buckling Curve 'a' -> alpha = 0.21 for EN 10210)
    # Effective length factors per BS EN 1993-1-1 Annex BB.1.3 & U-frame action
    if grp in ("Web_Diagonal", "Web_Vertical", "Plan_Brace"):
        k_eff = 0.75
    elif grp == "Top_Chord":
        k_eff = 0.95  # Half-through chord stabilized by U-frame verticals/transoms
    else:
        k_eff = 0.90

    L_cr = k_eff * L
    I_min = min(sec["Iy"], sec["Iz"])
    N_cr = (np.pi**2 * E * I_min) / (L_cr**2)
    lambda_bar = float(np.sqrt((sec["A"] * fy) / N_cr))

    alpha_imp = 0.21
    phi = 0.5 * (1.0 + alpha_imp * (lambda_bar - 0.2) + lambda_bar**2)
    chi = float(min(1.0, 1.0 / (phi + np.sqrt(max(1e-9, phi**2 - lambda_bar**2)))))
    Nb_Rd = (chi * sec["A"] * fy) / gamma_M1

    if N_Ed < 0:  # Compression member
        ur_buckling = (abs(N_Ed) / Nb_Rd) + (My_Ed / My_Rd) + (Mz_Ed / Mz_Rd)
    else:
        ur_buckling = 0.0

    checks = {
        "Cross_Section_M_N": ur_cross_section,
        "Flexural_Buckling": ur_buckling,
        "Shear": ur_shear
    }
    gov_name = max(checks, key=checks.get)
    return {
        "id": mf["id"],
        "group": grp,
        "section": sec["name"],
        "class": sec_class,
        "length_m": L,
        "N_Ed_kN": N_Ed / 1e3,
        "N_Rd_kN": (Nb_Rd if N_Ed < 0 else N_Rd) / 1e3,
        "My_Ed_kNm": My_Ed / 1e3,
        "My_Rd_kNm": My_Rd / 1e3,
        "V_Ed_kN": V_Ed / 1e3,
        "sigma_total_MPa": mf["sigma_total_MPa"],
        "lambda_bar": lambda_bar,
        "chi": chi,
        "ur_gov": checks[gov_name],
        "governing_check": gov_name
    }