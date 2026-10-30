"""
optimise.py
Discrete Cross-Sectional Sizing Optimizer for BS EN 10210 S355J2H Hollow Sections.
Evaluates all five ULS load cases (Full-Span, Half-Span, Concentrated Q_fwk, Parapet,
and Transverse Wind) plus SLS Deflection, Linear Elastic Buckling (alpha_cr), and
BS EN 1993-1-8 Chapter 7 welded joint resistance with exact geometric angle theta.
"""
from geometry import build_footbridge_geometry
from loads import assemble_load_vector
from solver import solve_3d_frame, solve_linear_buckling
from checks_ec3 import check_member_ec3
from joints import compute_angle_between_members_deg, check_shs_k_joint
from carbon import compute_bridge_mass_and_carbon

ULS_LOAD_CASES = [
    "ULS_Full_Span",
    "ULS_Half_Span",
    "ULS_Concentrated_Load",
    "ULS_Parapet",
    "ULS_Wind"
]


def evaluate_design(section_map, catalog, config, include_mid_vertical=True, include_plan_bracing=True):
    geom_cfg = config["geometry"]
    mat_cfg = config["materials"]
    b_deck = float(geom_cfg["clear_deck_width_m"])

    nodes, members = build_footbridge_geometry(
        span_m=geom_cfg["span_m"], width_cc_m=geom_cfg["width_cc_m"],
        height_m=geom_cfg["height_m"], top_chord_length_m=geom_cfg["top_chord_length_m"],
        section_map=section_map, include_mid_vertical=include_mid_vertical,
        include_plan_bracing=include_plan_bracing
    )

    # Solve all 5 ULS load cases
    uls_solutions = {}
    for lc in ULS_LOAD_CASES:
        F_lc, udl_lc = assemble_load_vector(nodes, members, catalog, config, lc, return_member_udls=True)
        uls_solutions[lc] = solve_3d_frame(
            nodes, members, catalog, F_lc, E=mat_cfg["E_Pa"], G=mat_cfg["G_Pa"],
            member_udl_z=udl_lc, clear_deck_width_m=b_deck
        )

    # Solve SLS characteristic deflection case
    F_sls, udl_sls = assemble_load_vector(nodes, members, catalog, config, "SLS_Char", return_member_udls=True)
    sol_sls = solve_3d_frame(
        nodes, members, catalog, F_sls, E=mat_cfg["E_Pa"], G=mat_cfg["G_Pa"],
        member_udl_z=udl_sls, clear_deck_width_m=b_deck
    )

    # Compute Linear Elastic Buckling factor alpha_cr under ULS_Full_Span
    F_full, udl_full = assemble_load_vector(nodes, members, catalog, config, "ULS_Full_Span", return_member_udls=True)
    buck_res = solve_linear_buckling(
        nodes, members, catalog, F_full, E=mat_cfg["E_Pa"], G=mat_cfg["G_Pa"],
        member_udl_z=udl_full, clear_deck_width_m=b_deck
    )

    # Evaluate governing EC3 member check across all 5 ULS load cases
    ec3_schedule = []
    for idx, m in enumerate(members):
        sec = catalog[m["section"]]
        gov_chk = None
        gov_case = None
        for lc in ULS_LOAD_CASES:
            chk = check_member_ec3(uls_solutions[lc]["member_forces"][idx], sec, config)
            if gov_chk is None or chk["ur_gov"] > gov_chk["ur_gov"]:
                gov_chk = chk
                gov_case = lc
        gov_chk["governing_case"] = gov_case
        ec3_schedule.append(gov_chk)

    # Compute true geometric brace-to-chord angle theta from 3D node coordinates
    # Top chord node 10 -> 11 vs Web diagonal node 10 -> 2
    v_chord = nodes[11] - nodes[10]
    v_diag = nodes[2] - nodes[10]
    theta_deg = compute_angle_between_members_deg(v_chord, v_diag)

    # Worst-case axial force in web diagonals across all ULS cases
    diag_indices = [i for i, m in enumerate(members) if m["group"] == "Web_Diagonal"]
    max_diag_N = 0.0
    for lc in ULS_LOAD_CASES:
        for d_idx in diag_indices:
            N_val = uls_solutions[lc]["member_forces"][d_idx]["N_Ed_N"]
            if abs(N_val) > abs(max_diag_N):
                max_diag_N = N_val

    joint_top = check_shs_k_joint(
        catalog[section_map["Top_Chord"]], catalog[section_map["Web_Diagonal"]],
        max_diag_N, theta_deg=theta_deg, fy0=mat_cfg["fy_Pa"], gamma_M5=mat_cfg["gamma_M5_joints"]
    )
    joint_bot = check_shs_k_joint(
        catalog[section_map["Bot_Chord"]], catalog[section_map["Web_Diagonal"]],
        max_diag_N, theta_deg=theta_deg, fy0=mat_cfg["fy_Pa"], gamma_M5=mat_cfg["gamma_M5_joints"]
    )

    carbon_data = compute_bridge_mass_and_carbon(members, catalog)
    defl_limit_mm = (geom_cfg["span_m"] * 1000.0) / float(config["loads"].get("project_sls_defl_ratio", 300.0))

    return {
        "nodes": nodes,
        "members": members,
        "section_map": section_map.copy(),
        "schedule": ec3_schedule,
        "sol_uls": uls_solutions["ULS_Full_Span"],
        "uls_solutions": uls_solutions,
        "sol_sls": sol_sls,
        "alpha_cr": buck_res["alpha_cr"],
        "theta_deg": theta_deg,
        "max_member_ur": max(r["ur_gov"] for r in ec3_schedule),
        "max_joint_ur": max(joint_top["ur_joint"], joint_bot["ur_joint"]),
        "joint_checks": (joint_top, joint_bot),
        "sls_defl_mm": sol_sls["max_vertical_defl_mm"],
        "sls_defl_limit_mm": defl_limit_mm,
        "carbon": carbon_data
    }


def run_discrete_section_optimiser(catalog, config, target_ur=0.85):
    """
    Discrete sizing loop: iterates over structural groups and selects the lightest
    published BS EN 10210 section satisfying all 5 ULS load cases, SLS deflection,
    global elastic buckling (alpha_cr >= 10), and BS EN 1993-1-8 welded joint checks.
    """
    shs_candidates = sorted([s for s in catalog.values() if s["family"] == "SHS"], key=lambda x: x["mass_per_m"])
    rhs_candidates = sorted([s for s in catalog.values() if s["family"] == "RHS"], key=lambda x: x["mass_per_m"])

    current_map = config["baseline_sections"].copy()

    for _ in range(3):
        for grp in ["Top_Chord", "Bot_Chord", "End_Raker", "Web_Vertical", "Web_Diagonal", "Transom", "Plan_Brace"]:
            pool = rhs_candidates if grp in ("Transom", "Plan_Brace") else shs_candidates
            for cand in pool:
                test_map = current_map.copy()
                test_map[grp] = cand["name"]

                b_chord_min = min(catalog[test_map["Top_Chord"]]["b_mm"], catalog[test_map["Bot_Chord"]]["b_mm"])
                b_web_max = max(catalog[test_map["Web_Vertical"]]["b_mm"], catalog[test_map["Web_Diagonal"]]["b_mm"])
                if b_web_max > b_chord_min or (catalog[test_map["Web_Diagonal"]]["b_mm"] / b_chord_min) < 0.60:
                    continue
                if catalog[test_map["End_Raker"]]["b_mm"] < b_chord_min:
                    continue

                ev = evaluate_design(test_map, catalog, config)
                grp_ur = max(r["ur_gov"] for r in ev["schedule"] if r["group"] == grp)
                if (
                    grp_ur <= target_ur
                    and ev["max_joint_ur"] <= target_ur
                    and ev["sls_defl_mm"] <= ev["sls_defl_limit_mm"]
                    and ev["alpha_cr"] >= 10.0
                ):
                    current_map[grp] = cand["name"]
                    break

    return evaluate_design(current_map, catalog, config)
