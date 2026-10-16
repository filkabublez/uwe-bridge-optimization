"""
optimise.py
Discrete Cross-Sectional Sizing Optimizer for BS EN 10210 S355J2H Hollow Sections.
Iterates through published catalog profiles sorted by mass per meter, updating structural
self-weight on every iteration and checking ULS Full-Span, ULS Half-Span, SLS Deflection,
and BS EN 1993-1-8 Chapter 7 welded K-joint capacity.
"""
from geometry import build_footbridge_geometry
from loads import assemble_load_vector
from solver import solve_3d_frame
from checks_ec3 import check_member_ec3
from joints import check_shs_k_joint
from carbon import compute_bridge_mass_and_carbon

def evaluate_design(section_map, catalog, config, include_mid_vertical=True, include_plan_bracing=True):
    geom_cfg = config["geometry"]
    nodes, members = build_footbridge_geometry(
        span_m=geom_cfg["span_m"], width_cc_m=geom_cfg["width_cc_m"],
        height_m=geom_cfg["height_m"], top_chord_length_m=geom_cfg["top_chord_length_m"],
        section_map=section_map, include_mid_vertical=include_mid_vertical,
        include_plan_bracing=include_plan_bracing
    )

    # Solve across all 3 load cases (ULS Full, ULS Half-Span Asymmetric, SLS Characteristic)
    F_uls_full = assemble_load_vector(nodes, members, catalog, config, "ULS_Full_Span")
    F_uls_half = assemble_load_vector(nodes, members, catalog, config, "ULS_Half_Span")
    F_sls = assemble_load_vector(nodes, members, catalog, config, "SLS_Char")

    sol_full = solve_3d_frame(nodes, members, catalog, F_uls_full)
    sol_half = solve_3d_frame(nodes, members, catalog, F_uls_half)
    sol_sls = solve_3d_frame(nodes, members, catalog, F_sls)

    # Evaluate governing EC3 member check across both ULS load cases
    ec3_schedule = []
    for idx, m in enumerate(members):
        sec = catalog[m["section"]]
        chk_full = check_member_ec3(sol_full["member_forces"][idx], sec, config)
        chk_half = check_member_ec3(sol_half["member_forces"][idx], sec, config)
        gov_chk = chk_full if chk_full["ur_gov"] >= chk_half["ur_gov"] else chk_half
        ec3_schedule.append(gov_chk)

    # Check welded SHS K-joints (BS EN 1993-1-8 Ch. 7) between web braces and chords
    max_diag_N = max(abs(r["N_Ed_kN"]) * 1e3 for r in ec3_schedule if r["group"] == "Web_Diagonal")
    joint_top = check_shs_k_joint(catalog[section_map["Top_Chord"]], catalog[section_map["Web_Diagonal"]], max_diag_N)
    joint_bot = check_shs_k_joint(catalog[section_map["Bot_Chord"]], catalog[section_map["Web_Diagonal"]], max_diag_N)

    carbon_data = compute_bridge_mass_and_carbon(members, catalog)
    defl_limit_mm = (geom_cfg["span_m"] * 1000.0) / 300.0

    return {
        "nodes": nodes,
        "members": members,
        "section_map": section_map.copy(),
        "schedule": ec3_schedule,
        "sol_uls": sol_full,
        "sol_sls": sol_sls,
        "max_member_ur": max(r["ur_gov"] for r in ec3_schedule),
        "max_joint_ur": max(joint_top["ur_joint"], joint_bot["ur_joint"]),
        "joint_checks": (joint_top, joint_bot),
        "sls_defl_mm": sol_sls["max_vertical_defl_mm"],
        "sls_defl_limit_mm": defl_limit_mm,
        "carbon": carbon_data
    }

def run_discrete_section_optimiser(catalog, config, target_ur=0.85):
    """
    Genuine discrete sizing loop: iterates over structural groups and selects the lightest
    published BS EN 10210 section satisfying ULS EC3 checks, SLS deflection, and EN 1993-1-8
    welded joint limits (beta >= 0.50 for practical fillet welding at fan nodes).
    """
    shs_candidates = sorted([s for s in catalog.values() if s["family"] == "SHS"], key=lambda x: x["mass_per_m"])
    rhs_candidates = sorted([s for s in catalog.values() if s["family"] == "RHS"], key=lambda x: x["mass_per_m"])

    current_map = config["baseline_sections"].copy()

    # 3 fixed-point passes to capture stiffness and self-weight redistribution (Fixes B3)
    for _ in range(3):
        for grp in ["Top_Chord", "Bot_Chord", "End_Raker", "Web_Vertical", "Web_Diagonal", "Transom", "Plan_Brace"]:
            pool = rhs_candidates if grp in ("Transom", "Plan_Brace") else shs_candidates
            for cand in pool:
                test_map = current_map.copy()
                test_map[grp] = cand["name"]

                # Enforce practical welded hollow-section width ratio beta = b_brace / b_chord >= 0.60
                b_chord_min = min(catalog[test_map["Top_Chord"]]["b_mm"], catalog[test_map["Bot_Chord"]]["b_mm"])
                b_web_max = max(catalog[test_map["Web_Vertical"]]["b_mm"], catalog[test_map["Web_Diagonal"]]["b_mm"])
                if b_web_max > b_chord_min or (catalog[test_map["Web_Diagonal"]]["b_mm"] / b_chord_min) < 0.60:
                    continue
                if catalog[test_map["End_Raker"]]["b_mm"] < b_chord_min:
                    continue

                ev = evaluate_design(test_map, catalog, config)
                grp_ur = max(r["ur_gov"] for r in ev["schedule"] if r["group"] == grp)
                if grp_ur <= target_ur and ev["max_joint_ur"] <= target_ur and ev["sls_defl_mm"] <= ev["sls_defl_limit_mm"]:
                    current_map[grp] = cand["name"]
                    break

    return evaluate_design(current_map, catalog, config)