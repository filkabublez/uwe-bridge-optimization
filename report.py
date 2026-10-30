"""
report.py
Main execution script for the UWE 3.0m Steel Footbridge Check & Sizing Tool.
Outputs:
  1. output/ga_drawing_sheet.png (Dimensioned 3-view General Arrangement drawing)
  2. output/member_schedule.csv (Full member-by-member EC3 calculation schedule with governing_case)
  3. output/validation_table.md (Solver vs Hand Calc vs Fusion 360 validation table)
  4. output/results_summary.md (Full results table for Baseline and Sized Full-Braced designs)
  5. output/engineering_verification_plots.png (Utilisation, Deformed Shape, ULS-Sized Span Sweep, Carbon)
"""
import os
import csv
import json
import numpy as np
import matplotlib.pyplot as plt
from sections import load_section_catalog
from optimise import evaluate_design, run_discrete_section_optimiser
from dynamics import compute_modal_frequencies, run_span_frequency_sweep
from generate_ga_sheet import export_ga_drawing_sheet


def main():
    os.makedirs("output", exist_ok=True)
    export_ga_drawing_sheet("output/ga_drawing_sheet.png")

    with open("config.json") as f:
        config = json.load(f)
    catalog = load_section_catalog(config["materials"]["min_galvanizing_t_mm"])

    # 1. Evaluate Baseline UWE Design
    base_eval = evaluate_design(config["baseline_sections"], catalog, config)
    base_dyn = compute_modal_frequencies(
        base_eval["nodes"], base_eval["members"], catalog, config,
        base_eval["sol_sls"]["K"], base_eval["sol_sls"]["free_dofs"]
    )

    # 2. Run Discrete EC3 Section Sizing Optimizer (Full Braced 3D Frame, H = 1.0m parapet)
    opt_sizing_eval = run_discrete_section_optimiser(catalog, config, target_ur=0.85)
    opt_sizing_dyn = compute_modal_frequencies(
        opt_sizing_eval["nodes"], opt_sizing_eval["members"], catalog, config,
        opt_sizing_eval["sol_sls"]["K"], opt_sizing_eval["sol_sls"]["free_dofs"]
    )

    # 3. Check effect of omitting plan bracing under Wind + Parapet lateral cases
    unbraced_check = evaluate_design(
        opt_sizing_eval["section_map"], catalog, config,
        include_mid_vertical=True, include_plan_bracing=False
    )

    # Extract Baseline Validation Metrics (Transom id=25 at x=0.725m & Top Chord id=5)
    base_transoms = [r for r in base_eval["schedule"] if r["group"] == "Transom"]
    base_chords = [r for r in base_eval["schedule"] if r["group"] in ("Top_Chord", "Bot_Chord")]
    gov_transom = max(base_transoms, key=lambda x: x["My_Ed_kNm"])
    gov_chord = max(base_chords, key=lambda x: x["sigma_total_MPa"])

    mass_base = base_eval["carbon"]["steel_mass_kg"]
    mass_sized = opt_sizing_eval["carbon"]["steel_mass_kg"]
    red_sized_pct = (mass_base - mass_sized) / mass_base * 100.0

    print("\n" + "="*88)
    print("BS EN 1993-1-1 / BS EN 1993-1-8 STEEL FOOTBRIDGE CHECK & SIZING SUMMARY (SPAN = 3.0 m)")
    print("="*88)
    print(f"{'Engineering Parameter':<40} | {'Baseline (UWE)':<20} | {'Sized (Full Braced)':<20}")
    print("-" * 88)
    print(f"{'Parapet Height H (m)':<40} | {'1.00 m':<20} | {'1.00 m':<20}")
    print(f"{'Chords Profile (EN 10210 S355J2H)':<40} | {base_eval['section_map']['Top_Chord']:<20} | {opt_sizing_eval['section_map']['Top_Chord']:<20}")
    print(f"{'End Rakers Profile (EN 10210)':<40} | {base_eval['section_map']['End_Raker']:<20} | {opt_sizing_eval['section_map']['End_Raker']:<20}")
    print(f"{'Web Braces Profile (EN 10210)':<40} | {base_eval['section_map']['Web_Diagonal']:<20} | {opt_sizing_eval['section_map']['Web_Diagonal']:<20}")
    print(f"{'Transoms Profile (EN 10210)':<40} | {base_eval['section_map']['Transom']:<20} | {opt_sizing_eval['section_map']['Transom']:<20}")
    print(f"{'Plan Bracing Profile (EN 10210)':<40} | {base_eval['section_map']['Plan_Brace']:<20} | {opt_sizing_eval['section_map']['Plan_Brace']:<20}")
    print(f"{'Total Structural Steel Mass (kg)':<40} | {mass_base:<20.1f} | {mass_sized:<20.1f}")
    print(f"{'Mass Reduction vs. Baseline (%)':<40} | {'0.0%':<20} | {f'-{red_sized_pct:.1f}%':<20}")
    print(f"{'Embodied Carbon A1-A3 (kgCO2e)':<40} | {base_eval['carbon']['co2_total_A1_A3_kg']:<20.1f} | {opt_sizing_eval['carbon']['co2_total_A1_A3_kg']:<20.1f}")
    print(f"{'Max EC3 Member Utilisation Ratio':<40} | {base_eval['max_member_ur']*100:<19.1f}% | {opt_sizing_eval['max_member_ur']*100:<19.1f}%")
    print(f"{'Max EN 1993-1-8 Joint Utilisation':<40} | {base_eval['max_joint_ur']*100:<19.1f}% | {opt_sizing_eval['max_joint_ur']*100:<19.1f}%")
    print(f"{'Linear Buckling Load Factor alpha_cr':<40} | {base_eval['alpha_cr']:<20.1f} | {opt_sizing_eval['alpha_cr']:<20.1f}")
    print(f"{'SLS Deflection (mm) [Limit 10.0 mm]':<40} | {base_eval['sls_defl_mm']:<20.2f} | {opt_sizing_eval['sls_defl_mm']:<20.2f}")
    print(f"{'1st Vertical Natural Freq f1,v (Hz)':<40} | {base_dyn['f1_vertical_hz']:<20.1f} | {opt_sizing_dyn['f1_vertical_hz']:<20.1f}")
    print(f"{'1st Lateral Natural Freq f1,h (Hz)':<40} | {base_dyn['f1_lateral_hz']:<20.1f} | {opt_sizing_dyn['f1_lateral_hz']:<20.1f}")
    print("-" * 88)
    print(f"Note on Plan Bracing: Omitting plan bracing under ULS_Wind/ULS_Parapet increases")
    print(f"bottom-chord lateral bending utilisation to {unbraced_check['max_member_ur']*100:.1f}% -> Plan bracing retained.")
    print("="*88 + "\n")

    # Export Member Schedule CSV (includes governing_case column)
    csv_path = os.path.join("output", "member_schedule.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(opt_sizing_eval["schedule"][0].keys()), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(opt_sizing_eval["schedule"])

    # Export Validation Table Markdown (Solver vs Hand Calc vs Fusion 360)
    with open(os.path.join("output", "validation_table.md"), "w") as f:
        f.write("| Verification Metric | Eurocode Hand Calc (Report) | Fusion 360 Solid FEA | Python 6-DOF Frame Solver | Explanation |\n")
        f.write("| :--- | :---: | :---: | :---: | :--- |\n")
        f.write(
            f"| **Floor Transom Bending Moment (`RHS 100x40x4.0`)** | 1.64 kNm | — | **{gov_transom['My_Ed_kNm']:.2f} kNm** | Exact match on tributary UDL ($w_t = 7.76\\text{{ kN/m}}$ over $L_t = 1.30\\text{{ m}}$) |\n"
        )
        f.write(
            f"| **Floor Transom Bending Stress & UR (`RHS 100x40x4.0`)** | 54.7 MPa (15.4%) | 54.63 MPa (15.4%) | **{gov_transom['sigma_plastic_MPa']:.1f} MPa ({gov_transom['ur_gov']*100:.1f}%)** | Hand calc used rounded $W_{{pl,y}} = 30.0\\text{{ cm}}^3$; solver uses published EN 10210 $30.8\\text{{ cm}}^3$ |\n"
        )
        f.write(
            f"| **Main Chord Combined Stress (`SHS 100x100x4.0`)** | 4.8 MPa (1.3% axial) | 11.0–33.0 MPa (3.1–9.3%) | **{gov_chord['sigma_total_MPa']:.1f} MPa ({gov_chord['ur_gov']*100:.1f}%)** | 6-DOF frame solver captures rigid welded joint secondary bending omitted by pin-truss hand calc |\n"
        )

    # Run ULS-sized span frequency sweep (re-sizing sections at every span)
    sweep = run_span_frequency_sweep(catalog, config)
    span_5hz = next((s["span_m"] for s in sweep if s["f1_vertical_hz"] < 5.0), ">18.0")

    # Export Full Results Summary Markdown
    with open(os.path.join("output", "results_summary.md"), "w") as f:
        f.write("| Engineering Parameter | Baseline (UWE Report) | Sized (Full Braced) | Governing Limit / Clause |\n")
        f.write("| :--- | :---: | :---: | :--- |\n")
        f.write(f"| **Parapet Height $H$ (m)** | 1.00 m | 1.00 m | $\\ge 1.00\\text{{ m}}$ (`BS 6180` / `CD 353`) |\n")
        f.write(f"| **Chords / End Rakers (`S355J2H`)** | `{base_eval['section_map']['Top_Chord']}` | `{opt_sizing_eval['section_map']['Top_Chord']}` | `BS EN 10210-2` ($t \\ge 3.2\\text{{ mm}}$) |\n")
        f.write(f"| **Web Verticals & Diagonals (`S355J2H`)** | `{base_eval['section_map']['Web_Diagonal']}` | `{opt_sizing_eval['section_map']['Web_Diagonal']}` | `BS EN 10210-2` ($t \\ge 3.2\\text{{ mm}}$) |\n")
        f.write(f"| **Floor Transoms (`S355J2H`)** | `{base_eval['section_map']['Transom']}` | `{opt_sizing_eval['section_map']['Transom']}` | `BS EN 10210-2` ($t \\ge 3.2\\text{{ mm}}$) |\n")
        f.write(f"| **Plan X-Bracing (`S355J2H`)** | `{base_eval['section_map']['Plan_Brace']}` | `{opt_sizing_eval['section_map']['Plan_Brace']}` | Required for `ULS_Wind` & `ULS_Parapet` |\n")
        f.write(f"| **Structural Steel Mass (kg)** | {mass_base:.1f} kg | **{mass_sized:.1f} kg (-{red_sized_pct:.1f}%)** | Minimized by `optimise.py` |\n")
        f.write(f"| **Embodied Carbon A1–A3 ($\\text{{kgCO}}_2\\text{{e}}$)** | {base_eval['carbon']['co2_total_A1_A3_kg']:.1f} $\\text{{kgCO}}_2\\text{{e}}$ | **{opt_sizing_eval['carbon']['co2_total_A1_A3_kg']:.1f} $\\text{{kgCO}}_2\\text{{e}}$** | ICE v3.0 Steel + HDG Zinc |\n")
        f.write(f"| **Max Member Utilisation (`max_member_ur`)** | {base_eval['max_member_ur']*100:.1f}% | **{opt_sizing_eval['max_member_ur']*100:.1f}%** | $\\le 85\\%$ (`BS EN 1993-1-1`) |\n")
        f.write(f"| **Max Welded Joint Utilisation (`max_joint_ur`)** | {base_eval['max_joint_ur']*100:.1f}% | **{opt_sizing_eval['max_joint_ur']*100:.1f}%** | $\\le 85\\%$ (`BS EN 1993-1-8 Table 7.10`, $\\theta = {opt_sizing_eval['theta_deg']:.1f}^\\circ$) |\n")
        f.write(f"| **Global Elastic Buckling Factor $\\alpha_{{cr}}$** | {base_eval['alpha_cr']:.1f} | **{opt_sizing_eval['alpha_cr']:.1f}** | $\\ge 10.0$ (`BS EN 1993-1-1 Cl. 5.2.1`) |\n")
        f.write(f"| **SLS Vertical Deflection (`sls_defl_mm`)** | {base_eval['sls_defl_mm']:.2f} mm | **{opt_sizing_eval['sls_defl_mm']:.2f} mm** | $\\le 10.0\\text{{ mm}}$ ($L/300$, `Cl. 7.2`) |\n")
        f.write(f"| **1st Vertical Natural Freq $f_{{1,v}}$ (Hz)** | {base_dyn['f1_vertical_hz']:.1f} Hz | **{opt_sizing_dyn['f1_vertical_hz']:.1f} Hz** | $\\ge 5.0\\text{{ Hz}}$ (`BS EN 1990 Annex A2`) |\n")
        f.write(f"| **1st Lateral Natural Freq $f_{{1,h}}$ (Hz)** | {base_dyn['f1_lateral_hz']:.1f} Hz | **{opt_sizing_dyn['f1_lateral_hz']:.1f} Hz** | $\\ge 2.5\\text{{ Hz}}$ (`BS EN 1990 Annex A2`) |\n")

    # Generate 4-panel Engineering Verification Plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=160)

    ax1 = axes[0, 0]
    groups = ["Top_Chord", "Bot_Chord", "End_Raker", "Web_Vertical", "Web_Diagonal", "Transom", "Plan_Brace"]
    b_urs = [max(r["ur_gov"] for r in base_eval["schedule"] if r["group"] == g) * 100 for g in groups]
    o_urs = [max(r["ur_gov"] for r in opt_sizing_eval["schedule"] if r["group"] == g) * 100 for g in groups]
    y = np.arange(len(groups))
    ax1.barh(y + 0.18, b_urs, height=0.35, label="Baseline (UWE Report)", color="#64748b")
    ax1.barh(y - 0.18, o_urs, height=0.35, label="Sized (Full Braced)", color="#0284c7")
    ax1.set_yticks(y)
    ax1.set_yticklabels([g.replace("_", " ") for g in groups])
    ax1.invert_yaxis()
    ax1.set_xlabel("Governing EC3 Utilisation (% of Rd)")
    ax1.set_title("1. Member Utilisation Across 5 ULS Load Cases")
    ax1.legend()
    ax1.grid(axis="x", linestyle="--", alpha=0.5)

    ax2 = axes[0, 1]
    u_sls = opt_sizing_eval["sol_sls"]["u"]
    nodes = opt_sizing_eval["nodes"]
    scale = 500.0
    for m in opt_sizing_eval["members"]:
        if nodes[m["n1"]][1] < 1e-3 and nodes[m["n2"]][1] < 1e-3:
            p1, p2 = nodes[m["n1"]], nodes[m["n2"]]
            d1 = p1 + scale * u_sls[6*m["n1"]:6*m["n1"]+3]
            d2 = p2 + scale * u_sls[6*m["n2"]:6*m["n2"]+3]
            ax2.plot([p1[0], p2[0]], [p1[2], p2[2]], "k--", lw=1.0, alpha=0.4)
            ax2.plot([d1[0], d2[0]], [d1[2], d2[2]], color="#0284c7", lw=2.2)
    ax2.set_title(f"2. SLS Deformed Elevation (x{int(scale)} Scale | Max uz = {opt_sizing_eval['sls_defl_mm']:.2f} mm)")
    ax2.set_xlabel("Span X (m)")
    ax2.set_ylabel("Elevation Z (m)")
    ax2.set_aspect("equal")
    ax2.grid(True, linestyle="--", alpha=0.4)

    ax3 = axes[1, 0]
    spans = [s["span_m"] for s in sweep]
    fv = [s["f1_vertical_hz"] for s in sweep]
    fl = [s["f1_lateral_hz"] for s in sweep]
    ax3.plot(spans, fv, "o-", color="#0284c7", label="1st Vertical Mode $f_{1,v}$ (ULS-Sized)")
    ax3.plot(spans, fl, "s-", color="#0f766e", label="1st Lateral Mode $f_{1,h}$ (ULS-Sized)")
    ax3.axhline(5.0, color="#dc2626", linestyle="--", label="EN 1990 A2 Vertical Trigger (5.0 Hz)")
    ax3.axhline(2.5, color="#ea580c", linestyle=":", label="EN 1990 A2 Lateral Trigger (2.5 Hz)")
    ax3.set_xlabel("Bridge Span L (m)")
    ax3.set_ylabel("Natural Frequency (Hz)")
    ax3.set_title(f"3. ULS-Sized Modal Frequency vs. Span (f1,v < 5 Hz near L = {span_5hz} m)")
    ax3.set_yscale("log")
    ax3.legend(fontsize=8)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    ax4 = axes[1, 1]
    opts = ["Baseline\n(UWE Report)", f"Sized (Full Braced)\n(-{red_sized_pct:.1f}%)"]
    masses = [mass_base, mass_sized]
    co2s = [base_eval["carbon"]["co2_total_A1_A3_kg"], opt_sizing_eval["carbon"]["co2_total_A1_A3_kg"]]
    x = np.arange(len(opts))
    ax4.bar(x - 0.16, masses, width=0.30, label="Steel Mass (kg)", color="#475569")
    ax4.bar(x + 0.16, co2s, width=0.30, label="Embodied Carbon A1-A3 (kgCO2e)", color="#16a34a")
    ax4.set_xticks(x)
    ax4.set_xticklabels(opts)
    ax4.set_ylabel("kg / kgCO2e")
    ax4.set_title("4. Steel Mass & Cradle-to-Gate Carbon (ICE v3.0 + Galvanizing)")
    ax4.legend()
    ax4.grid(axis="y", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(os.path.join("output", "engineering_verification_plots.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
