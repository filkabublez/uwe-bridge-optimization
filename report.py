"""
report.py
Main execution script for the UWE 3.0m Steel Footbridge Check & Sizing Tool.
Outputs:
  1. output/member_schedule.csv (Full member-by-member EC3 calculation schedule)
  2. output/engineering_verification_plots.png (Utilization, Deformed Shape, Span-Frequency Sweep, Carbon)
"""
import os
import csv
import json
import numpy as np
import matplotlib.pyplot as plt
from sections import load_section_catalog
from optimise import evaluate_design, run_discrete_section_optimiser
from dynamics import compute_modal_frequencies, run_span_frequency_sweep

def main():
    os.makedirs("output", exist_ok=True)
    with open("config.json") as f:
        config = json.load(f)
    catalog = load_section_catalog(config["materials"]["min_galvanizing_t_mm"])

    # 1. Evaluate Baseline UWE Design
    base_eval = evaluate_design(config["baseline_sections"], catalog, config)
    base_dyn = compute_modal_frequencies(
        base_eval["nodes"], base_eval["members"], catalog, config,
        base_eval["sol_uls"]["K"], base_eval["sol_uls"]["free_dofs"]
    )

    # 2. Run Discrete EC3 Section Sizing Optimizer (Keeping full 3D topology & H = 1.0m parapet)
    opt_sizing_eval = run_discrete_section_optimiser(catalog, config, target_ur=0.85)
    opt_sizing_dyn = compute_modal_frequencies(
        opt_sizing_eval["nodes"], opt_sizing_eval["members"], catalog, config,
        opt_sizing_eval["sol_uls"]["K"], opt_sizing_eval["sol_uls"]["free_dofs"]
    )

    # 3. Evaluate Rationalized Topology + Sizing (Removes redundant mid-vertical & plan bracing)
    opt_rational_eval = evaluate_design(
        opt_sizing_eval["section_map"], catalog, config,
        include_mid_vertical=False, include_plan_bracing=False
    )

    mass_base = base_eval["carbon"]["steel_mass_kg"]
    mass_sized = opt_sizing_eval["carbon"]["steel_mass_kg"]
    mass_rat = opt_rational_eval["carbon"]["steel_mass_kg"]

    red_sized_pct = (mass_base - mass_sized) / mass_base * 100.0
    red_rat_pct = (mass_base - mass_rat) / mass_base * 100.0

    print("\n" + "="*96)
    print("BS EN 1993-1-1 / BS EN 1993-1-8 STEEL FOOTBRIDGE CHECK & SIZING SUMMARY (SPAN = 3.0 m)")
    print("="*96)
    print(f"{'Engineering Parameter':<36} | {'Baseline (UWE)':<18} | {'Sized (Full Braced)':<18} | {'Sized (Rationalized)':<18}")
    print("-" * 96)
    print(f"{'Parapet Height H (m)':<36} | {'1.00 m':<18} | {'1.00 m':<18} | {'1.00 m':<18}")
    print(f"{'Chords Profile (EN 10210 S355J2H)':<36} | {base_eval['section_map']['Top_Chord']:<18} | {opt_sizing_eval['section_map']['Top_Chord']:<18} | {opt_rational_eval['section_map']['Top_Chord']:<18}")
    print(f"{'Web Braces Profile (EN 10210)':<36} | {base_eval['section_map']['Web_Diagonal']:<18} | {opt_sizing_eval['section_map']['Web_Diagonal']:<18} | {opt_rational_eval['section_map']['Web_Diagonal']:<18}")
    print(f"{'Transoms Profile (EN 10210)':<36} | {base_eval['section_map']['Transom']:<18} | {opt_sizing_eval['section_map']['Transom']:<18} | {opt_rational_eval['section_map']['Transom']:<18}")
    print(f"{'Total Structural Steel Mass (kg)':<36} | {mass_base:<18.1f} | {mass_sized:<18.1f} | {mass_rat:<18.1f}")
    print(f"{'Mass Reduction vs. Baseline (%)':<36} | {'0.0%':<18} | {f'-{red_sized_pct:.1f}%':<18} | {f'-{red_rat_pct:.1f}%':<18}")
    print(f"{'Embodied Carbon A1-A3 (kgCO2e)':<36} | {base_eval['carbon']['co2_total_A1_A3_kg']:<18.1f} | {opt_sizing_eval['carbon']['co2_total_A1_A3_kg']:<18.1f} | {opt_rational_eval['carbon']['co2_total_A1_A3_kg']:<18.1f}")
    print(f"{'Max EC3 Member Utilisation Ratio':<36} | {base_eval['max_member_ur']*100:<17.1f}% | {opt_sizing_eval['max_member_ur']*100:<17.1f}% | {opt_rational_eval['max_member_ur']*100:<17.1f}%")
    print(f"{'Max EN 1993-1-8 K-Joint Utilisation':<36} | {base_eval['max_joint_ur']*100:<17.1f}% | {opt_sizing_eval['max_joint_ur']*100:<17.1f}% | {opt_rational_eval['max_joint_ur']*100:<17.1f}%")
    print(f"{'SLS Deflection (mm) [Limit 10.0 mm]':<36} | {base_eval['sls_defl_mm']:<18.3f} | {opt_sizing_eval['sls_defl_mm']:<18.3f} | {opt_rational_eval['sls_defl_mm']:<18.3f}")
    print(f"{'1st Vertical Natural Freq f1 (Hz)':<36} | {base_dyn['f1_vertical_hz']:<18.1f} | {opt_sizing_dyn['f1_vertical_hz']:<18.1f} | {opt_sizing_dyn['f1_vertical_hz']:<18.1f}")
    print("="*96 + "\n")

    # Export Member Schedule CSV (Fixes C6)
    csv_path = os.path.join("output", "member_schedule.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(opt_sizing_eval["schedule"][0].keys()))
        writer.writeheader()
        writer.writerows(opt_sizing_eval["schedule"])

    # Generate Engineering Calculation Plots (Fixes C8 & C9)
    sweep = run_span_frequency_sweep(catalog, config, opt_sizing_eval["section_map"])
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=160)

    # Plot 1: Governing Member Utilisation by Structural Family
    ax1 = axes[0, 0]
    groups = ["Top_Chord", "Bot_Chord", "End_Raker", "Web_Vertical", "Web_Diagonal", "Transom", "Plan_Brace"]
    b_urs = [max(r["ur_gov"] for r in base_eval["schedule"] if r["group"] == g) * 100 for g in groups]
    o_urs = [max(r["ur_gov"] for r in opt_sizing_eval["schedule"] if r["group"] == g) * 100 for g in groups]
    y = np.arange(len(groups))
    ax1.barh(y + 0.18, b_urs, height=0.35, label="Baseline (100x100 SHS)", color="#64748b")
    ax1.barh(y - 0.18, o_urs, height=0.35, label="Sized (EN 10210 Min Pool)", color="#0284c7")
    ax1.set_yticks(y)
    ax1.set_yticklabels([g.replace("_", " ") for g in groups])
    ax1.invert_yaxis()
    ax1.set_xlabel("Governing EC3 Utilisation (% of Rd)")
    ax1.set_title("1. Member Utilisation by Structural Group (ULS)")
    ax1.legend()
    ax1.grid(axis="x", linestyle="--", alpha=0.5)

    # Plot 2: 2D Side-Elevation Deformed Shape (Exaggerated x500)
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

    # Plot 3: Modal Frequency vs. Span Sweep (EN 1990 Annex A2 Comfort Triggers)
    ax3 = axes[1, 0]
    spans = [s["span_m"] for s in sweep]
    fv = [s["f1_vertical_hz"] for s in sweep]
    fl = [s["f1_lateral_hz"] for s in sweep]
    ax3.plot(spans, fv, "o-", color="#0284c7", label="1st Vertical Mode $f_{1,v}$")
    ax3.plot(spans, fl, "s-", color="#0f766e", label="1st Lateral Mode $f_{1,h}$")
    ax3.axhline(5.0, color="#dc2626", linestyle="--", label="EN 1990 A2 Vertical Trigger (5.0 Hz)")
    ax3.axhline(2.5, color="#ea580c", linestyle=":", label="EN 1990 A2 Lateral Trigger (2.5 Hz)")
    ax3.set_xlabel("Bridge Span L (m)")
    ax3.set_ylabel("Natural Frequency (Hz)")
    ax3.set_title("3. Eigenvalue Frequency vs. Span (EN 1990 Annex A2)")
    ax3.set_yscale("log")
    ax3.legend(fontsize=8)
    ax3.grid(True, which="both", linestyle="--", alpha=0.4)

    # Plot 4: Steel Mass & Embodied Carbon (A1-A3) Across Options
    ax4 = axes[1, 1]
    opts = ["Baseline\n(UWE Report)", f"Sized (Full)\n(-{red_sized_pct:.1f}%)", f"Sized (Rationalized)\n(-{red_rat_pct:.1f}%)"]
    masses = [mass_base, mass_sized, mass_rat]
    co2s = [base_eval["carbon"]["co2_total_A1_A3_kg"], opt_sizing_eval["carbon"]["co2_total_A1_A3_kg"], opt_rational_eval["carbon"]["co2_total_A1_A3_kg"]]
    x = np.arange(len(opts))
    ax4.bar(x - 0.18, masses, width=0.35, label="Steel Mass (kg)", color="#475569")
    ax4.bar(x + 0.18, co2s, width=0.35, label="Embodied Carbon A1-A3 (kgCO2e)", color="#16a34a")
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