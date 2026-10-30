"""
generate_ga_sheet.py
Generates a dimensioned 3-view CAD-style General Arrangement (GA) drawing sheet
(Elevation, Plan View, and Cross-Section) of the 3.0m UWE Bristol footbridge.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def export_ga_drawing_sheet(output_path="output/ga_drawing_sheet.png"):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig = plt.figure(figsize=(14, 8), dpi=180, facecolor="#0f172a")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0], height_ratios=[1.0, 1.0], hspace=0.32, wspace=0.22)

    def _style_cad_axes(ax, title):
        ax.set_facecolor("#0f172a")
        ax.set_title(title, color="#f8fafc", fontsize=10.5, fontweight="bold", pad=8)
        ax.tick_params(colors="#94a3b8", labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#334155")
        ax.grid(True, linestyle=":", color="#1e293b", alpha=0.8)

    # 1. Side Elevation (Top-Left, matches Report Figure 2 by Kirill Filatov)
    ax_elev = fig.add_subplot(gs[0, 0])
    _style_cad_axes(ax_elev, "1. SIDE ELEVATION (SPAN = 3000 mm, HEIGHT = 1000 mm, TOP CHORD = 1550 mm)")
    xb = [0, 725, 1500, 2275, 3000]
    xt = [725, 1500, 2275]
    # Chords & Rakers
    ax_elev.plot([0, 3000], [0, 0], color="#38bdf8", lw=3.0, label="Chords / Rakers")
    ax_elev.plot([725, 2275], [1000, 1000], color="#38bdf8", lw=3.0)
    ax_elev.plot([0, 725], [0, 1000], color="#38bdf8", lw=3.0)
    ax_elev.plot([3000, 2275], [0, 1000], color="#38bdf8", lw=3.0)
    # Verticals & Fan Diagonals
    for x_v in xt:
        ax_elev.plot([x_v, x_v], [0, 1000], color="#e2e8f0", lw=2.0)
    ax_elev.plot([725, 1500], [1000, 0], color="#facc15", lw=2.0, label="Fan Web Diagonals")
    ax_elev.plot([2275, 1500], [1000, 0], color="#facc15", lw=2.0)
    # Dimension annotations
    ax_elev.annotate("", xy=(0, -180), xytext=(3000, -180), arrowprops=dict(arrowstyle="<->", color="#f43f5e", lw=1.2))
    ax_elev.text(1500, -260, "3000 mm (SPAN)", color="#f43f5e", ha="center", fontsize=8.5, fontweight="bold")
    ax_elev.annotate("", xy=(725, 1140), xytext=(2275, 1140), arrowprops=dict(arrowstyle="<->", color="#f43f5e", lw=1.2))
    ax_elev.text(1500, 1180, "1550 mm", color="#f43f5e", ha="center", fontsize=8.5, fontweight="bold")
    ax_elev.set_xlim(-250, 3250)
    ax_elev.set_ylim(-380, 1350)
    ax_elev.set_aspect("equal")
    ax_elev.legend(loc="upper right", facecolor="#1e293b", edgecolor="#334155", labelcolor="#f8fafc", fontsize=7.5)

    # 2. Plan View (Bottom-Left, matches Report Figure 4 by Bruce Hirwa)
    ax_plan = fig.add_subplot(gs[1, 0])
    _style_cad_axes(ax_plan, "2. PLAN VIEW — TRANSOMS (s = 1000 mm NOMINAL) & PLAN X-BRACING")
    ax_plan.plot([0, 3000], [0, 0], color="#38bdf8", lw=3.0)
    ax_plan.plot([0, 3000], [1400, 1400], color="#38bdf8", lw=3.0)
    for x_t in xb:
        ax_plan.plot([x_t, x_t], [0, 1400], color="#4ade80", lw=2.4, label="Floor Transoms" if x_t == 0 else "")
    for x1, x2 in [(0, 725), (725, 2275), (2275, 3000)]:
        ax_plan.plot([x1, x2], [0, 1400], color="#94a3b8", lw=1.4, linestyle="--", label="Plan Bracing" if x1 == 0 else "")
        ax_plan.plot([x1, x2], [1400, 0], color="#94a3b8", lw=1.4, linestyle="--")
    ax_plan.set_xlim(-250, 3250)
    ax_plan.set_ylim(-250, 1650)
    ax_plan.set_aspect("equal")
    ax_plan.legend(loc="upper right", facecolor="#1e293b", edgecolor="#334155", labelcolor="#f8fafc", fontsize=7.5)

    # 3. Cross-Section View (Right, matches Report Figure 3 by Ratanak Borat)
    ax_sec = fig.add_subplot(gs[:, 1])
    _style_cad_axes(ax_sec, "3. TYPICAL HALF-THROUGH CROSS-SECTION (U-FRAME ACTION)")
    # Left and Right Vertical Posts + Chords
    ax_sec.add_patch(patches.Rectangle((-50, 0), 100, 1000, edgecolor="#38bdf8", facecolor="#1e293b", lw=2))
    ax_sec.add_patch(patches.Rectangle((1350, 0), 100, 1000, edgecolor="#38bdf8", facecolor="#1e293b", lw=2))
    # Floor Transom (1400 mm c/c, 1300 mm clear width)
    ax_sec.add_patch(patches.Rectangle((50, -100), 1300, 100, edgecolor="#4ade80", facecolor="#14532d", lw=2, label="RHS Floor Transom"))
    # FRP Grating Deck (25 mm thick on top of transom)
    ax_sec.add_patch(patches.Rectangle((50, 0), 1300, 30, edgecolor="#facc15", facecolor="#713f12", hatch="xx", lw=1.5, label="FRP Grating Deck (1300 mm clear)"))
    ax_sec.annotate("", xy=(0, -220), xytext=(1400, -220), arrowprops=dict(arrowstyle="<->", color="#f43f5e", lw=1.2))
    ax_sec.text(700, -290, "1400 mm c/c (1500 mm OA / 1300 mm Clear)", color="#f43f5e", ha="center", fontsize=8.5, fontweight="bold")
    ax_sec.annotate("", xy=(-160, 0), xytext=(-160, 1000), arrowprops=dict(arrowstyle="<->", color="#f43f5e", lw=1.2))
    ax_sec.text(-220, 500, "1000 mm Parapet", color="#f43f5e", rotation=90, va="center", ha="center", fontsize=8.5, fontweight="bold")
    ax_sec.set_xlim(-350, 1650)
    ax_sec.set_ylim(-400, 1250)
    ax_sec.set_aspect("equal")
    ax_sec.legend(loc="upper center", facecolor="#1e293b", edgecolor="#334155", labelcolor="#f8fafc", fontsize=8)

    plt.savefig(output_path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == "__main__":
    export_ga_drawing_sheet()
