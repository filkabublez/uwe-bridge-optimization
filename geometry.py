"""
geometry.py
Parametric 3D node and member topology generator for the UWE 3.0m span footbridge.
Enforces pedestrian parapet height >= 1.0m (CD 353 / BS 6180) and brace angle >= 30 deg (EN 1993-1-8).
"""
import numpy as np

def build_footbridge_geometry(span_m=3.0, width_cc_m=1.4, height_m=1.0,
                              top_chord_length_m=1.55, section_map=None,
                              include_mid_vertical=True, include_plan_bracing=True):
    if height_m < 1.0:
        raise ValueError(f"Height {height_m}m violates minimum 1.0m pedestrian parapet constraint.")

    x_sh = (span_m - top_chord_length_m) / 2.0
    x_mid = span_m / 2.0
    x_rt = span_m - x_sh

    # Check brace-to-chord angle >= 30 deg per BS EN 1993-1-8 Cl. 7.1.2
    raker_angle_deg = np.degrees(np.arctan2(height_m, x_sh))
    diag_angle_deg = np.degrees(np.arctan2(height_m, x_mid - x_sh))
    if min(raker_angle_deg, diag_angle_deg) < 30.0:
        raise ValueError(f"Brace angle {min(raker_angle_deg, diag_angle_deg):.1f} deg < 30 deg EN 1993-1-8 limit.")

    x_bot = [0.0, x_sh, x_mid, x_rt, span_m]
    x_top = [x_sh, x_mid, x_rt]

    nodes = []
    # 0..4: Left Bottom (y=0); 5..9: Right Bottom (y=width_cc_m)
    for y in [0.0, width_cc_m]:
        for xb in x_bot:
            nodes.append([xb, y, 0.0])
    # 10..12: Left Top (y=0, z=height_m); 13..15: Right Top (y=width_cc_m, z=height_m)
    for y in [0.0, width_cc_m]:
        for xt in x_top:
            nodes.append([xt, y, height_m])
    nodes = np.array(nodes, dtype=float)

    members = []
    for s in [0, 1]:
        b_off = s * 5
        t_off = 10 + s * 3
        for i in range(4):
            members.append({"n1": b_off+i, "n2": b_off+i+1, "group": "Bot_Chord", "section": section_map["Bot_Chord"]})
        for i in range(2):
            members.append({"n1": t_off+i, "n2": t_off+i+1, "group": "Top_Chord", "section": section_map["Top_Chord"]})
        members.append({"n1": b_off+0, "n2": t_off+0, "group": "End_Raker", "section": section_map["End_Raker"]})
        members.append({"n1": b_off+4, "n2": t_off+2, "group": "End_Raker", "section": section_map["End_Raker"]})

        # Outer vertical posts (form U-frames with transoms at x_sh and x_rt)
        members.append({"n1": b_off+1, "n2": t_off+0, "group": "Web_Vertical", "section": section_map["Web_Vertical"]})
        if include_mid_vertical:
            members.append({"n1": b_off+2, "n2": t_off+1, "group": "Web_Vertical", "section": section_map["Web_Vertical"]})
        members.append({"n1": b_off+3, "n2": t_off+2, "group": "Web_Vertical", "section": section_map["Web_Vertical"]})

        # Fan diagonals converging at mid-span bottom node
        members.append({"n1": t_off+0, "n2": b_off+2, "group": "Web_Diagonal", "section": section_map["Web_Diagonal"]})
        members.append({"n1": t_off+2, "n2": b_off+2, "group": "Web_Diagonal", "section": section_map["Web_Diagonal"]})

    # 4 Transverse Floor Transoms (at x=0, x_sh, x_rt, span_m; plus mid-span node 2-7)
    for nL, nR in [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]:
        members.append({"n1": nL, "n2": nR, "group": "Transom", "section": section_map["Transom"]})

    if include_plan_bracing:
        for bL1, bL2, bR1, bR2 in [(0, 1, 5, 6), (1, 3, 6, 8), (3, 4, 8, 9)]:
            members.append({"n1": bL1, "n2": bR2, "group": "Plan_Brace", "section": section_map["Plan_Brace"]})
            members.append({"n1": bR1, "n2": bL2, "group": "Plan_Brace", "section": section_map["Plan_Brace"]})

    for idx, m in enumerate(members):
        m["id"] = idx + 1
        p1, p2 = nodes[m["n1"]], nodes[m["n2"]]
        m["length_m"] = float(np.linalg.norm(p2 - p1))

    return nodes, members