import adsk.core, adsk.fusion, traceback

# Choose which optimized bridge to generate in Fusion 360:
# Options: "TIER_1_CV_14PCT" (H=1000mm, 60x60x3.2 webs) or "TIER_2_WARREN_600MM" (H=600mm, Warren topology)
DESIGN_MODE = "TIER_1_CV_14PCT"

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        sketches = rootComp.sketches
        xyPlane = rootComp.xYConstructionPlane

        # Note: Fusion 360 internal API units are centimeters (cm)!
        # 3000 mm = 300.0 cm, 1400 mm = 140.0 cm
        L = 300.0
        W = 140.0

        if DESIGN_MODE == "TIER_1_CV_14PCT":
            H = 100.0       # 1000 mm height
            L_top = 144.0   # 1440 mm optimized top chord
            x_sh = (L - L_top) / 2.0
            x_mid = L / 2.0
            x_rt = L - x_sh
            x_bot = [0.0, x_sh, x_mid, x_rt, L]
            x_top = [x_sh, x_mid, x_rt]

            nodes = []
            for y in [0.0, W]:
                for xb in x_bot: nodes.append((xb, y, 0.0))
            for y in [0.0, W]:
                for xt in x_top: nodes.append((xt, y, H))

            # (n1, n2, outer_b_cm, outer_h_cm, wall_t_cm)
            bars = []
            for s in [0, 1]:
                bo, to = s*5, 10 + s*3
                for i in range(4): bars.append((bo+i, bo+i+1, 10.0, 10.0, 0.40)) # Bot chord 100x100x4
                for i in range(2): bars.append((to+i, to+i+1, 10.0, 10.0, 0.40)) # Top chord 100x100x4
                bars.append((bo+0, to+0, 9.0, 9.0, 0.36))                        # End raker 90x90x3.6
                bars.append((bo+4, to+2, 9.0, 9.0, 0.36))
                bars.append((bo+1, to+0, 7.0, 7.0, 0.32))                        # Verticals 70x70x3.2
                bars.append((bo+2, to+1, 7.0, 7.0, 0.32))
                bars.append((bo+3, to+2, 7.0, 7.0, 0.32))
                bars.append((to+0, bo+2, 6.0, 6.0, 0.32))                        # Diagonals 60x60x3.2
                bars.append((to+2, bo+2, 6.0, 6.0, 0.32))
            for nL, nR in [(0,5), (1,6), (3,8), (4,9)]:
                bars.append((nL, nR, 4.0, 10.0, 0.40))                           # Transoms 100x40x4
        else:
            H = 60.0        # 600 mm optimized height (L/H = 5:1)
            L_top = 156.0
            x_sh = (L - L_top) / 2.0
            x_mid = L / 2.0
            x_rt = L - x_sh
            x_bot = [0.0, L/3.0, 2.0*L/3.0, L]
            x_top = [x_sh, x_mid, x_rt]

            nodes = []
            for y in [0.0, W]:
                for xb in x_bot: nodes.append((xb, y, 0.0))
            for y in [0.0, W]:
                for xt in x_top: nodes.append((xt, y, H))

            bars = []
            for s in [0, 1]:
                bo, to = s*4, 8 + s*3
                for i in range(3): bars.append((bo+i, bo+i+1, 8.0, 8.0, 0.36))   # Bot chord 80x80x3.6
                for i in range(2): bars.append((to+i, to+i+1, 8.0, 8.0, 0.36))   # Top chord 80x80x3.6
                bars.append((bo+0, to+0, 8.0, 8.0, 0.36))
                bars.append((bo+3, to+2, 8.0, 8.0, 0.36))
                bars.append((to+0, bo+1, 6.0, 6.0, 0.32))                        # Warren Diags 60x60x3.2
                bars.append((bo+1, to+1, 6.0, 6.0, 0.32))
                bars.append((to+1, bo+2, 6.0, 6.0, 0.32))
                bars.append((bo+2, to+2, 6.0, 6.0, 0.32))
            for i in range(4):
                bars.append((i, i+4, 5.0, 9.0, 0.36))                            # Transoms 90x50x3.6

        # Draw 3D Centerline Sketch in Fusion 360
        sketch3D = sketches.add(xyPlane)
        sketch3D.name = f"Optimized_Bridge_{DESIGN_MODE}"
        lines = sketch3D.sketchCurves.sketchLines
        for n1, n2, b_cm, h_cm, t_cm in bars:
            p1 = adsk.core.Point3D.create(*nodes[n1])
            p2 = adsk.core.Point3D.create(*nodes[n2])
            lines.addByTwoPoints(p1, p2)

        ui.messageBox(f"Generated {len(bars)} parametric centerlines for {DESIGN_MODE}.\n"
                      f"Ready for Pipe/Sweep and Static Stress Simulation!")
    except:
        if ui:
            ui.messageBox("Failed:\n{}".format(traceback.format_exc()))
