| Engineering Parameter | Baseline (UWE Report) | Sized (Full Braced) | Governing Limit / Clause |
| :--- | :---: | :---: | :--- |
| **Parapet Height $H$ (m)** | 1.00 m | 1.00 m | $\ge 1.00\text{ m}$ (`BS 6180` / `CD 353`) |
| **Chords / End Rakers (`S355J2H`)** | `SHS 100x100x4.0` | `SHS 50x50x3.2` | `BS EN 10210-2` ($t \ge 3.2\text{ mm}$) |
| **Web Verticals & Diagonals (`S355J2H`)** | `SHS 100x100x3.2` | `SHS 50x50x3.2` | `BS EN 10210-2` ($t \ge 3.2\text{ mm}$) |
| **Floor Transoms (`S355J2H`)** | `RHS 100x40x4.0` | `RHS 80x40x3.2` | `BS EN 10210-2` ($t \ge 3.2\text{ mm}$) |
| **Plan X-Bracing (`S355J2H`)** | `RHS 100x40x4.0` | `RHS 80x40x3.2` | Required for `ULS_Wind` & `ULS_Parapet` |
| **Structural Steel Mass (kg)** | 417.1 kg | **214.1 kg (-48.7%)** | Minimized by `optimise.py` |
| **Embodied Carbon A1–A3 ($\text{kgCO}_2\text{e}$)** | 606.5 $\text{kgCO}_2\text{e}$ | **315.1 $\text{kgCO}_2\text{e}$** | ICE v3.0 Steel + HDG Zinc |
| **Max Member Utilisation (`max_member_ur`)** | 14.4% | **25.2%** | $\le 85\%$ (`BS EN 1993-1-1`) |
| **Simplified Indicative Joint UR (`max_joint_ur`)** | 0.9% | **2.0%** | $\le 85\%$ (Indicative check, $\theta = 52.2^\circ$) |
| **Global Elastic Buckling Factor $\alpha_{cr}$** | 367.8 | **50.4** | $\ge 10.0$ (`BS EN 1993-1-1 Cl. 5.2.1`) |
| **SLS Vertical Deflection (`sls_defl_mm`)** | 0.08 mm | **0.19 mm** | $\le 10.0\text{ mm}$ ($L/300$, `Cl. 7.2`) |
| **1st Vertical Natural Freq $f_{1,v}$ (Hz)** | 134.7 Hz | **110.8 Hz** | $\ge 5.0\text{ Hz}$ (`BS EN 1990 Annex A2`) |
| **1st Lateral Natural Freq $f_{1,h}$ (Hz)** | 31.5 Hz | **21.4 Hz** | $\ge 2.5\text{ Hz}$ (`BS EN 1990 Annex A2`) |
