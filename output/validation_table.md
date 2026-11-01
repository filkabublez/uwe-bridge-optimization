| Verification Metric | Eurocode Hand Calc (Report) | Fusion 360 Solid FEA | Python 6-DOF Frame Solver | Explanation |
| :--- | :---: | :---: | :---: | :--- |
| **Floor Transom Bending Moment (`RHS 100x40x4.0`)** | 1.64 kNm | — | **1.50 kNm** | 8.5% difference: hand calc included higher self-weight allowance in $w_t = 7.76\text{ kN/m}$, whereas solver derives $w_t = 7.15\text{ kN/m}$ from exact deck + transom mass and models end rotational restraint |
| **Floor Transom Bending Stress & UR (`RHS 100x40x4.0`)** | 54.7 MPa (15.4%) | 54.63 MPa (15.4%) | **50.2 MPa (14.4%)** | Within 8%: ~5.5% from lower derived $M_{t,Ed}$ (1.50 vs 1.64 kNm) and ~2.6% from published EN 10210 $W_{pl,y} = 30.8\text{ cm}^3$ vs rounded $30.0\text{ cm}^3$ |
| **Main Chord Combined Stress (`SHS 100x100x4.0`)** | 4.8 MPa (1.3% axial) | 11.0–33.0 MPa (3.1–9.3%) | **14.0 MPa (3.5%)** | 6-DOF frame solver captures rigid welded joint secondary bending omitted by pin-truss hand calc |
