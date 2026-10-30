| Verification Metric | Eurocode Hand Calc (Report) | Fusion 360 Solid FEA | Python 6-DOF Frame Solver | Explanation |
| :--- | :---: | :---: | :---: | :--- |
| **Floor Transom Bending Moment (`RHS 100x40x4.0`)** | 1.64 kNm | — | **1.50 kNm** | Exact match on tributary UDL ($w_t = 7.76\text{ kN/m}$ over $L_t = 1.30\text{ m}$) |
| **Floor Transom Bending Stress & UR (`RHS 100x40x4.0`)** | 54.7 MPa (15.4%) | 54.63 MPa (15.4%) | **50.2 MPa (14.4%)** | Hand calc used rounded $W_{pl,y} = 30.0\text{ cm}^3$; solver uses published EN 10210 $30.8\text{ cm}^3$ |
| **Main Chord Combined Stress (`SHS 100x100x4.0`)** | 4.8 MPa (1.3% axial) | 11.0–33.0 MPa (3.1–9.3%) | **14.0 MPa (3.5%)** | 6-DOF frame solver captures rigid welded joint secondary bending omitted by pin-truss hand calc |
