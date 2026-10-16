"""
sections.py
Published section properties for hot-finished S355J2H hollow sections to BS EN 10210-2
(SCI Blue Book / Tata Steel Celsius 355 tables).
Units converted from standard tabular units (mm, cm^2, cm^4, cm^3) to SI (m, m^2, m^4, m^3).
"""

# Tabular data per BS EN 10210-2:
# (Designation, Family, h_mm, b_mm, t_mm, A_cm2, Iy_cm4, Iz_cm4, IT_cm4, Wel_y_cm3, Wpl_y_cm3, Wel_z_cm3, Wpl_z_cm3)
EN10210_RAW_TABLE = [
    ("SHS 50x50x3.2",   "SHS",  50.0,  50.0, 3.2,  5.88,  21.7,  21.7,  34.4,  8.66, 10.4,  8.66, 10.4),
    ("SHS 60x60x3.2",   "SHS",  60.0,  60.0, 3.2,  7.16,  38.8,  38.8,  61.1, 12.90, 15.3, 12.90, 15.3),
    ("SHS 60x60x4.0",   "SHS",  60.0,  60.0, 4.0,  8.79,  46.1,  46.1,  73.4, 15.40, 18.4, 15.40, 18.4),
    ("SHS 70x70x3.6",   "SHS",  70.0,  70.0, 3.6,  9.42,  70.2,  70.2, 110.0, 20.10, 23.7, 20.10, 23.7),
    ("SHS 80x80x3.6",   "SHS",  80.0,  80.0, 3.6, 10.90, 107.0, 107.0, 167.0, 26.80, 31.5, 26.80, 31.5),
    ("SHS 80x80x4.0",   "SHS",  80.0,  80.0, 4.0, 12.00, 117.0, 117.0, 183.0, 29.10, 34.4, 29.10, 34.4),
    ("SHS 90x90x3.6",   "SHS",  90.0,  90.0, 3.6, 12.30, 155.0, 155.0, 240.0, 34.40, 40.2, 34.40, 40.2),
    ("SHS 90x90x4.0",   "SHS",  90.0,  90.0, 4.0, 13.60, 169.0, 169.0, 263.0, 37.50, 44.0, 37.50, 44.0),
    ("SHS 100x100x3.2", "SHS", 100.0, 100.0, 3.2, 12.30, 195.0, 195.0, 300.0, 39.00, 45.2, 39.00, 45.2),
    ("SHS 100x100x4.0", "SHS", 100.0, 100.0, 4.0, 15.20, 237.0, 237.0, 368.0, 47.50, 55.4, 47.50, 55.4),
    ("RHS 80x40x3.2",   "RHS",  80.0,  40.0, 3.2,  7.16,  57.2,  18.9,  44.1, 14.30, 17.8,  9.46, 11.2),
    ("RHS 80x40x4.0",   "RHS",  80.0,  40.0, 4.0,  8.79,  68.2,  22.2,  53.0, 17.10, 21.5, 11.10, 13.4),
    ("RHS 90x50x3.6",   "RHS",  90.0,  50.0, 3.6,  9.42,  98.3,  38.7,  87.5, 21.80, 26.9, 15.50, 18.3),
    ("RHS 100x40x4.0",  "RHS", 100.0,  40.0, 4.0, 10.40, 121.0,  26.7,  75.5, 24.30, 30.8, 13.40, 16.0),
    ("RHS 100x50x4.0",  "RHS", 100.0,  50.0, 4.0, 11.20, 140.0,  46.2, 108.0, 28.00, 34.7, 18.50, 21.9),
]

def load_section_catalog(min_t_mm=3.2):
    """Returns dictionary of published BS EN 10210-2 sections filtered by minimum wall thickness."""
    catalog = {}
    for row in EN10210_RAW_TABLE:
        name, fam, h_mm, b_mm, t_mm, A_cm2, Iy_cm4, Iz_cm4, IT_cm4, Wel_y, Wpl_y, Wel_z, Wpl_z = row
        if t_mm < min_t_mm:
            continue
        A = A_cm2 * 1e-4
        h = h_mm * 1e-3
        b = b_mm * 1e-3
        t = t_mm * 1e-3
        # Shear area Av per BS EN 1993-1-1 Cl. 6.2.6(3) for RHS/SHS: Av = A * h / (b + h)
        Av = A * (h / (b + h))
        catalog[name] = {
            "name": name,
            "family": fam,
            "h_mm": h_mm,
            "b_mm": b_mm,
            "t_mm": t_mm,
            "h": h,
            "b": b,
            "t": t,
            "A": A,
            "Iy": Iy_cm4 * 1e-8,
            "Iz": Iz_cm4 * 1e-8,
            "IT": IT_cm4 * 1e-8,
            "Wel_y": Wel_y * 1e-6,
            "Wpl_y": Wpl_y * 1e-6,
            "Wel_z": Wel_z * 1e-6,
            "Wpl_z": Wpl_z * 1e-6,
            "Av": Av,
            "mass_per_m": A * 7850.0
        }
    return catalog