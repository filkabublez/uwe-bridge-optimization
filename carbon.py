"""
carbon.py
Cradle-to-Gate (Modules A1-A3) Embodied Carbon Calculator.
Sources:
  - Inventory of Carbon & Energy (ICE) Database v3.0 (Hammond & Jones, 2019):
      * UK Hot-Rolled Structural Hollow Sections (S355J2H): 1.37 kgCO2e / kg steel (UK EAF/BOS mix)
      * Hot-Dip Galvanizing Zinc Coating (BS EN ISO 1461, min 85 um): 3.86 kgCO2e / kg zinc
"""

ICE_V3_STEEL_SHS_FACTOR = 1.37      # kgCO2e / kg steel (Modules A1-A3)
ICE_V3_ZINC_GALV_FACTOR = 3.86      # kgCO2e / kg zinc deposited
ZINC_COATING_KG_PER_M2 = 0.61       # 85 um zinc thickness @ 7140 kg/m^3 density

def compute_bridge_mass_and_carbon(members, catalog):
    total_vol_m3 = 0.0
    total_steel_mass_kg = 0.0
    total_surface_area_m2 = 0.0

    for m in members:
        sec = catalog[m["section"]]
        L = m["length_m"]
        total_vol_m3 += sec["A"] * L
        total_steel_mass_kg += sec["mass_per_m"] * L
        perimeter_m = 2.0 * (sec["h"] + sec["b"])
        total_surface_area_m2 += perimeter_m * L

    zinc_mass_kg = total_surface_area_m2 * ZINC_COATING_KG_PER_M2
    co2_steel_kg = total_steel_mass_kg * ICE_V3_STEEL_SHS_FACTOR
    co2_galv_kg = zinc_mass_kg * ICE_V3_ZINC_GALV_FACTOR

    return {
        "steel_volume_cm3": total_vol_m3 * 1e6,
        "steel_mass_kg": total_steel_mass_kg,
        "zinc_mass_kg": zinc_mass_kg,
        "co2_steel_A1_A3_kg": co2_steel_kg,
        "co2_galv_A1_A3_kg": co2_galv_kg,
        "co2_total_A1_A3_kg": co2_steel_kg + co2_galv_kg
    }