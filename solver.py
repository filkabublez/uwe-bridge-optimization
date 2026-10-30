"""
solver.py
3D 6-DOF Euler-Bernoulli Frame Finite Element Solver and Linear Elastic Buckling Solver.
Checks stiffness matrix rank and condition number to detect mechanisms without artificial
diagonal regularization, recovers exact parabolic mid-span bending moments on UDL-loaded
transoms, and computes the critical elastic buckling load factor alpha_cr.
"""
import numpy as np
from scipy.linalg import eig


def _element_rotation_matrix_3d(p1, p2):
    """Computes 3x3 direction cosine rotation matrix R and 12x12 transformation T."""
    vec = p2 - p1
    L = float(np.linalg.norm(vec))
    cx = vec / L
    if abs(cx[0]) < 1e-6 and abs(cx[1]) < 1e-6:
        cz = np.array([1.0, 0.0, 0.0])
    else:
        up = np.array([0.0, 0.0, 1.0])
        cy = np.cross(up, cx)
        cy /= np.linalg.norm(cy)
        cz = np.cross(cx, cy)
    cy = np.cross(cz, cx)
    R = np.vstack([cx, cy, cz])
    T = np.zeros((12, 12), dtype=float)
    for i in range(4):
        T[3*i:3*(i+1), 3*i:3*(i+1)] = R
    return L, R, T


def _local_frame_stiffness_12x12(E, G, A, Iy, Iz, IT, L):
    """Standard 12x12 3D Euler-Bernoulli frame element stiffness matrix in local coordinates."""
    k = np.zeros((12, 12), dtype=float)
    ea_l = E * A / L
    k[0, 0] = k[6, 6] = ea_l
    k[0, 6] = k[6, 0] = -ea_l

    gj_l = G * IT / L
    k[3, 3] = k[9, 9] = gj_l
    k[3, 9] = k[9, 3] = -gj_l

    eiz = E * Iz
    k[1, 1] = k[7, 7] = 12.0 * eiz / L**3
    k[1, 7] = k[7, 1] = -12.0 * eiz / L**3
    k[1, 5] = k[5, 1] = k[1, 11] = k[11, 1] = 6.0 * eiz / L**2
    k[7, 5] = k[5, 7] = k[7, 11] = k[11, 7] = -6.0 * eiz / L**2
    k[5, 5] = k[11, 11] = 4.0 * eiz / L
    k[5, 11] = k[11, 5] = 2.0 * eiz / L

    eiy = E * Iy
    k[2, 2] = k[8, 8] = 12.0 * eiy / L**3
    k[2, 8] = k[8, 2] = -12.0 * eiy / L**3
    k[2, 4] = k[4, 2] = k[2, 10] = k[10, 2] = -6.0 * eiy / L**2
    k[8, 4] = k[4, 8] = k[8, 10] = k[10, 8] = 6.0 * eiy / L**2
    k[4, 4] = k[10, 10] = 4.0 * eiy / L
    k[4, 10] = k[10, 4] = 2.0 * eiy / L
    return k


def _local_geometric_stiffness_12x12(N_axial, L):
    """12x12 3D frame geometric stiffness matrix kg_loc (positive N_axial = tension)."""
    kg = np.zeros((12, 12), dtype=float)
    f = N_axial / L
    # Transverse y-deflection & theta_z coupling
    kg[1, 1] = kg[7, 7] = 1.2 * f
    kg[1, 7] = kg[7, 1] = -1.2 * f
    kg[1, 5] = kg[5, 1] = kg[1, 11] = k_11 = 0.1 * L * f
    kg[11, 1] = k_11
    kg[7, 5] = kg[5, 7] = kg[7, 11] = kg[11, 7] = -0.1 * L * f
    kg[5, 5] = kg[11, 11] = (2.0 / 15.0) * (L**2) * f
    kg[5, 11] = kg[11, 5] = -(1.0 / 30.0) * (L**2) * f
    # Transverse z-deflection & theta_y coupling
    kg[2, 2] = kg[8, 8] = 1.2 * f
    kg[2, 8] = kg[8, 2] = -1.2 * f
    kg[2, 4] = kg[4, 2] = kg[2, 10] = kg[10, 2] = -0.1 * L * f
    kg[8, 4] = kg[4, 8] = kg[8, 10] = kg[10, 8] = 0.1 * L * f
    kg[4, 4] = kg[10, 10] = (2.0 / 15.0) * (L**2) * f
    kg[4, 10] = kg[10, 4] = -(1.0 / 30.0) * (L**2) * f
    return kg


def solve_3d_frame(nodes, members, catalog, F, E=210e9, G=80.77e9, support_dofs=None,
                   member_udl_z=None, clear_deck_width_m=1.3):
    """
    Assembles and solves [K]{u} = {F} using 6-DOF 3D Frame Elements.
    Raises RuntimeError if [K_ff] is rank-deficient or ill-conditioned (mechanism).
    """
    n_nodes = len(nodes)
    ndof = 6 * n_nodes
    K = np.zeros((ndof, ndof), dtype=float)

    if support_dofs is None:
        support_dofs = [
            6*0+0, 6*0+1, 6*0+2, 6*0+3,  # Node 0: Pin
            6*5+0, 6*5+2,                # Node 5: Pin
            6*4+1, 6*4+2,                # Node 4: Roller
            6*9+2                        # Node 9: Roller
        ]

    elem_cache = []
    for m in members:
        p1, p2 = nodes[m["n1"]], nodes[m["n2"]]
        sec = catalog[m["section"]]
        L, R, T = _element_rotation_matrix_3d(p1, p2)
        k_loc = _local_frame_stiffness_12x12(E, G, sec["A"], sec["Iy"], sec["Iz"], sec["IT"], L)
        k_glob = T.T @ k_loc @ T
        dofs = [6*m["n1"]+j for j in range(6)] + [6*m["n2"]+j for j in range(6)]
        K[np.ix_(dofs, dofs)] += k_glob
        elem_cache.append((m, L, T, k_loc, dofs, sec))

    free_dofs = [d for d in range(ndof) if d not in support_dofs]
    K_ff = K[np.ix_(free_dofs, free_dofs)]

    rank = np.linalg.matrix_rank(K_ff)
    if rank < len(free_dofs):
        raise RuntimeError(f"Structural mechanism detected! Rank {rank} < {len(free_dofs)} free DOFs.")
    cond_num = np.linalg.cond(K_ff)
    if cond_num > 1e12:
        raise RuntimeError(f"Ill-conditioned stiffness matrix (cond = {cond_num:.2e}): mechanism detected.")

    u = np.zeros(ndof, dtype=float)
    u[free_dofs] = np.linalg.solve(K_ff, F[free_dofs])
    reactions = K @ u - F

    member_forces = []
    for m, L, T, k_loc, dofs, sec in elem_cache:
        u_glob = u[dofs]
        f_loc = k_loc @ (T @ u_glob)
        N_axial = float(f_loc[6])
        Vy_max = float(max(abs(f_loc[1]), abs(f_loc[7])))
        Vz_max = float(max(abs(f_loc[2]), abs(f_loc[8])))
        My_max = float(max(abs(f_loc[4]), abs(f_loc[10])))
        Mz_max = float(max(abs(f_loc[5]), abs(f_loc[11])))

        # Superpose distributed UDL mid-span bending moment & end shear for loaded transoms
        if member_udl_z and m["id"] in member_udl_z:
            w_z = float(member_udl_z[m["id"]])
            L_eff = clear_deck_width_m if m["group"] == "Transom" else L
            M_simply_supported = (w_z * (L_eff ** 2)) / 8.0
            My_max = max(My_max, M_simply_supported)
            Vz_max = max(Vz_max, 0.5 * w_z * L_eff)

        sigma_axial_MPa = (abs(N_axial) / sec["A"]) / 1e6
        sigma_bend_el_MPa = ((My_max / sec["Wel_y"]) + (Mz_max / sec["Wel_z"])) / 1e6
        sigma_bend_pl_MPa = ((My_max / sec["Wpl_y"]) + (Mz_max / sec["Wpl_z"])) / 1e6
        member_forces.append({
            "id": m["id"],
            "group": m["group"],
            "section": m["section"],
            "length_m": L,
            "N_Ed_N": N_axial,
            "V_Ed_N": max(Vy_max, Vz_max),
            "My_Ed_Nm": My_max,
            "Mz_Ed_Nm": Mz_max,
            "sigma_axial_MPa": sigma_axial_MPa,
            "sigma_total_MPa": sigma_axial_MPa + sigma_bend_el_MPa,
            "sigma_plastic_MPa": sigma_axial_MPa + sigma_bend_pl_MPa
        })

    return {
        "u": u,
        "K": K,
        "free_dofs": free_dofs,
        "reactions": reactions,
        "condition_number": cond_num,
        "member_forces": member_forces,
        "max_vertical_defl_mm": float(np.max(np.abs(u[2::6])) * 1000.0)
    }


def solve_linear_buckling(nodes, members, catalog, F, E=210e9, G=80.77e9, support_dofs=None,
                          member_udl_z=None, clear_deck_width_m=1.3):
    """
    Computes the lowest positive linear elastic critical buckling load factor alpha_cr
    by solving the generalized eigenvalue problem [K_ff]{v} = alpha_cr * [-Kg_ff]{v}.
    """
    sol = solve_3d_frame(
        nodes, members, catalog, F, E=E, G=G, support_dofs=support_dofs,
        member_udl_z=member_udl_z, clear_deck_width_m=clear_deck_width_m
    )
    ndof = 6 * len(nodes)
    Kg = np.zeros((ndof, ndof), dtype=float)

    for idx, m in enumerate(members):
        p1, p2 = nodes[m["n1"]], nodes[m["n2"]]
        L, R, T = _element_rotation_matrix_3d(p1, p2)
        N_axial = sol["member_forces"][idx]["N_Ed_N"]
        kg_loc = _local_geometric_stiffness_12x12(N_axial, L)
        kg_glob = T.T @ kg_loc @ T
        dofs = [6*m["n1"]+j for j in range(6)] + [6*m["n2"]+j for j in range(6)]
        Kg[np.ix_(dofs, dofs)] += kg_glob

    free_dofs = sol["free_dofs"]
    K_ff = sol["K"][np.ix_(free_dofs, free_dofs)]
    Kg_ff = Kg[np.ix_(free_dofs, free_dofs)]

    eigvals, _ = eig(K_ff, -Kg_ff)
    real_pos = [float(np.real(val)) for val in eigvals if np.isreal(val) and np.real(val) > 1e-3 and np.isfinite(np.real(val))]
    alpha_cr = min(real_pos) if real_pos else float("inf")
    return {"alpha_cr": alpha_cr, "static_solution": sol}
