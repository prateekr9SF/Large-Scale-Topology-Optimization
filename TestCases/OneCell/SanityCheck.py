import numpy as np

# ===========================
# Tetrahedral element B-matrix and volume
# ===========================
def tet_B_matrix(nodes):
    """
    nodes: 4x3 array of node coordinates
    Returns:
        B: 6x12 strain-displacement matrix
        V: volume of tetrahedron
    """
    x = nodes[:,0]
    y = nodes[:,1]
    z = nodes[:,2]

    J = np.array([
        [x[1]-x[0], x[2]-x[0], x[3]-x[0]],
        [y[1]-y[0], y[2]-y[0], y[3]-y[0]],
        [z[1]-z[0], z[2]-z[0], z[3]-z[0]]
    ])
    V = np.linalg.det(J)/6.0
    if V <= 0:
        raise ValueError("Negative tetrahedron volume!")

    Jinv = np.linalg.inv(J)

    # Gradients of shape functions in reference tetrahedron
    dN_hat = np.array([
        [-1,-1,-1],
        [ 1, 0, 0],
        [ 0, 1, 0],
        [ 0, 0, 1]
    ])
    dN = dN_hat @ Jinv.T  # 4x3

    # Build B-matrix 6x12
    B = np.zeros((6,12))
    for i in range(4):
        B[0,3*i]     = dN[i,0]
        B[1,3*i+1]   = dN[i,1]
        B[2,3*i+2]   = dN[i,2]
        B[3,3*i]     = dN[i,1]
        B[3,3*i+1]   = dN[i,0]
        B[4,3*i+1]   = dN[i,2]
        B[4,3*i+2]   = dN[i,1]
        B[5,3*i]     = dN[i,2]
        B[5,3*i+2]   = dN[i,0]
    return B, V

# ===========================
# Material matrix for isotropic linear elastic
# ===========================
def material_D(E, nu):
    lam = E * nu / ((1+nu)*(1-2*nu))
    mu  = E / (2*(1+nu))
    D = np.zeros((6,6))
    # normal terms
    D[0:3,0:3] = lam
    np.fill_diagonal(D[0:3,0:3], lam+2*mu)
    # shear terms (engineering strain)
    D[3:,3:] = mu*np.eye(3)
    return D

# ===========================
# FEM solver
# ===========================
def fem_2tet(nodes, elements, E0, rho, nu, load, fixed_nodes):
    n_nodes = nodes.shape[0]
    ndof = n_nodes * 3
    K = np.zeros((ndof, ndof))
    F = np.zeros(ndof)
    F += load

    # Assembly
    for elemI, elem in enumerate(elements):
        E = E0 * rho[elemI]
        coords = nodes[elem]
        B, V = tet_B_matrix(coords)
        D = material_D(E, nu)
        Ke = B.T @ D @ B *V

        # Map DOFs
        dofs = np.zeros(12, dtype=int)
        for i,n in enumerate(elem):
            dofs[3*i:3*i+3] = [3*n, 3*n+1, 3*n+2]

        # Assemble
        for i in range(12):
            for j in range(12):
                K[dofs[i],dofs[j]] += Ke[i,j]

    # Apply BCs
    for n in fixed_nodes:
        for d in range(3):
            fd = 3*n + d
            K[fd,:] = 0
            K[:,fd] = 0
            K[fd,fd] = 1
            F[fd] = 0

    # Solve
    u = np.linalg.solve(K,F)

    # Postprocess von Mises per element
    vm = np.zeros(len(elements))
    for elemI, elem in enumerate(elements):
        E = E0 * rho[elemI]
        coords = nodes[elements[elemI]]
        B, V = tet_B_matrix(coords)
        D = material_D(E, nu)

        # DOFs
        dofs = np.zeros(12, dtype=int)
        for i,n in enumerate(elements[elemI]):
            dofs[3*i:3*i+3] = [3*n, 3*n+1, 3*n+2]
        u_e = u[dofs]

        strain = B @ u_e
        stress = D @ strain
        sxx, syy, szz, sxy, syz, sxz = stress

        vm[elemI] = np.sqrt(
            0.5*((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2)
            + 3*(sxy**2 + syz**2 + sxz**2)
        )

    return u, vm

# ===========================
# Example input
# ===========================
nodes = np.array([
    [0,0,0],   # 0 fixed
    [1,0,0],   # 1
    [0,1,0],   # 2
    [0,0,1],   # 3
    [0.6666666666,0.6666666666,0.6666666666],   # 4 loaded
])

elements = np.array([
    [0,1,2,3],
    [1,2,3,4]
])

E0 = 4.0e6
rho = [1,1]
nu  = 0.33

# Unit load at node 4
load = np.zeros(3*len(nodes))
load[3*4:3*4+3] = [-0.57735,-0.57735,-0.57735]  # normal to face of nodes 2,3,4

# Node 0 fixed
fixed_nodes = [0]

# Solve
u, vm = fem_2tet(nodes, elements, E0, rho, nu, load, fixed_nodes)

print("Nodal displacements (u,v,w):\n", u.reshape(-1,3))
print("Von Mises stress per element:\n", vm)
