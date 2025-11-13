import numpy as np

# ---------------------------
# Step 0: Node coordinates
# ---------------------------
nodes = np.array([
    [0.0, 0.0, 0.0],  # Node 1
    [1.0, 0.0, 0.0],  # Node 2
    [0.0, 1.0, 0.0],  # Node 3
    [0.0, 0.0, 1.0]   # Node 4
])

# ---------------------------
# Step 1: Compute tetrahedron volume
# ---------------------------
v = np.linalg.det(np.array([
    nodes[1] - nodes[0],
    nodes[2] - nodes[0],
    nodes[3] - nodes[0]
])) / 6.0
print("Tetrahedron volume V =", v)

# ---------------------------
# Step 2: Compute B matrix
# ---------------------------
# Coefficients for linear tetrahedron
x = nodes[:, 0]
y = nodes[:, 1]
z = nodes[:, 2]

# Bi, Ci, Di for shape function derivatives
b = np.array([
    y[1]*(z[2]-z[3]) + y[2]*(z[3]-z[1]) + y[3]*(z[1]-z[2]),
    y[1]*(z[3]-z[0]) + y[2]*(z[0]-z[1]) + y[0]*(z[1]-z[3]),
    y[2]*(z[3]-z[0]) + y[3]*(z[0]-z[2]) + y[0]*(z[2]-z[3]),
    y[3]*(z[1]-z[0]) + y[1]*(z[0]-z[3]) + y[0]*(z[3]-z[1])
])

c = np.array([
    z[1]*(x[2]-x[3]) + z[2]*(x[3]-x[1]) + z[3]*(x[1]-x[2]),
    z[1]*(x[3]-x[0]) + z[2]*(x[0]-x[1]) + z[0]*(x[1]-x[3]),
    z[2]*(x[3]-x[0]) + z[3]*(x[0]-x[2]) + z[0]*(x[2]-x[3]),
    z[3]*(x[1]-x[0]) + z[1]*(x[0]-x[3]) + z[0]*(x[3]-x[1])
])

d = np.array([
    x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2]),
    x[1]*(y[3]-y[0]) + x[2]*(y[0]-y[1]) + x[0]*(y[1]-y[3]),
    x[2]*(y[3]-y[0]) + x[3]*(y[0]-y[2]) + x[0]*(y[2]-y[3]),
    x[3]*(y[1]-y[0]) + x[1]*(y[0]-y[3]) + x[0]*(y[3]-y[1])
])

B = np.zeros((6, 12))
for i in range(4):
    B[0, i*3+0] = b[i]
    B[1, i*3+1] = c[i]
    B[2, i*3+2] = d[i]
    B[3, i*3+0] = c[i]
    B[3, i*3+1] = b[i]
    B[4, i*3+1] = d[i]
    B[4, i*3+2] = c[i]
    B[5, i*3+0] = d[i]
    B[5, i*3+2] = b[i]

B /= (6 * v)  # Correct 1/(6V) factor

# ---------------------------
# Step 3: Material D matrix (isotropic)
# ---------------------------
E = 4e6
nu = 0.3
factor = E / ((1+nu)*(1-2*nu))
D = factor * np.array([
    [1-nu, nu, nu, 0, 0, 0],
    [nu, 1-nu, nu, 0, 0, 0],
    [nu, nu, 1-nu, 0, 0, 0],
    [0, 0, 0, (1-2*nu)/2, 0, 0],
    [0, 0, 0, 0, (1-2*nu)/2, 0],
    [0, 0, 0, 0, 0, (1-2*nu)/2]
])

# ---------------------------
# Step 4: Element stiffness matrix
# ---------------------------
K = B.T @ D @ B * v  # 12x12

# ---------------------------
# Step 5: Apply BCs and load
# ---------------------------
# Fixed nodes 2,3,4 => DOFs 3:12
free_dofs = np.array([0,1,2])  # Node 1 DOFs
F = np.array([-0.57735, -0.57735, -0.57735])  # N

K_reduced = K[np.ix_(free_dofs, free_dofs)]
u = np.linalg.solve(K_reduced, F)

# Build full displacement vector
U = np.zeros(12)
U[free_dofs] = u

# ---------------------------
# Step 6: Compute stress
# ---------------------------
epsilon = B @ U  # strain vector
sigma = D @ epsilon  # stress vector: [xx,yy,zz,xy,yz,xz]

# ---------------------------
# Step 7: Von Mises stress
# ---------------------------
sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz = sigma
sigma_vm = np.sqrt(
    0.5*((sigma_xx - sigma_yy)**2 + (sigma_yy - sigma_zz)**2 + (sigma_zz - sigma_xx)**2)
    + 3*(sigma_xy**2 + sigma_yz**2 + sigma_xz**2)
)

# ---------------------------
# Output results
# ---------------------------
print("Nodal displacements U1 =", U[0:3])
print("Stress [xx,yy,zz,xy,yz,xz] =", sigma)
print("Von Mises stress =", sigma_vm)

V = np.array([
    [2/3, -1/3, -1/3, 0, 0, 0],
    [-1/3, 2/3, -1/3, 0, 0, 0],
    [-1/3, -1/3, 2/3, 0, 0, 0],
    [0, 0, 0, 3, 0, 0],
    [0, 0, 0, 0, 3, 0],
    [0, 0, 0, 0, 0, 3]
])

# T = D * B
T = D @ B

# Energy-based matrix M
M = T.T @ V @ T  # multiply by volume

# Extract only the free DOFs (Node 1)
M_reduced = M[np.ix_(free_dofs, free_dofs)]

# Von Mises from u^T M u
sigma_vm_energy = np.sqrt(u.T @ M_reduced @ u)

print("Von Mises stress (energy-based) =", sigma_vm_energy)