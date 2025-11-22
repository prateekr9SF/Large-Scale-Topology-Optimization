# -*- coding: utf-8 -*-
"""
Created on Sat Nov 15 19:19:35 2025

@author: wz10
"""

import numpy as np
import copy

def C_matrix(E, nu):
    """3D isotropic elasticity matrix"""
    lam = E*nu / ((1+nu)*(1-2*nu))
    mu  = E / (2*(1+nu))
    C = np.array([
        [lam+2*mu, lam,      lam,      0, 0, 0],
        [lam,      lam+2*mu, lam,      0, 0, 0],
        [lam,      lam,      lam+2*mu, 0, 0, 0],
        [0, 0, 0, mu, 0, 0],
        [0, 0, 0, 0, mu, 0],
        [0, 0, 0, 0, 0, mu]
    ])
    return C

def B_matrix(coords):
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]

    # Volume
    V = np.linalg.det(np.array([
        [1,x[0],y[0],z[0]],
        [1,x[1],y[1],z[1]],
        [1,x[2],y[2],z[2]],
        [1,x[3],y[3],z[3]]
    ])) / 6.0

    # Coefficients
    b = np.zeros(4)
    c = np.zeros(4)
    d = np.zeros(4)
    node_idx = [[0,1,2,3],[1,0,2,3],[2,0,1,3],[3,0,1,2]]
    for i in range(4):
        ii,j,k,l = node_idx[i]
        b[i] = y[j]*(z[k]-z[l]) + y[k]*(z[l]-z[j]) + y[l]*(z[j]-z[k])
        c[i] = z[j]*(x[k]-x[l]) + z[k]*(x[l]-x[j]) + z[l]*(x[j]-x[k])
        d[i] = x[j]*(y[k]-y[l]) + x[k]*(y[l]-y[j]) + x[l]*(y[j]-y[k])

    # Construct B matrix
    B = np.zeros((6,12))
    for i in range(4):
        B[0, 3*i]   = b[i]
        B[1, 3*i+1] = c[i]
        B[2, 3*i+2] = d[i]
        B[3, 3*i]   = c[i]
        B[3, 3*i+1] = b[i]
        B[4, 3*i+1] = d[i]
        B[4, 3*i+2] = c[i]
        B[5, 3*i]   = d[i]
        B[5, 3*i+2] = b[i]

    B /= (6*V)
    return B, V

def assemble_K(nodes, elems, E0, density, nu,penal):
    n_dof = 3*len(nodes)
    K = np.zeros((n_dof, n_dof))
    for e, conn in enumerate(elems):
        coords = nodes[conn]
        E = (density[e]**penal) * E0
        C = C_matrix(E, nu)
        B, V = B_matrix(coords)
        Ke = B.T @ C @ B * V
        dofs = []
        for n in conn:
            dofs.extend([3*n,3*n+1,3*n+2])
        for i in range(12):
            for j in range(12):
                K[dofs[i], dofs[j]] += Ke[i,j]
    return K

def apply_bc(K, F, bc):
    """bc: list of tuples (node_index, dof_index, value)"""
    for node,dof,val in bc:
        idx = 3*node + dof
        K[idx,:] = 0
        K[:,idx] = 0
        K[idx,idx] = 1
        F[idx] = val
    return K, F

def FEM_solver(nodes, elems, E0, density, nu, F, bc,penal):
    K = assemble_K(nodes, elems, E0, density, nu,penal)
    K, F = apply_bc(K, F, bc)
    U = np.linalg.solve(K,F)
    return U,K
def dKdrho(nodes, elems, E0, density, nu, F, bc,penal):
    K = assemble_K(nodes, elems, E0*penal, density, nu, (penal-1))
    K, F = apply_bc(K, F, bc)
    return K
def von_mises_energy(nodes, elems, density, E0, nu, U,penal):
    vm = np.zeros(len(elems))
    for e, conn in enumerate(elems):
        coords = nodes[conn]
        E = (density[e]**penal)*E0
        C = C_matrix(E,nu)
        B, V = B_matrix(coords)
        Vvm = np.array([
        [1, -0.5, -0.5, 0, 0, 0],
        [-0.5, 1, -0.5, 0, 0, 0],
        [-0.5, -0.5, 1, 0, 0, 0],
        [0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 3, 0],
        [0, 0, 0, 0, 0, 3]
                            ])
        dofs = []
        for n in conn:
            dofs.extend([3*n,3*n+1,3*n+2])
        ue = U[dofs]
        M = B.T @ C.T @Vvm@ C @ B 
        vm[e] = np.sqrt(ue.T @ M @ ue)
    return vm, M
def constructM_elem(nodes, elem, density, E0, nu, penal):
    coords = nodes[elem]
    E = (density**penal)*E0
    C = C_matrix(E,nu)
    B, V = B_matrix(coords)
    Vvm = np.array([
        [1, -0.5, -0.5, 0, 0, 0],
        [-0.5, 1, -0.5, 0, 0, 0],
        [-0.5, -0.5, 1, 0, 0, 0],
        [0, 0, 0, 3, 0, 0],
        [0, 0, 0, 0, 3, 0],
        [0, 0, 0, 0, 0, 3]
                            ])
    Me = B.T @ C.T @Vvm@ C @ B 
    return Me
    
# -------------------------------
# Example usage:
nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[2/3,2/3,2/3]]) # 5 nodes
elems = np.array([[0,1,2,3],[1,2,3,4]])                        # 2 tetra
density = np.array([1.0,1.0])
F = np.zeros(3*len(nodes))
F[4*3] = -0.57  # Apply force in x-direction on node 1
F[4*3+1] = -0.57  # Apply force in x-direction on node 1
F[4*3+2] = -0.57  # Apply force in x-direction on node 1
bc = [(0,0,0),(0,1,0),(0,2,0)]  # Fix node 0
rho = np.array([1.0,1.0])
rho1 = np.array([1.0,1.0])
########################    1 element problem
#nodes = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]) # 5 nodes
#elems = np.array([[0,1,2,3]])                        # 2 tetra
#F = np.zeros(3*len(nodes))
#F[0*3] = -0.57  # Apply force in x-direction on node 1
#F[0*3+1] = -0.57  # Apply force in x-direction on node 1
#F[0*3+2] = -0.57  # Apply force in x-direction on node 1
#bc = [(1,0,0),(1,1,0),(1,2,0),
#      (2,0,0),(2,1,0),(2,2,0),
#      (3,0,0),(3,1,0),(3,2,0)]  # Fix node 0
#rho = np.array([0.8])
#rho1 = np.array([1.0])
###################################################################
###Here are the variables to check
E0 = 4e6
nu = 0.3
sigmin=0.001
relax=0.001
Pexp = 1       #checked
penal = 4
#Do solid calculations

U,K = FEM_solver(nodes, elems, E0, rho, nu, F, bc,penal)
vm,M = von_mises_energy(nodes, elems,rho, E0, nu, U,penal)
sigE = vm/sigmin+np.ones(len(elems))*relax-np.ones(len(elems))*relax/rho
#correspond to eqn 18: Me0 should be M0

#print("Nodal displacements:\n", U.reshape(-1,3))

Pnorm = sum(sigE**Pexp)**(1/Pexp)
print("Pnorm:\n", Pnorm)
dPdrho_FD_his=np.zeros(len(elems))
dPdrho_ADJ_his=np.zeros(len(elems))
######################FD Calculation of d Pnorm d rho:#########################
h=1e-7
for i in range(len(elems)):
    rhoP= copy.deepcopy(rho)
    rhoP[i]=rhoP[i]-h
    Up,Kp = FEM_solver(nodes, elems, E0, rhoP, nu, F, bc,penal)
    vmP,Mp = von_mises_energy(nodes, elems, rho, E0, nu, Up,penal) 
    sigEP = vmP/sigmin+np.ones(len(elems))*relax-np.ones(len(elems))*relax/rhoP   #rho or rhoP?
    #Note here we uses curent Me and U=Up  for eqn 18
    PnormP=sum(sigEP**Pexp)**(1/Pexp)
    dPdrho_FD= (Pnorm-PnormP)/h
    dPdrho_FD_his[i]=dPdrho_FD
    print("dPnormd_rho_FD: ", dPdrho_FD)
######################ADJ Calculation of d Pnorm d rho:#######################
rhs=np.zeros(len(nodes)*3)
for i in range(len(elems)):
    node = elems[i]
    ue = np.zeros(12)
    for j in range(len(node)):
        N = node[j]
        uN =np.array( [U[N*3],U[N*3+1],U[N*3+2]])
        ue[j*3:j*3+3] = uN
    Me = constructM_elem(nodes, node, rho[i], E0, nu,penal)
    rhsT = (sigE[i]**(Pexp-1))/(sigmin*np.sqrt(ue.T@Me@ue))* (Me @ ue)   
    for k,j in enumerate(node):
        rhs[j*3]=rhs[j*3]+rhsT[k*3]
        rhs[j*3+1]=rhs[j*3+1]+rhsT[k*3+1]
        rhs[j*3+2]=rhs[j*3+2]+rhsT[k*3+2]
#This means Eqn 23 M0 are also current M, not M@rho=1
qb = np.linalg.solve(K,rhs)
for i in range(len(rhs)):
    print(f"RHS[{i}] = {rhs[i]:.4e}")
for i in range(len(rhs)):
    print(f"Q_tilda[{i}] = {qb[i]:.4e}")
for i in range(len(elems)):
    rhoP= copy.deepcopy(rho)*0.0
    rhoP[i]=rho[i]    #Corresponding to def rho
    dKdr = dKdrho(nodes, elems, E0, rhoP, nu, F, bc,penal)
    dPdrho_ADJ=Pnorm/(Pnorm**Pexp)*(-qb.T@dKdr@U)
    dPdrho_ADJ_his[i]=dPdrho_ADJ
    print("dPnormd_rho_ADJ: ", dPdrho_ADJ)    
relError = (dPdrho_ADJ_his-dPdrho_FD_his)/dPdrho_FD_his *100
print("Relative Error:\n", relError) 