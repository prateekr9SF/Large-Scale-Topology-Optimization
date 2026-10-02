from FADO import *
import subprocess
import numpy as np
import csv
import os
import glob
import shutil
import pandas as pd
import ipyopt
from numpy import ones, array, zeros

n_mpi=6
n_omp=4

penalty = 5
rmin = 0.012
volfrac=0.11
nDV_struct = 674115

nnz = 3000

# Config Master files
Elastic_ConfigMaster = "UCRMW_adj.inp"

# Fiter address
FilterAddress = "/home/prateek/Documents/scratch/UCRMW/Filter"


# Stress constraints
sigYield = 500*(10**6) # Yield stress for Al
FOS = 2.5              # Factor of Safety
sigmax = sigYield/FOS  # Max. allowable stress (sigma_l)
pexp = 10              # P-Norm exponent
sigrelax = 0.01       # Epsilon relaxation

with open("density.dat", "w") as file:
    file.write("__X__")

rho0 = volfrac * np.ones(nDV_struct, dtype=float)

rho = InputVariable(rho0,ArrayLabelReplacer("__X__", "\n"),0,np.ones(nDV_struct, dtype=float),lb=0.001,ub=1.0)

# Angle of attack
nDV_aoa = 1

aoa0 = np.array([0.0])

aoa = InputVariable(aoa0,PreStringHandler("AOA="),0,np.ones(nDV_aoa),lb=np.array([-4.0]),ub=np.array([2.0]))


# Combined optimizer vector
x0 = np.concatenate((rho0, aoa0))


nDV = len(x0)


# Replace %__DIRECT__% with an empty string when using enable_direct
enable_direct = Parameter([""], LabelReplacer("%__DIRECT__"))

# Replace *_def.cfg with *_FFD.cfg as required
mesh_in = Parameter(["MESH_FILENAME= UCRM9_1p2M_Fluid.su2"],\
       LabelReplacer("MESH_FILENAME= UCRM9_1p2M_Fluid.su2"))


# Use this wrapper for FSI problems
#evalFun1 = ExternalRun("Direct",f"python3 run_direct_MDA.py {n_mpi} {n_omp} {penalty} {rmin} {nnz}",True)

# MULTIDISCIPLIANRY ANALYSIS ---------------------------------------#
evalFun1 = ExternalRun("DIRECT",f"python3 run_direct_MDA.py "f"{n_mpi} {n_omp} {penalty} {rmin} {nnz} "f"{sigmax} {pexp} {sigrelax}",True)
evalFun1.addConfig("FSI_driver.py")
evalFun1.addConfig("run_direct_MDA.py")
evalFun1.addConfig("density.dat")
evalFun1.addData(f"{FilterAddress}/dnnz.bin")
evalFun1.addData(f"{FilterAddress}/drow.bin")
evalFun1.addData(f"{FilterAddress}/dval.bin")
evalFun1.addData(f"{FilterAddress}/dcol.bin")
evalFun1.addData(f"{FilterAddress}/dsum.bin")
evalFun1.addData("precice-config.xml")
evalFun1.addData("config.yml")
evalFun1.addConfig("Euler_UCRMW.cfg")
evalFun1.addData("UCRM9_1p2M_Fluid.su2")
evalFun1.addData("UCRM9_700K_Solid.su2")
evalFun1.addData("UCRMW.inp")


# SENSITIVITY ANALYSIS ---------------------------------------------#
evalSens = ExternalRun("SENS", f"OMP_NUM_THREADS={n_omp} "f"calFSI_ADJ.exe "f"-i UCRMW_adj "f"-p {penalty} "f"-r {rmin} "f"-f {nnz} "f"--pexp {pexp} "f"--sigmin {sigmax} "f"--sigrelax {sigrelax} "f"-precice-particiapant Calculix",True)
evalSens.addConfig(Elastic_ConfigMaster)
evalSens.addData("DIRECT/convergedRHS.dat")
evalSens.addData("DIRECT/Solid/skinElementList.nam")
evalSens.addData("DIRECT/Solid/mesh.nam")
evalSens.addData("DIRECT/Solid/NSurface.nam")
evalSens.addData("DIRECT/Solid/Nfix1.nam")
evalSens.addData(f"{FilterAddress}/dnnz.bin")
evalSens.addData(f"{FilterAddress}/drow.bin")
evalSens.addData(f"{FilterAddress}/dval.bin")
evalSens.addData(f"{FilterAddress}/dcol.bin")
evalSens.addData(f"{FilterAddress}/dsum.bin")
evalSens.addConfig("density.dat")




#evalFun1.addConfig("density.dat")
#evalFun1.addData("skinElementList.nam")



# FUNCTIONS ------------------------------------------------------------ #


# Drag --------------------------------------- #
#drag = Function("drag","Direct/Fluid/history.csv",LabeledTableReader('"CD"'))

#drag.addInputVariable(aoa,"Direct/Fluid/aero_derivatives.csv",TableReader(None, 7, (1,0), (None,None), ","))

#drag.addValueEvalStep(evalFun1)

#----------------------------------------------#

# Lift ----------------------------------------#
lift = Function("lift","DIRECT/Fluid/history.csv",LabeledTableReader('"CL"'))

lift.addInputVariable(aoa,"DIRECT/Fluid/aero_derivatives.csv",TableReader(None, 4, (1,0), (None,None), ","))

lift.addValueEvalStep(evalFun1)

#----------------------------------------------#


# Compliance ----------------------------------------#
#Compliance = Function("Compliance","Direct/Solid/objectives.csv",TableReader(0,0,(1,0),(None,None),","))

#Compliance.addInputVariable(rho,"Direct/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))

#Compliance.addValueEvalStep(evalFun1)


Compliance = Function("Compliance","DIRECT/Solid/objectives.csv",TableReader(0,0,(1,0),(None,None),","))
Compliance.addValueEvalStep(evalFun1)

Compliance.addInputVariable(rho,"SENS/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))
Compliance.addGradientEvalStep(evalSens)


#----------------------------------------------------#

# Volume --------------------------------------------#
volumeFraction = Function("volumeFraction","DIRECT/Solid/objectives.csv",TableReader(0,3,(1,0),(None,None),","))
volumeFraction.addValueEvalStep(evalFun1)

volumeFraction.addInputVariable(rho,"DIRECT/volume_sens.csv",TableReader(None,2,(1,0),(None,None),","))
volumeFraction.addGradientEvalStep(evalSens)

#----------------------------------------------------#

# Stress constraint ---------------------------------$
Stress = Function("Stress","DIRECT/Solid/objectives.csv",TableReader(0,9,(1,0),(None,None),","))
Stress.addValueEvalStep(evalFun1)

Stress.addInputVariable(rho,"SENS/stress_sens.csv",TableReader(None,0,(1,0),(None,None),","))
Stress.addGradientEvalStep(evalSens)



ncon = 4
driver = IpoptDriver()

# Minimize compliance
driver.addObjective("min", Compliance, 1)

# Subject to upper bound on CL @ 0.4 (2)
cl_limit_upper = 0.4
driver.addUpperBound(lift,cl_limit_upper)

# Subject to lower bound on CL @ 0.2 (3)
cl_limit_lower = 0.2
driver.addLowerBound(lift,cl_limit_lower)

# Subject to fixed volumefraction (4)
driver.addUpperBound(volumeFraction, volfrac)

# Subject to P-norm gradient (5)
driver.addUpperBound(Stress,1,1)

optIter = 0

driver.setEvaluationMode(False,2.0)
driver.setStorageMode(True, "DSN_")
driver.setFailureMode("HARD")

nlp = driver.getNLP()

print("Expected number of design variables:", nDV)
print("Initial x0 size:", x0.size)
print("rho0 size:", rho0.size)
print("aoa0 size:", aoa0.size)

lbMult = np.zeros(nDV)
ubMult= np.ones(nDV)
conMult = np.zeros(ncon)

nlp.set(warm_start_init_point = 'no' ,
            nlp_scaling_method = "none",    
            nlp_scaling_max_gradient=0.01,
            accept_every_trial_step = "yes",
            limited_memory_max_history = 50,
            max_iter = optIter,
            tol = 1e-4,                     
            acceptable_iter = optIter,
            acceptable_tol = 1e-6,
            acceptable_obj_change_tol=1e-5,
            dual_inf_tol=1e-06,
            mu_strategy = "adaptive",
            mu_oracle = "loqo", 
            mu_min = 1e-9,
            adaptive_mu_globalization="kkt-error",
            adaptive_mu_kkterror_red_iters = 4, 
            adaptive_mu_kkterror_red_fact = 0.999,
            adaptive_mu_kkt_norm_type="max-norm",
            fixed_mu_oracle="average_compl",
            print_timing_statistics = "yes",
            alpha_for_y = "primal",                
            output_file = 'ipopt_output.txt')  

x, obj, status = nlp.solve(x0, mult_g = conMult, mult_x_L = lbMult, mult_x_U = ubMult)

driver.update()
print(status)
# Print the optimized results---->

print("Primal variables solution")
print("x: ", x)

print("Bound multipliers solution: Lower bound")
print("lbMult: ", lbMult)

print("Bound multipliers solution: Upper bound")
print("ubMult: ", ubMult)

print("Constraint multipliers solution")
print("lambda:",conMult)