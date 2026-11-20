#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 00:35:56 2025

@author: wz10
"""
import sys
sys.path.append("/home/prateek/Documents/GitHub/FADO_pyoptsparse")
import csv
import numpy as np
from FADO import *
import ipyopt
import pickle
#import matplotlib.pyplot as plt
#from matplotlib import colors
import pandas as pd
import os
       
# Design variables of the problem
# this defines initial value and how they are written to an arbitrary file
NCPU = 16
penalty = 5
rmin = 0.25
volfrac=0.30
InputFileName="SCB"
MeshDir = "/Meshes"
nDV =118511
nnz = 1000
RESTART = False
optIter = 40
# Set number of OpenMP threads
#os.environ["OMP_NUM_THREADS"] = str(NCPU)

with open('density.dat', 'w') as file:
   file.write("__X__")

if RESTART:
   with open("restart_var.pkl", "rb") as f:
      loaded_var = pickle.load(f)
   x0 = loaded_var["x"]
   conMult = loaded_var["conMult"]
   lbMult = loaded_var["lbMult"]
   ubMult = loaded_var["ubMult"]
else: 
   ncon = 7
   x0=volfrac * np.ones(nDV,dtype=float)
   lbMult = np.zeros(nDV)
   ubMult= np.ones(nDV)
   conMult = np.zeros(ncon)
    
    
var = InputVariable(x0, ArrayLabelReplacer("__X__", "\n" ), 0, np.ones(nDV,dtype=float), lb=0.001, ub=1.0)
evalFun1 = ExternalRun("Direct",f"calTop.exe {InputFileName} -p {penalty} -r {rmin} -f {nnz} --cg_eval YES", True)
evalFun1.addConfig("density.dat")

# Add filter files in advance
evalFun1.addData(MeshDir + "dnnz.bin")
evalFun1.addData(MeshDir +"drow.bin")
evalFun1.addData(MeshDir +"dval.bin")
evalFun1.addData(MeshDir +"dcol.bin")
evalFun1.addData(MeshDir +"dsum.bin")
evalFun1.addData(InputFileName+".inp")

fun1 = Function("Topop","Direct/objectives.csv",TableReader(0,0,(1,0),(None,None),","))
fun1.addInputVariable(var,"Direct/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))
fun1.addValueEvalStep(evalFun1)

# Volume fraction constraint
con1 = Function("Topop","Direct/objectives.csv",TableReader(0,3,(1,0),(None,None),","))
con1.addInputVariable(var,"Direct/volume_sens.csv",TableReader(None,2,(1,0),(None,None),","))

# CGx constraint
con2 = Function("Topop","Direct/objectives.csv",TableReader(0,6,(1,0),(None,None),","))
con2.addInputVariable(var,"Direct/cg_sens.csv",TableReader(None,0,(1,0),(None,None),","))

# CGy constraint
con3 = Function("Topop","Direct/objectives.csv",TableReader(0,7,(1,0),(None,None),","))
con3.addInputVariable(var,"Direct/cg_sens.csv",TableReader(None,1,(1,0),(None,None),","))

# CGz constraint
con4 = Function("Topop","Direct/objectives.csv",TableReader(0,8,(1,0),(None,None),","))
con4.addInputVariable(var,"Direct/cg_sens.csv",TableReader(None,2,(1,0),(None,None),","))

# CGy constraint
#con3 = Function("Topop","Direct/objectives.csv",TableReader(0,7,(1,0),(None,None),","))
#con3.addInputVariable(var,"Direct/cg_sens.csv",TableReader(None,1,(1,0),(None,None),","))

# CGz constraint
#con4 = Function("Topop","Direct/objectives.csv",TableReader(0,8,(1,0),(None,None),","))
#con4.addInputVariable(var,"Direct/cg_sens.csv",TableReader(None,2,(1,0),(None,None),","))

driver = IpoptDriver()
driver.addObjective("min", fun1, 1)
driver.addUpperBound(con1,volfrac)

# Constrain CG-x between 0.01 and 5.01 m
driver.addUpperBound(con2,5.01, 1)
driver.addLowerBound(con2,4.99, 1)

# Constrain CG-y between 4.99 and 5.01 m
driver.addUpperBound(con3,5.01, 1)
driver.addLowerBound(con3,4.99, 1)

# Constrain CG-z between 4.99 and 5.01 m
driver.addUpperBound(con4,0.5, 1)
driver.addLowerBound(con4,0.48, 1)

driver.setEvaluationMode(False,2.0)
driver.setStorageMode(False, "DSN_")
driver.setFailureMode("HARD")

nlp = driver.getNLP()

#ncon = 4
#lbMult = np.zeros(nDV)
#ubMult= np.ones(nDV)
#conMult = np.zeros(ncon)


nlp.set(warm_start_init_point = 'yes' if RESTART else 'no',
            warm_start_entire_iterate = 'yes' if RESTART else 'no',
            nlp_scaling_method = "gradient-based",
            nlp_scaling_max_gradient = 1e-4,
            accept_every_trial_step = "no",
            limited_memory_max_history = 20,
            max_iter = optIter,
            tol = 1e-12,                     
            acceptable_iter = optIter,
            acceptable_tol = 1e-8,
            acceptable_obj_change_tol=1e-5,
            dual_inf_tol=1e-06,
            mu_strategy = "adaptive",
            mu_oracle = "quality-function", 
            mu_min = 1e-9,
            adaptive_mu_globalization="kkt-error",
            adaptive_mu_kkterror_red_iters = 6, 
            adaptive_mu_kkterror_red_fact = 0.999,
            adaptive_mu_kkt_norm_type="max-norm",
            fixed_mu_oracle="average_compl",
            print_timing_statistics = "yes",
            alpha_for_y = "primal",
            print_level = 5,                 
            output_file = 'ipopt_output.txt',
            warm_start_bound_push=1e-4,
            warm_start_mult_bound_push=1e-4,
            warm_start_slack_bound_push=1e-4,
            warm_start_bound_frac=1e-4)   

x, obj, status = nlp.solve(x0, mult_g = conMult, mult_x_L = lbMult, mult_x_U = ubMult)


driver.update()
# Print the optimized results---->
print("Writing Results in to restart_var.pkl....\n")
restart_var= {
    "x": x,
    "conMult": conMult,
    "lbMult": lbMult,
    "ubMult":ubMult
}
with open("restart_var.pkl", "wb") as f:
    pickle.dump(restart_var, f)

print("Primal variables solution")
print("x: ", x)

print("Bound multipliers solution: Lower bound")
print("lbMult: ", lbMult)

print("Bound multipliers solution: Upper bound")
print("ubMult: ", ubMult)

print("Constraint multipliers solution")
print("lambda:",conMult)
