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
#import matplotlib.pyplot as plt
#from matplotlib import colors
#import pandas as pd
import os
import subprocess    
# Design variables of the problem


########### Filter
rmin = 0.018
nnz = 1000
MeshDir = "//scratch/bcrx/pranjan2/TopOpt/LBEAM/12MM/Filter/R_0.018/"

###########Constraint
sigmax = 47*(10**6)
pexp = 5   
sigrelax = 0.0001
volfrac = 1
compstar = 160
############Problem
InputFileName="LBeam_12mm"
nDV = 1028278



####NO continuation for this one, it's mass minimization problem#####
optIter = 500     #This is per-increment 
########################

penalty = 5

#Now remove the garbage density.dat, and write the new formatted one
with open('density.dat', 'w') as file:
   file.write("__X__")
   x0= 0.95 * np.ones(nDV,dtype=float)
var = InputVariable(x0, ArrayLabelReplacer("__X__", "\n"), 0, np.ones(nDV,dtype=float), lb=0.001, ub=1.0)

evalFun1 = ExternalRun("Direct",f"calTop.exe {InputFileName} -p {penalty} -r {rmin} -f {nnz} --pexp {pexp} --sigmin {sigmax} --sigrelax {sigrelax} --cg_eval NO", True)
evalFun1.addConfig("density.dat")
evalFun1.addData(InputFileName+".inp")
evalFun1.addData(MeshDir + "dnnz.bin")
evalFun1.addData(MeshDir +"drow.bin")
evalFun1.addData(MeshDir +"dval.bin")
evalFun1.addData(MeshDir +"dcol.bin")
evalFun1.addData(MeshDir +"dsum.bin")

fun1 = Function("Topop","Direct/objectives.csv",TableReader(0,3,(1,0),(None,None),","))
fun1.addInputVariable(var,"Direct/volume_sens.csv",TableReader(None,2,(1,0),(None,None),","))
fun1.addValueEvalStep(evalFun1)

# Stress constraint
con = Function("Topop","Direct/objectives.csv",TableReader(0,9,(1,0),(None,None),","))
con.addInputVariable(var,"Direct/stress_sens.csv",TableReader(None,0,(1,0),(None,None),","))

# Compliance constraint
con2 = Function("Topop","Direct/objectives.csv",TableReader(0,0,(1,0),(None,None),","))
con2.addInputVariable(var,"Direct/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))

driver = IpoptDriver()
driver.addObjective("min", fun1, 1)

driver.addUpperBound(con,1,1)

driver.addUpperBound(con2,compstar,1)
driver.setEvaluationMode(False,2.0)
driver.setStorageMode(False, "DSN_")
driver.setFailureMode("HARD")

nlp = driver.getNLP()

ncon = 2
lbMult = np.zeros(nDV)
ubMult= np.ones(nDV)
conMult = np.zeros(ncon)

nlp.set(warm_start_init_point = 'no' ,
            nlp_scaling_method = "gradient-based",    
            nlp_scaling_max_gradient=0.1,
            accept_every_trial_step = "no",
            limited_memory_max_history = 20,
            max_iter = optIter,
            tol = 1e-4,                     
            acceptable_iter = optIter,
            acceptable_tol = 1e-6,
            acceptable_obj_change_tol=1e-5,
            dual_inf_tol=1e-12,
            mu_strategy = "adaptive",
            mu_oracle = "quality-function", 
            mu_min = 1e-9,
            adaptive_mu_globalization="kkt-error",
            adaptive_mu_kkterror_red_iters = 2, 
            adaptive_mu_kkterror_red_fact = 0.999,
            adaptive_mu_kkt_norm_type="max-norm",
            fixed_mu_oracle="average_compl",
            print_timing_statistics = "yes",
            alpha_for_y = "primal",                
            output_file = 'ipopt_output.txt',
            print_level = 6)  

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

