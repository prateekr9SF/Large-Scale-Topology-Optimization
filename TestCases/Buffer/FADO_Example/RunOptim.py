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
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import os
       
# Design variables of the problem
# this defines initial value and how they are written to an arbitrary file
NCPU = 4
penalty = 3
rmin = 0.15
volfrac=0.3
InputFileName="Block"
nDV = 5412




x0=volfrac * np.ones(nDV,dtype=float)
with open("Constants.txt", mode="w+") as file:
    file.write(str(nDV)+"\n")
    file.write(str(volfrac)+"\n")
    file.write(str(rmin)+"\n")
    file.write(str(penalty)+"\n")
    file.write(InputFileName+"\n")
    file.write(str(NCPU)+"\n")
    
    
var = InputVariable(x0, ArrayLabelReplacer("__X__"), 0, np.ones(nDV,dtype=float), lb=0.0, ub=1.0)
parData1 = Parameter(["Constants.txt"],LabelReplacer("__DATA_FILE__"))

evalFun1 = ExternalRun("Direct","python3 ../../direct.py config_temp.txt")
evalFun1.addConfig("config_temp.txt")
evalFun1.addData("Constants.txt")
evalFun1.addData(InputFileName+".inp")
evalFun1.addParameter(parData1)

fun1 = Function("Topop","Direct/objectives.csv",TableReader(0,0,(1,0),(None,None),","))
fun1.addInputVariable(var,"Direct/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))
fun1.addValueEvalStep(evalFun1)
con = Function("Topop","Direct/objectives.csv",TableReader(0,2,(1,0),(None,None),","))
con.addInputVariable(var,"Direct/volume_sens.csv",TableReader(None,2,(1,0),(None,None),","))

driver = IpoptDriver()
driver.addObjective("min", fun1, 1)
driver.addUpperBound(con,volfrac)

optIter = 1000

driver.setEvaluationMode(False,2.0)
driver.setStorageMode(True, "DSN_")
driver.setFailureMode("HARD")

nlp = driver.getNLP()

ncon = 1
lbMult = np.zeros(nDV)
ubMult= np.ones(nDV)
conMult = np.zeros(ncon)


nlp.set(warm_start_init_point = 'no' ,
            nlp_scaling_method = "none",    
            accept_every_trial_step = "no",
            limited_memory_max_history = 50,
            max_iter = optIter,
            tol = 1e-4,                     
            acceptable_iter = optIter,
            acceptable_tol = 1e-6,
            acceptable_obj_change_tol=1e-5,
            dual_inf_tol=1e-06,
            mu_strategy = "adaptive",
            mu_oracle = "probing", 
            mu_min = 1e-9,
            adaptive_mu_globalization="kkt-error",
            adaptive_mu_kkterror_red_iters = 7, 
            adaptive_mu_kkterror_red_fact = 0.999,
            adaptive_mu_kkt_norm_type="max-norm",
            fixed_mu_oracle="average_compl",
            print_timing_statistics = "yes",
            alpha_for_y = "primal",                
            output_file = 'ipopt_output.txt')  

x, obj, status = nlp.solve(x0, mult_g = conMult, mult_x_L = lbMult, mult_x_U = ubMult)


driver.update()
# Print the optimized results---->

print("Primal variables solution")
print("x: ", x)

print("Bound multipliers solution: Lower bound")
print("lbMult: ", lbMult)

print("Bound multipliers solution: Upper bound")
print("ubMult: ", ubMult)

print("Constraint multipliers solution")
print("lambda:",conMult)