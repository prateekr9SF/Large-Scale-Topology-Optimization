# Messerschmitt-Bölkow-Blohm (MBB) Beam Problem
## Problem Description 
This is a demonstration of Compliance minimization with volume fraction constraint and the Mass minimization with compliance constraint topology optimization problem. The original description of the problem can be found in _On CAD-integrated structural topology and design optimization_ by Niels Olhoff et.al. 

CalTop currently only works with tetrahedral elements, therefore the original problem is extended in the depth dimension and the elastic modulus of the material is scaled down such that the reference structure has the same stiffness as that described by Olhoff.

This demonstration uses a mesh with 192,955 elements with average edge length of 20mm. In light of the original description of the problem that "Upper and lower surface must be planar", a 1mm skin is added on the top and bottom surface of the beam and is considered as passive element to the topology optimization problem.

The results demonstrated are solved using filter radius `rmin=0.04 ` and non-zero element `nnz=1000`. The SIMP factor is chosen to be 3, and 300 iterations are solved.
## File Structure
The files are structured as follows:
```
MBB
|-Mesh
|  |-MBB_20MM-C3D4.nam
|  |-skinElementList.nam
|  |-RunFilt.inp
|-MBB_Comp
|  |-RunOptim.py
|  |-MBB.inp
|-MBB_Comp
|  |-RunOptim.py
|  |-MBB.inp
```
Each `MBB_*` is the root directory of the respective optimization problem.
`Mesh` is the location where the mesh and filtering matrix will be stored.

## Running The Problem
### Filtering
The first step is to construct the filtering matrix for the topology optimization problems. 
1. Navigate to  `Mesh`
2. Open `RunFilt.inp` and make sure that the second line points to the location of the mesh file.
3. Make sure `skinElementList.nam` is present in the same directory, `calFilt.exe` will automatically detect and create skinsets if there is a file named `skinElementList.nam` in the directory.
4. Run `calFilt.exe <INPNAME> -r <RMIN> -f <NNZ> -Passive <0 or 1>`
   
   4.1 For this case, run `calFilt.exe RunFilt -r 0.04 -f 1000 -Passive 1`
   
   4.2 For parallel filtering, run `export OMP_NUM_THREADS=<NUMofTHREADS>` prior to running `calFilt`.
   
   4.3 Once done, `dnnz.bin, drow.bin, dval.bin, dcol.bin, dsun.bin` will appear in the mesh directory.
   
   4.4 The argument `-Passive` argument will govern if the "passive element will be excluded during filtering". If `-Passive 1` ONLY active elements are considered during filtering, i.e. there will be a clean cut between the skin and internal structure. If `-Passive 0` then passive elements are considered during filtering, there will be a layer of partial density elements at the passive-active element interface due to filtering.
   

### RunOptim.py
1. Before running the problem configure `RunOptim.py`. The input section is between the `Line22-Line41`. A few things to note is MeshDir should be the ABSOLUTE directory pointing to the Mesh directory that we are working with in the previous section.
2. Depending on the problem we are running, the objective function and constraints are formatted as `Functions`. For instance for the volume fraction minimization problem:
```
   fun1 = Function("Topop","Direct/objectives.csv",TableReader(0,3,(1,0),(None,None),","))
   fun1.addInputVariable(var,"Direct/volume_sens.csv",TableReader(None,2,(1,0),(None,None),","))
   fun1.addValueEvalStep(evalFun1)


   con1 = Function("Topop","Direct/objectives.csv",TableReader(0,0,(1,0),(None,None),","))
   con1.addInputVariable(var,"Direct/compliance_sens.csv",TableReader(None,1,(1,0),(None,None),","))

   driver = IpoptDriver()
   driver.addObjective("min", fun1, 1)
   driver.addUpperBound(con1,CompStar,1)
```
The `objectives.csv` is formatted as such:
  `COMPLIANCE	 FULL SOLID VOLUME	 DESIGN VOLUME	 VOLUME_FRACTION	 DISCRETENESS	 MASS	 CGx	 CGy	 CGz	 PNORM` 
  And therefore `TableReader(0,3...)` points to `VOLUME_FRACTION` and `TableReader(0,0...)` points to `COMPLIANCE`


3. The optimizer setup is toward the later part of the code using `nlp.set()`. For more references, checkout https://coin-or.github.io/Ipopt/OPTIONS.html

### MBB.inp 
  The `*.inp` file are written similar in style of NASTRAN input files. Before running, make sure the 2nd line `*INCLUDE, INPUT=` points to the absolute location of the mesh file in `\Mesh`

### Running:
  After modifications and configurations outlined by the previous steps, run:
  `python RunOptim.py >log.opt &` to run the optimization problem. The `>log.opt` is optional, but very helpful in observing optimizer behaviour during the run.
  Similarly, you may run `export OMP_NUM_THREADS=<NUMofTHREADS>` before to use multiple threads.

## Results:
After running both problem, another folder `FinalCounter1` should appear. (If you are doing continuation, a `FinalCounter*` should appear for each step taken during continuation. This directory stores the last iteration's output from `calTop.exe`. Inside this director `FinalCounter1\Direct` There are the following that is worth investigating:
```
elastic_Field_active.vtu is the post-processing of the design domain
elastic_Field_passive.vtu is the post-processing of the passive domain (skin)
objectives.csv is the summary of objective value and constraint values
density.dat is the material distribution. It might be useful when doing other forms of post-processing or restarting the problem. It has to be used together with the meshes and skinelements as it constains no geometry definition and element index information
```

The compliance minimization is run first. The compliance resulted from the compliance-minimization problem is strategically applied to the volume fraction minimization problem as the compliance constraint. As we would expect, the two resulted in the almost the same material distribution. 

<img width="1514" height="625" alt="image" src="https://github.com/user-attachments/assets/9d655a61-a9a2-49eb-a34c-d91ff57287a3" />

It is also interesting to observe the convergence history of the two problems. A helpful python code is provided with calTop under `HelperFiles\ReadIPOPT.py`. Running that function with the log file as input with translate the log file into a csv file containing objective, optimality, feasibility, step size, and other informations for each iteration. Plotting these variables will provide insights on convergence behavior of the optimizer. 
