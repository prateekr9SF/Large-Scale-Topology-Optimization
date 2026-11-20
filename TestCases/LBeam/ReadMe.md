# LBeam Problem with Stress Constraints
## Problem Description 
This is a demonstration of mass minimization with and without stress constraints on a LBeam
## File Structure
The files are structured as follows:
```
LBeam
|-Mesh
|  |-Solid.su2
|  |-RunFilt.inp
|-LBeam_No_StxCon
|  |-RunOptim.py
|  |-LBeam_12mm.inp
|-LBeam_StxCon
|  |-RunOptim.py
|  |-LBeam_12mm.inp
```
Each `LBeam_*` is the root directory of the respective optimization problem.
`Mesh` is the location where the mesh and filtering matrix will be stored.

## Running The Problem
### Filtering
The first step is to construct the filtering matrix for the topology optimization problems. 
1. Navigate to  `Mesh`
2. Open `RunFilt.inp` and make sure that the second line points to the location of the mesh file.
3. Run ` Flags: -i jobname -r FILTER RADIUS -f FILTERNNZ -Passive 0/1 -SU2_MESH True/False `
   
   3.1 For this case, run `calFilt.exe RunFilt -r 0.018 -f 1000 -SU_MESH True`
   
   3.2 For parallel filtering, run `export OMP_NUM_THREADS=<NUMofTHREADS>` prior to running `calFilt`.

      3.2.1 If `-SU2_MESH True`, then `CalFilt.exe` will detect for a .su2 mesh file called `Solid.su2` within the directory and parse it into `mesh.nam, NSurface.nam, Nfix1.nam`. These mesh files will be read into the CCX input deck as`*INCLUDE`. 
   
   3.3 Once done, `dnnz.bin, drow.bin, dval.bin, dcol.bin, dsun.bin` will appear in the mesh directory.
   
   

### RunOptim.py
1. Before running the problem configure `RunOptim.py`. The input section is between the `Line8-Line40`. A few things to note is MeshDir should be the ABSOLUTE directory pointing to the Mesh directory that we are working with in the previous section.
2. Depending on the problem we are running, the objective function and constraints are formatted as `Functions`. For instance for the mass fraction minimization problem with stress and compliance constraint:
```
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
```
The `objectives.csv` is formatted as such:
  `COMPLIANCE	 FULL SOLID VOLUME	 DESIGN VOLUME	 VOLUME_FRACTION	 DISCRETENESS	 MASS	 CGx	 CGy	 CGz	 PNORM` 
  And therefore `TableReader(0,3...)` points to `VOLUME_FRACTION` and `TableReader(0,0...)` points to `COMPLIANCE`


3. The optimizer setup is toward the later part of the code using `nlp.set()`. For more references, checkout https://coin-or.github.io/Ipopt/OPTIONS.html

### LBeam.inp 
  The `*.inp` file are written similar in style of NASTRAN input files. Before running, make sure the 2-4 lines `*INCLUDE, INPUT=` points to the absolute location of the mesh file in `\Mesh`

### Running:
  After modifications and configurations outlined by the previous steps, run:
  `python RunOptim.py >log.opt &` to run the optimization problem. The `>log.opt` is optional, but very helpful in observing optimizer behaviour during the run.
  Similarly, you may run `export OMP_NUM_THREADS=<NUMofTHREADS>` before to use multiple threads.

## Results:
After running both problem, another folder `__WORKDIR__` should appear. This directory stores the last iteration's output from `calTop.exe`. There are the following that is worth investigating:
```
elastic_Field_active.vtu is the post-processing of the design domain
objectives.csv is the summary of objective value and constraint values
density.dat is the material distribution. It might be useful when doing other forms of post-processing or restarting the problem. It has to be used together with the meshes and skinelements as it constains no geometry definition and element index information
```

For a mass minimization problem, the target compliance is set to 160 Joule, and if stress constraint is used, the maximum allowable stress is set to 156.66MPa. The resultant structure using LOQO, with and without stress constaint, are shown below:
<img width="752" height="477" alt="image" src="https://github.com/user-attachments/assets/74a5ae6f-53ab-4fbe-ac6b-bc503154cbf0" />






