# Helper Files 

## Overview
HelperFiles contains useful utility python files that helps with writing the .nam or .inp files for specific cases

## Before use
Before running calGeo, create an alias in your `.bashrc`:

``` sh
alias LineLoad='python3 some_path/Large-Scale-Topology-Optimization/HelperFiles/LineLoad.py'
alias DistributeLoad='python3 some_path/Large-Scale-Topology-Optimization/HelperFiles/DistributeLoad.py'
```
Thereafter source `.bashrc` and you may use the helperfiles at any directory.

## LineLoad.py
This file parses the .su2 file and find the nodes located on specific X, Y, and Z coordindates. ONLY 2 out of the 3 coordinates should be provided to define a line in the mesh where load is going to be applied. 
As an example 
``` sh
LineLoad CB.su2 NSurface.nam  --ycoord 1 --zcoord 2 
```
Will read `CB.su2`, find all nodes that lines on the line Y==1, Z==2, and write them in `NSurface.nam`. 

##DistributeLoad.py
This file takes in a total load magnitude and distribute it along the loaded points and modifies the inp file to reflect the distributed load
As an example
``` sh
DistributeLoad NSurface.nam CB.inp -10  2
```
Reads in NSurface.nam to record the number of loaded points. Then calculates the load per node by dividing the total load of (-10) by the number of loaded points. It then modifies CB.inp to find the line beneath `*CLOAD` and inputs the distributed load magnitude and the loading direction (2 stands for Y direction).

## PlotSens
``` sh
plotSens --compliance compliance_sens.csv --volume volume_sens.csv --cg cg_sens.csv
```

Parses the sensitivity CSVs from CalTop and plots a stacked bar chart