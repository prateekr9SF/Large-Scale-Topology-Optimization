import argparse
import csv
import numpy as np
def SineLoad(filename_in, INPfile,Load,Vector):
    """
    Modifies the .inp file such that load is distributed in a sinusoidal fashion with the total magnitude specified
    Parameters:
    -----------
    filename_in: str
        Path to the loaded points file, Default is Nsurface.nam
    INPfile : str
        path to the .inp file to modify loads
    Load: float, default []
        Magnitude of the load
    Vector: int
        Vector of the loading direction
    Returns:
    --------
    nodes_of_interest: Int 
        Number of nodes that fits the line criterion 
    None
        Saves a .inp file that have modified loading condition
    """
    
    data = np.loadtxt(filename_in, delimiter=",")
    idx = data[:, 0]
    coord = data[:, 1]
    minC = min(coord)
    maxC = max(coord)
    L = maxC-minC
    SLoad = np.sin(coord/L*np.pi)
    SLoad = SLoad/np.sum(SLoad)*Load
    with open('SineLoad.inp','w')as f:
        for i in range(len(SLoad)):
            f.write(f"{int(idx[i]+1)}, {int(Vector)},{round(SLoad[i],4)}\n")
    with open(INPfile, "r") as f:
        lines=f.readlines()
    for i, line in enumerate(lines):
        if line.strip().upper() =="*CLOAD":
            lines[i+1] = "*INCLUDE INPUT=SineLoad.inp\n"

            break
    with open(INPfile, "w") as f:
        f.writelines(lines)

    return SLoad

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Distribute load along loading points")
    parser.add_argument("filename_in", type=str, help=".dat file containing loaded nodes")
    parser.add_argument("INPfile", type=str, help="INP file ")
    parser.add_argument("Load", type=float, help="Magnitude of the Load")
    parser.add_argument("Vector", type=int, help="Direction of the Load")
    args = parser.parse_args()

    NumNodes = SineLoad(args.filename_in, args.INPfile,args.Load,args.Vector)