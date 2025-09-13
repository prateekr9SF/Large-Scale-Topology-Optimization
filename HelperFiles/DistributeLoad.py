# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 17:30:14 2025

@author: zheng
"""
import argparse
def distribute_load(filename_in, INPfile,Load,Vector):
    """
    Modifies the .inp file such that load is uniformly distributed along the loaded surface/line/points
    Parameters:
    -----------
    filename_in: str
        Path to the loaded points file, Default is Nsurface.nam
    fINPfile : str
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

    nodes_of_interest = []
    count = 0
    with open(filename_in, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0].isdigit():
                count+=1
    Dload = Load/count

    with open(INPfile, "r") as f:
        lines=f.readlines()
    for i, line in enumerate(lines):
        if line.strip().upper() =="*CLOAD":
            node_set = lines[i+1].strip().split(",")[0].strip()
            lines[i+1] = f"Nsurface, {int(Vector)}, {Dload:.4f}\n"

            break
    with open(INPfile, "w") as f:
        f.writelines(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Distribute load along loading points")
    parser.add_argument("filename_in", type=str, help=".nam file containing loaded nodes")
    parser.add_argument("INPfile", type=str, help="INP file ")
    parser.add_argument("Load", type=int, help="Direction of the Load")
    parser.add_argument("Vector", type=float, help="Magnitude of the Load")
    args = parser.parse_args()

    NumNodes = distribute_load(args.filename_in, args.INPfile,args.Load,args.Vector)