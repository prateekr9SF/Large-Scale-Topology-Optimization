# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 17:30:14 2025

@author: zheng
"""
import argparse
def parse_su2_and_record_nodes(filename_in, filename_out,xcoord,ycoord,zcoord):
    """
    parse through a .su2 mesh file and record the node indexes that meets X, Y and Z criterion
    the expected criterion is defined as a line, i.e. one of the xcoord, ycoord and zcoord should not be defined
    Example use:
    python3 LineLoad.py CB.su2 NSurface.nam  --ycoord 1 --zcoord 2 
    this gets you all the nodes that is on the line of Y==1 Z ==2
    Parameters:
    -----------
    filename_in: str
        Path to the .su2 file to parse
    filename_out : str
        path to the .nam file that record the node indicies, usually Nsurface.nam to meet CalTop syntax
    --xcoord: float, default []
        the X definition of the line 
    --ycoord: float, default []
        the Y definition of the line 
    --zcoord: float, default []
        the Z definition of the line 
    Returns:
    --------
    nodes_of_interest: Int 
        Number of nodes that fits the line criterion 
    None
        Saves a .nam file that have all the node indices of the line
    """

    nodes_of_interest = []
    coord_of_interest = []
    with open(filename_in, 'r') as f:
        lines = f.readlines()

    in_node_section = False
    for line in lines:
        if line.strip().startswith('NPOIN= '):
            in_node_section = True
            continue
        if in_node_section:
            if line.strip().startswith('%') or line.strip() == '':
                in_node_section = False
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            x = float(parts[0])
            y = float(parts[1])
            z = float(parts[2])
            idx = int(parts[3])
            if xcoord==[]:
                if abs(y - ycoord[0])<=1e-4 and abs(z- zcoord[0])<= 1e-4:
                    nodes_of_interest.append(idx)
                    coord_of_interest.append(x)
            elif ycoord==[]:
                if abs(x - xcoord[0])<=1e-4 and abs(z- zcoord[0])<= 1e-4:
                    nodes_of_interest.append(idx)
                    coord_of_interest.append(y)
            elif zcoord==[]:
                if abs(x - xcoord[0])<=1e-4 and abs(y- ycoord[0])<= 1e-4:
                    nodes_of_interest.append(idx)
                    coord_of_interest.append(z)

    with open(filename_out, 'w') as f:
        f.write("*NSET, NSET=Nsurface\n")
        for idx in nodes_of_interest:
            f.write(f"{idx+1}\n")
    with open('LoadCoord.dat','w')as f:
        for i in range(len(nodes_of_interest)):
            f.write(f"{nodes_of_interest[i]+1}, {coord_of_interest[i]}\n")
            

    print(f"Number of loaded Point is {len(nodes_of_interest)}")
    return nodes_of_interest

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse .su2 file for nodes on a line")
    parser.add_argument("filename_in", type=str, help="Path to su2 file")
    parser.add_argument("filename_out", type=str, help="Path to .nam file")
    parser.add_argument("--xcoord", nargs='*', type=float, default=[], help="X coord definition")
    parser.add_argument("--ycoord", nargs='*', type=float, default=[], help="Y coord definition")
    parser.add_argument("--zcoord", nargs='*', type=float, default=[], help="Z coord definition")
    args = parser.parse_args()

    NumNodes = parse_su2_and_record_nodes(args.filename_in, args.filename_out,args.xcoord,args.ycoord,args.zcoord)