# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 17:30:14 2025

@author: zheng
"""

def parse_su2_and_record_nodes(filename_in, filename_out,ycoord,zcoord):
    nodes_of_interest = []

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
            if abs(y - ycoord)<=1e-4 and abs(z- zcoord)<= 1e-4:
                nodes_of_interest.append(idx)

    with open(filename_out, 'w') as f:
        f.write("*NSET, NSET=Nsurface\n")
        for idx in nodes_of_interest:
            f.write(f"{idx+1}\n")

    print(f"Number of loaded Point is {len(nodes_of_interest)}")

# Usage:
# Place your SU2 snippet in "example.su2"
# Adjust the Z target as needed
parse_su2_and_record_nodes("CB.su2", "NSurface.nam",0.0,2.0)