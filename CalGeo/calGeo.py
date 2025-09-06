#!/usr/bin/env python3
import argparse
import numpy as np

# -------- Pretty banner (unchanged) --------
def print_banner():
    print("\n")
    print("   #####     ###    #       #####   #######   #####   ")
    print("  #         #   #   #      #       #         #     #  ")
    print("  #        #     #  #      #       #         #     #  ")
    print("  #        #######  #      #  ###  #####     #     #  ")
    print("  #        #     #  #      #     # #         #     #  ")
    print("  #     #  #     #  #      #     # #         #     #  ")
    print("   #####   #     #  #####   #####  #######   ##### #  ")
    print("\n")
    print(" A mesh pre-processing tool for CalTop.exe ")
    print("\n")
    print("* Contributors:")
    print("* Prateek Ranjan, Dept. of Aeronautics & Astronautics,")
    print("* Massachusetts Institute of Technology")
    print("* Wanzheng Zheng, Dept. of Aerospace Engineering,")
    print("* University of Illinois at Urbana-Champaign")
    print("* Kai A. James, Dept. of Aerospace Engineering,")
    print("* Georgia Institute of Technology")
    print("\n")

# -------- Utils --------
def signed_tet_volume(a, b, c, d):
    # a,b,c,d are 3-vectors
    return float(np.dot(np.cross(b - a, c - a), d - a) / 6.0)

def _dedupe_preserve_order(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def _parse_marker_spec(s):
    """
    Accepts 'name' or 'name:layers'. If layers omitted, defaults to 1.
    Examples: 'skin', 'skin:2', 'myTag:3'
    """
    if ":" in s:
        name, layer = s.split(":", 1)
        layer = int(layer.strip()) if layer.strip() else 1
        return name.strip(), max(1, layer)
    return s.strip(), 1

# -------- Core: stream-based mesh extraction --------
def extract_su2_mesh_data_with_element_index(su2_filepath, output_filepath="mesh.nam"):
    """
    Memory-light:
      - Stream file; keep only node array (float32).
      - Stream elements directly to output; don't store them.
    """
    print("Extracting mesh coordinates and element connectivity information ...", end='')

    # ---- Read nodes (single pass over just the node block) ----
    node_count = 0
    node_data = []
    with open(su2_filepath, "r") as f:
        node_section = False
        for line in f:
            if "NPOIN=" in line:
                node_count = int(line.split("=")[1].strip())
                node_section = True
                continue
            if node_section:
                if len(node_data) < node_count:
                    toks = line.split()
                    if len(toks) >= 3:
                        node_data.append([float(toks[0]), float(toks[1]), float(toks[2])])
                if len(node_data) == node_count:
                    break

    if not node_data or len(node_data) != node_count:
        raise ValueError("No node data or count mismatch. Check SU2 file format.")

    # shrink memory footprint: float32 is enough for coords
    node_data = np.asarray(node_data, dtype=np.float32)

    # ---- Open output & write nodes ----
    with open(output_filepath, "w") as out:
        out.write("**This file contains mesh information\n")
        out.write("*NODE, NSET=NALL\n")
        for i, (x, y, z) in enumerate(node_data, start=1):
            out.write(f"      {i}, {x}, {y}, {z}\n")

        out.write("\n*ELEMENT, TYPE=C3D4, ELSET=EALL\n")

        # ---- Stream elements (no in-memory list) ----
        nDV = 0
        element_section = False
        expected = None

        with open(su2_filepath, "r") as f:
            for line in f:
                if "NELEM=" in line:
                    expected = int(line.split("=")[1].strip())
                    element_section = True
                    continue
                if not element_section:
                    continue

                parts = line.strip().split()
                # Expect: <etype> n1 n2 n3 n4 <tag>
                if len(parts) == 6 and parts[0].isdigit():
                    n0 = [int(v) for v in parts[1:5]]  # SU2 indices assumed 0-based
                    a, b, c, d = (node_data[n0[0]], node_data[n0[1]], node_data[n0[2]], node_data[n0[3]])
                    vol = signed_tet_volume(a, b, c, d)
                    if vol < 0.0:
                        n0[1], n0[2] = n0[2], n0[1]
                    elif vol == 0.0:
                        raise ValueError("Zero-volume tetrahedron encountered.")

                    n1 = [k + 1 for k in n0]  # write 1-based
                    nDV += 1
                    out.write(f"    {nDV}, {n1[0]}, {n1[1]}, {n1[2]}, {n1[3]}\n")

                if expected is not None and nDV >= expected:
                    break

    print("done!")
    return nDV

# -------- Boundary node extraction (streamed) --------
def get_skin_nodes(su2_filepath, output_filepath="NSurface.nam"):
    """
    Marker 'skin' → writes *NSET, NSET=Nsurface with 1-based nodes (order-preserved).
    """
    print("Extracting skin node indices ...", end=' ')
    skin_section = False
    nodes = []
    with open(su2_filepath, "r") as f:
        for line in f:
            if "MARKER_TAG= skin" in line:
                skin_section = True
                continue
            if skin_section:
                if "MARKER_ELEMS=" in line:
                    continue
                if "MARKER_TAG=" in line:
                    break
                parts = line.strip().split()
                if len(parts) == 4 and parts[0].isdigit():
                    nodes.extend(int(x) + 1 for x in parts[1:])  # 1-based

    if not nodes:
        raise ValueError("Error: No nodes found for marker 'skin'.")

    ordered = _dedupe_preserve_order(nodes)
    with open(output_filepath, "w") as out:
        out.write("*NSET, NSET=Nsurface\n")
        for v in ordered:
            out.write(f"{v}\n")
    print("done!")

def get_traction_nodes(su2_filepath, output_filepath="NSurface.nam"):
    """
    Marker 'surface' → writes *NSET, NSET=Nsurface with 1-based nodes (order-preserved).
    """
    print("Extracting traction nodes ...", end=' ')
    skin_section = False
    nodes = []
    with open(su2_filepath, "r") as f:
        for line in f:
            if "MARKER_TAG= surface" in line:
                skin_section = True
                continue
            if skin_section:
                if "MARKER_ELEMS=" in line:
                    continue
                if "MARKER_TAG=" in line:
                    break
                parts = line.strip().split()
                if len(parts) == 4 and parts[0].isdigit():
                    nodes.extend(int(x) + 1 for x in parts[1:])  # 1-based

    if not nodes:
        raise ValueError("Error: No nodes found for marker 'surface'.")

    ordered = _dedupe_preserve_order(nodes)
    with open(output_filepath, "w") as out:
        out.write("*NSET, NSET=Nsurface\n")
        for v in ordered:
            out.write(f"{v}\n")
    print("done!")

def get_fixed_nodes(su2_filepath, output_filepath="Nfix1.nam"):
    """
    Marker 'fixed' → writes *NSET, NSET=Nfix1 with 1-based nodes (order-preserved).
    """
    print("Extracting fixed nodes ...", end='')
    root_section = False
    nodes = []
    with open(su2_filepath, "r") as f:
        for line in f:
            if "MARKER_TAG= fixed" in line:
                root_section = True
                continue
            if root_section:
                if "MARKER_ELEMS=" in line:
                    continue
                if "MARKER_TAG=" in line:
                    break
                parts = line.strip().split()
                if len(parts) == 4 and parts[0].isdigit():
                    nodes.extend(int(x) + 1 for x in parts[1:])  # 1-based

    if not nodes:
        raise ValueError("Error: No nodes found for marker 'fixed'.")

    ordered = _dedupe_preserve_order(nodes)
    with open(output_filepath, "w") as out:
        out.write("*NSET, NSET=Nfix1\n")
        for v in ordered:
            out.write(f"{v}\n")
    print("done!")

# -------- Skin element growth (stream per layer; no giant dict) --------
def extract_skin(su2_filepath, markername, layers=1):
    """
    Returns a set of *0-based* element ids touching the given marker,
    grown outward by 'layers' adjacency (including the first layer).
    Streams the element list each layer to keep memory flat.
    """
    # Collect base skin nodes (0-based internally)
    skin_nodes = set()
    skin_section = False
    with open(su2_filepath, "r") as f:
        for line in f:
            if f"MARKER_TAG= {markername}" in line:
                skin_section = True
                continue
            if skin_section:
                if "MARKER_ELEMS=" in line:
                    continue
                if "MARKER_TAG=" in line:
                    break
                parts = line.strip().split()
                if len(parts) == 4 and parts[0].isdigit():
                    skin_nodes.update(int(x) for x in parts[1:])  # 0-based

    if not skin_nodes:
        return set()

    selected = set()
    frontier = set(skin_nodes)
    layers = max(1, int(layers))

    for _ in range(layers):
        new_nodes = set()
        with open(su2_filepath, "r") as f:
            element_section = False
            for line in f:
                if "NELEM=" in line:
                    element_section = True
                    continue
                if not element_section:
                    continue
                parts = line.strip().split()
                if len(parts) == 6 and parts[0].isdigit():
                    eid = int(parts[-1])                # 0-based id from SU2
                    nodes = {int(x) for x in parts[1:5]}
                    if not nodes.isdisjoint(frontier):
                        selected.add(eid)
                        new_nodes.update(nodes)
        frontier = new_nodes

    return selected  # 0-based ids

def write_skin_element(skin_element_ids_zero_based, output_filepath="skinElementList.nam"):
    """
    Writes 1-based element ids, one per line.
    """
    with open(output_filepath, "w") as out:
        for eid0 in sorted(skin_element_ids_zero_based):
            out.write(f"{eid0 + 1}\n")

# -------- CLI --------
def main():
    parser = argparse.ArgumentParser(
        description="Convert SU2 mesh to CalculiX NAM format and extract boundary nodes."
    )
    parser.add_argument("su2_file", type=str, help="Path to the input SU2 mesh file")
    parser.add_argument(
        "SkinMarkerList",
        nargs='+',
        type=str,
        help="List of skin marker names or name:layers (e.g., skin or skin:2)."
    )
    parser.add_argument("--mesh_out", type=str, default="mesh.nam",
                        help="Output file for CalculiX mesh (default: mesh.nam)")
    parser.add_argument("--fixed_out", type=str, default="Nfix1.nam",
                        help="Output file for fixed node set (default: Nfix1.nam)")
    parser.add_argument("--surface_out", type=str, default="NSurface.nam",
                        help="Output file for surface node set (default: NSurface.nam)")
    args = parser.parse_args()

    print_banner()

    # Fixed/surface NSETs
    get_fixed_nodes(args.su2_file, args.fixed_out)
    get_traction_nodes(args.su2_file, args.surface_out)

    # Mesh (nodes + streamed elements)
    nDV = extract_su2_mesh_data_with_element_index(args.su2_file, args.mesh_out)

    # Skin elements (streamed, per layer)
    skinelement = set()
    for spec in args.SkinMarkerList:
        name, layers = _parse_marker_spec(spec)
        skinelement.update(extract_skin(args.su2_file, name, layers))

    if skinelement:
        write_skin_element(skinelement, output_filepath="skinElementList.nam")
        print("skin marker detected, writing skinElementList.nam")
    else:
        print("There are no skin markers in the SU2, no skinElementList.nam written")

    print("Number of design variables: ", nDV)

if __name__ == "__main__":
    main()
