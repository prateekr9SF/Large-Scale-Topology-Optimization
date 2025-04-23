
import os
from collections import defaultdict

# Final correction: Add element index in the first column while removing the first and last columns from the SU2 file


def extract_su2_mesh_data_with_element_index(su2_filepath, output_filepath="mesh.nam"):
    """
    Extracts node coordinates and element connectivity from a SU2 file,
    writes them to mesh.nam in a structured format with a header and footer,
    and ensures correct formatting for element data by adding an index as the first column 
    while removing the original first and last columns.
    
    Parameters:
    su2_filepath (str): Path to the input SU2 file.
    output_filepath (str): Path to the output text file (default: mesh.nam).
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()
    
    # Extract node coordinates
    node_section = False
    node_count = 0
    node_data = []
    
    for line in lines:
        if "NPOIN=" in line:
            node_count = int(line.split("=")[1].strip())
            node_section = True
            continue
        
        if node_section:
            if len(node_data) < node_count:
                node_data.append(line.strip())
            else:
                node_section = False  # Stop after collecting all nodes
    
    # Extract element connectivity
    element_section = False
    element_count = 0
    element_data = []
    
    for line in lines:
        if "NELEM=" in line:
            element_count = int(line.split("=")[1].strip())
            element_section = True
            continue
        
        if element_section:
            if len(element_data) < element_count:
                elements = line.strip().split()
                if len(elements) > 2:  # Ensure valid element data with more than two entries
                    intelem = [int(x)+1 for x in elements[1:-1]]
                    formatted_element = ",".join(str(x) for x in intelem)  # Exclude first and last column
                    element_data.append(formatted_element)
            else:
                break  # Stop after collecting all elements
    
    # Ensure node data was extracted correctly
    if not node_data:
        return "Error: No node data extracted. Check SU2 file format."

    # Write extracted data to the output file
    with open(output_filepath, "w") as output_file:
        output_file.write("**This file contains mesh information\n")
        output_file.write("*NODE, NSET=NALL\n")
        
        # Writing nodal coordinates with node index starting from 1
        for idx, node in enumerate(node_data, start=1):
            coords = node.split()  # Assume space-separated coordinates
            formatted_line = f"{idx}, {coords[0]}, {coords[1]}, {coords[2]}"
            output_file.write(formatted_line + "\n")
        
        # Add a blank line before the element section
        output_file.write("\n")
        output_file.write("*ELEMENT, TYPE=C3D4, ELSET=EALL\n")
        
        # Writing element connectivity, adding element index as the first column
        for idx, elem in enumerate(element_data, start=1):
            formatted_element = f"{idx},{elem}"
            output_file.write(formatted_element + "\n")

    #print (f"Successfully formatted {len(node_data)} nodes and {len(element_data)} elements with indexed elements, saved to {output_filepath}"_

# Function to extract node indices associated with the "root" marker without sorting and with a comma per line
def get_skin_nodes(su2_filepath, output_filepath="NSurface.nam"):
    """
    Extracts node indices associated with the marker 'skin' from a SU2 file,
    offsets the node indices by 1, and writes them to Nsurface.nam with one entry per line,
    ensuring that each value is separated by a comma and maintaining the original order.
    
    Ignores only the first column in the marker element data.
    
    Parameters:
    su2_filepath (str): Path to the input SU2 file.
    output_filepath (str): Path to the output text file (default: Nfix1.nam).
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()

    skin_section = False
    skin_nodes = []  # Use a list to maintain original order

    for line in lines:
        # Identify the marker for "skin" (lowercase as found in the file)
        if "MARKER_TAG= skin" in line:
            skin_section = True
            continue  # Move to the next line where elements start
        
        # If in ROOT section, read node indices from element connectivity
        if skin_section:
            if "MARKER_ELEMS=" in line:  # Found number of elements in the marker
                continue
            elif "MARKER_TAG=" in line:  # Reached a new marker, stop reading
                break
            else:
                elements = line.strip().split()
                if len(elements) == 4:  # Ensure it's a valid element line
                    node_indices = elements[1:]  # Exclude only the first column
                    skin_nodes.extend(map(lambda x: int(x) + 1, node_indices))  # Offset by 1, maintain order

    # Ensure nodes were extracted
    if not skin_nodes:
        return "Error: No nodes found for marker 'skin'. Check SU2 file format."

    # Write the extracted ROOT node indices to the output file with one entry per line, comma-separated
    with open(output_filepath, "w") as output_file:
        output_file.write("*NSET, NSET=Nfix1\n")
        
        # Write nodes one per line, with a comma after each entry
        for node in skin_nodes:
            output_file.write(f"{node},\n")

    #return f"Successfully extracted {len(skin_nodes)} skin nodes with offset, written one per line with comma (original order preserved), saved to {output_filepath}"

def get_traction_nodes(su2_filepath, output_filepath="NSurface.nam"):
    """
    Extracts node indices associated with the marker 'surface' from a SU2 file,
    offsets the node indices by 1, and writes them to Nsurface.nam with one entry per line,
    ensuring that each value is separated by a comma and maintaining the original order.
    
    Ignores only the first column in the marker element data.
    
    Parameters:
    su2_filepath (str): Path to the input SU2 file.
    output_filepath (str): Path to the output text file (default: Nfix1.nam).
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()

    skin_section = False
    skin_nodes = []  # Use a list to maintain original order

    for line in lines:
        # Identify the marker for "skin" (lowercase as found in the file)
        if "MARKER_TAG= surface" in line:
            skin_section = True
            continue  # Move to the next line where elements start
        
        # If in ROOT section, read node indices from element connectivity
        if skin_section:
            if "MARKER_ELEMS=" in line:  # Found number of elements in the marker
                continue
            elif "MARKER_TAG=" in line:  # Reached a new marker, stop reading
                break
            else:
                elements = line.strip().split()
                if len(elements) == 4:  # Ensure it's a valid element line
                    node_indices = elements[1:]  # Exclude only the first column
                    skin_nodes.extend(map(lambda x: int(x) + 1, node_indices))  # Offset by 1, maintain order

    # Ensure nodes were extracted
    if not skin_nodes:
        return "Error: No nodes found for marker 'surface'. Check SU2 file format."
    new_skin_nodes = list(set(skin_nodes))
    # Write the extracted surface node indices to the output file with one entry per line, comma-separated
    with open(output_filepath, "w") as output_file:
        output_file.write("*NSET, NSET=Nsurface\n")
        
        # Write nodes one per line, with a comma after each entry
        for node in new_skin_nodes:
            output_file.write(f"{node},\n")

    #return f"Successfully extracted {len(skin_nodes)} skin nodes with offset, written one per line with comma (original order preserved), saved to {output_filepath}"


# Function to extract node indices associated with the "fixed" marker without sorting and with a comma per line
def get_fixed_nodes(su2_filepath, output_filepath="Nfix1.nam"):
    """
    Extracts node indices associated with the marker 'fixed' from a SU2 file,
    offsets the node indices by 1, and writes them to Nfix1.nam with one entry per line,
    ensuring that each value is separated by a comma and maintaining the original order.
    
    Ignores only the first column in the marker element data.
    
    Parameters:
    su2_filepath (str): Path to the input SU2 file.
    output_filepath (str): Path to the output text file (default: Nfix1.nam).
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()

    root_section = False
    root_nodes = []  # Use a list to maintain original order

    for line in lines:
        # Identify the marker for "fixed" (lowercase as found in the file)
        if "MARKER_TAG= fixed" in line:
            root_section = True
            continue  # Move to the next line where elements start
        
        # If in fixed section, read node indices from element connectivity
        if root_section:
            if "MARKER_ELEMS=" in line:  # Found number of elements in the marker
                continue
            elif "MARKER_TAG=" in line:  # Reached a new marker, stop reading
                break
            else:
                elements = line.strip().split()
                if len(elements) == 4:  # Ensure it's a valid element line
                    node_indices = elements[1:]  # Exclude only the first column
                    root_nodes.extend(map(lambda x: int(x) + 1, node_indices))  # Offset by 1, maintain order

    # Ensure nodes were extracted
    if not root_nodes:
        return "Error: No nodes found for marker 'fixed'. Check SU2 file format."

    # Write the extracted fixed node indices to the output file with one entry per line, comma-separated
    new_root_nodes = list(set(root_nodes))
    with open(output_filepath, "w") as output_file:
        output_file.write("*NSET, NSET=Nfix1\n")
        
        # Write nodes one per line, with a comma after each entry
        for node in new_root_nodes:
            output_file.write(f"{node},\n")

    #return f"Successfully extracted {len(root_nodes)} ROOT nodes with offset, written one per line with comma (original order preserved), saved to {output_filepath}"

def extract_skin_node_triplets(su2_filepath):
    """
    Extracts node triplets associated with the marker 'skin' from a SU2 file.
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()
    
    skin_section = False
    skin_triplets = []
    
    for line in lines:
        if "MARKER_TAG= skin" in line:
            skin_section = True
            continue 
        
        if skin_section:
            if "MARKER_ELEMS=" in line:
                continue
            elif "MARKER_TAG=" in line:  
                break
            else:
                elements = line.strip().split()
                if len(elements) == 4:
                    node_triplet = tuple(map(int, elements[1:]))
                    skin_triplets.append(node_triplet)
    return list(set(skin_triplets))


def find_all_tetrahedral_elements_for_skin_optimized(su2_filepath, skin_triplets, output_filepath="skin_tetra_elements.nam"):
    """
    Finds the tetrahedral element IDs that contain the given node triplets from the skin marker.
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()
    
    tetrahedral_elements = {}
    node_to_elements = defaultdict(set)
    element_section = False
    
    for line in lines:
        if "NELEM=" in line:
            element_section = True
            continue
        
        if element_section:
            elements = line.strip().split()
            if len(elements) == 6 and elements[0].isdigit():
                element_id = int(elements[-1])
                node_indices = set(map(int, elements[1:5]))
                tetrahedral_elements[element_id] = node_indices
                for node in node_indices:
                    node_to_elements[node].add(element_id)
    
    skin_to_tetra_mapping = {}
    
    for triplet in skin_triplets:
        possible_elements = set.intersection(*(node_to_elements[node] for node in triplet))
        for elem_id in possible_elements:
            if set(triplet).issubset(tetrahedral_elements[elem_id]):
                skin_to_tetra_mapping[triplet] = elem_id
                break  
    
    with open(output_filepath, "w") as output_file:
        output_file.write("*ELEMENTS_MAPPING, ELSET=EALL\n")
        for triplet, elem_id in skin_to_tetra_mapping.items():
            output_file.write(f"{elem_id}, {triplet[0]}, {triplet[1]}, {triplet[2]}\n")
    
    return skin_to_tetra_mapping

su2path = "TestCube/TestCube.su2"
# Extract all "fixed" nodes and write to Nfix1.nam
get_fixed_nodes(su2path)

#get_skin_nodes("SU2_MESH/CRMWS_WINGBOX.su2")

# Extract all "surface" nodes and write to Nsurface.nam
get_traction_nodes(su2path)


# Read all coordinates and connectivity matrix and write to mesh.nam
extract_su2_mesh_data_with_element_index(su2path)



#skin_triplets = extract_skin_node_triplets("SU2_MESH/CRMWS_WINGBOX.su2")


#output_filepath = "skin_tetra_elements.nam"
#find_all_tetrahedral_elements_for_skin_optimized("SU2_MESH/CRMWS_WINGBOX.su2", skin_triplets, output_filepath)
