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
                    formatted_element = ",".join(elements[1:-1])  # Exclude first and last column
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
            formatted_line = f"{idx},{coords[0]},{coords[1]},{coords[2]}"
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
    offsets the node indices by 1, and writes them to Nfix1.nam with one entry per line,
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
        # Identify the marker for "root" (lowercase as found in the file)
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


# Function to extract node indices associated with the "root" marker without sorting and with a comma per line
def get_root_nodes(su2_filepath, output_filepath="Nfix1.nam"):
    """
    Extracts node indices associated with the marker 'root' from a SU2 file,
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
        # Identify the marker for "root" (lowercase as found in the file)
        if "MARKER_TAG= root" in line:
            root_section = True
            continue  # Move to the next line where elements start
        
        # If in ROOT section, read node indices from element connectivity
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
        return "Error: No nodes found for marker 'root'. Check SU2 file format."

    # Write the extracted ROOT node indices to the output file with one entry per line, comma-separated
    with open(output_filepath, "w") as output_file:
        output_file.write("*NSET, NSET=Nfix1\n")
        
        # Write nodes one per line, with a comma after each entry
        for node in root_nodes:
            output_file.write(f"{node},\n")

    #return f"Successfully extracted {len(root_nodes)} ROOT nodes with offset, written one per line with comma (original order preserved), saved to {output_filepath}"

# Run the function to extract ROOT node indices without sorting, with a comma per line
get_root_nodes("SU2_MESH/CRMWS_WINGBOX.su2")
get_skin_nodes("SU2_MESH/CRMWS_WINGBOX.su2")
extract_su2_mesh_data_with_element_index("SU2_MESH/CRMWS_WINGBOX.su2")

