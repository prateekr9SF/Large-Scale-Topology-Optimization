def extract_su2_coordinates(su2_filepath, output_filepath="mesh.nam"):
    """
    Extracts node coordinates from a SU2 file and writes them to mesh.nam.
    
    Parameters:
    su2_filepath (str): Path to the input SU2 file.
    output_filepath (str): Path to the output text file (default: mesh.nam).
    """
    with open(su2_filepath, "r") as file:
        lines = file.readlines()
    
    # Find the section where node coordinates start
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
                break
    
    # Write extracted coordinates to the output file
    with open(output_filepath, "w") as output_file:
        for node in node_data:
            output_file.write(node + "\n")
    
    print(f"Extracted {len(node_data)} nodes and saved to {output_filepath}")

# Example usage
extract_su2_coordinates("/mnt/data/CRMWS_WINGBOX.su2")

