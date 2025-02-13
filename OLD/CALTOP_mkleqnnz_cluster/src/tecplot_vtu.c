#include <stdio.h>
#include <stdlib.h>

void tecplot_vtu(int nk, int ne, double *co, int *kon, int *ipkon, double *v) 
    
    {
    FILE *fp = fopen("elastic_Field.vtu", "w");
    if (fp == NULL) {
        perror("Error opening file");
        return;
    }

    // Write VTU XML header
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, ne);

    // Write nodal coordinates
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n", co[3 * node], co[3 * node + 1], co[3 * node + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // Write connectivity
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        for (int j = 0; j < 4; j++) {  // Assuming quadrilateral elements
            fprintf(fp, " %d", kon[ipkon[ielem] + j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Write offsets
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        fprintf(fp, " %d\n", (ielem + 1) * 4);
    }
    fprintf(fp, "        </DataArray>\n");

    // Write cell types
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        fprintf(fp, " 10\n");  // Assuming VTU type 10 for quadrilaterals
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // Write nodal displacements
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.8f %.8f %.8f\n", v[3 * node], v[3 * node + 1], v[3 * node + 2]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // Close XML tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}
