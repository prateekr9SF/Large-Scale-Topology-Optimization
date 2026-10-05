#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @file tecplot_vtu.c
 * @brief Exports an elastic field solution to a VTU (VTK Unstructured Grid) file.
 *
 * This function writes nodal coordinates, connectivity information, and 
 * displacement field data to a VTU file (`elastic_Field.vtu`) for visualization 
 * in ParaView or other VTU-compatible software.
 *
 * @param nk     Number of nodes in the mesh.
 * @param ne     Number of elements in the mesh.
 * @param co     Array of nodal coordinates (size: 3 * nk).
 * @param kon    Element connectivity array, storing node indices per element.
 * @param ipkon  Index array mapping elements to the connectivity array.
 * @param lakon  Element Type
 * @param mi0    Element Index
 * @param v      Nodal displacement array (size: 3 * nk).
 * @param stx    Stress vector of elements
 * @param rhoPhy  Physical density of elements
 *
 * @details
 * The function writes the data in VTK Unstructured Grid (.vtu) format, including:
 * - **Nodal coordinates**: Each node's (x, y, z) position.
 * - **Element connectivity**: Defines elements by listing their node indices.
 * - **Offsets**: Specifies the end position of each element's connectivity in the list.
 * - **Cell types**: Assumes quadrilateral elements (VTK type 10).
 * - **Nodal displacements**: Exports displacement vectors at each node.
 *
 * @note 
 * - The function assumes 4-node quadrilateral elements.
 * - The generated file uses ASCII format for readability.
 * - The function does not return any value but writes the data to `elastic_Field.vtu`.
 * - If the file cannot be opened, an error is printed, and the function returns.
 *
 * @usage
 * ```c
 * int nk = 100; // Example number of nodes
 * int ne = 50;  // Example number of elements
 * double co[300];  // Node coordinates (3 per node)
 * int kon[200];  // Connectivity array (4 nodes per element)
 * int ipkon[50]; // Index mapping elements to connectivity array
 * char lakon = 'C3D4'; // Index mapping elements to connectivity array
 * int mi0[50];			// Element Index Array
 * double v[300];  // Nodal displacements (3 per node)
 * double stx[300];		// Stress per element, 6 stress for each element
 * double rhoPhys[50];  // Element Density, 1 per element
 * 
 * tecplot_vtu(nk, ne, co, kon, ipkon, lakon, mi0, v, stx, rhoPhys);
 * ```
 */

void tecplot_vtu(int nk, int ne, double *co, int *kon, int *ipkon, char *lakon, int mi0, double *v, double *stx, double *rhoPhy)    
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

    // Write displaced nodal coordinates
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    
    for (int node = 0; node < nk; node++) 
    {
        fprintf(fp, "        %.5f %.5f %.5f\n",
                co[3 * node]     + v[4 * node +1],
                co[3 * node + 1] + v[4 * node + 2],
                co[3 * node + 2] + v[4 * node +3]);
    }

    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // Write connectivity
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk= &lakon[8*ielem];
        int nope = 4;
        if (lk[3]=='1' && lk[4]=='0') nope = 10;
        for (int j = 0; j < nope; j++) {  // Node # based on nope (C3D4or C3D10)
            fprintf(fp, " %d", kon[ipkon[ielem] + j]-1);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Write offsets
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    
    int offset = 0;
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk = &lakon[8*ielem];
        offset += (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        fprintf(fp, " %d\n", offset);
    }
    fprintf(fp, "        </DataArray>\n");

    // Write cell types
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int ielem = 0; ielem < ne; ielem++) {
        char *lk = &lakon[8*ielem];
        int vtktype = (lk[3]=='1' && lk[4]=='0') ? 24 : 10;
        fprintf(fp, " %d\n", vtktype);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // Write nodal displacements
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");

    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    for (int node = 0; node < nk; node++) 
    {
        fprintf(fp, "        %.8f %.8f %.8f\n", v[4*node+1], v[4*node+2], v[4* node+3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // Write Cell Stress and Density data
    // Set default scalar to VonMises (nice for coloring in ParaView)
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // --- Von Mises scalar per cell ---
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises(Pa)\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int cell = 0; cell < ne; cell++) {
        char *lk = &lakon[8*cell];
        int ngp = (lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0, syy=0, szz=0;
        double txy=0, txz=0, tyz=0;

        for (int gp = 0; gp < ngp; gp++) {
            sxx += stx[6*mi0*cell + 6*gp    ];
            syy += stx[6*mi0*cell + 6*gp + 1];
            szz += stx[6*mi0*cell + 6*gp + 2];
            txy += stx[6*mi0*cell + 6*gp + 3];
            txz += stx[6*mi0*cell + 6*gp + 4];
            tyz += stx[6*mi0*cell + 6*gp + 5];
        }
        sxx/=ngp; syy/=ngp; szz/=ngp; txy/=ngp; tyz/=ngp; txz/=ngp;

        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
                + 3.0 * (txy*txy + txz*txz + tyz*tyz)
        );

        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");

    // Sigma_xx
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_xx(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0.0;

        for(int gp=0; gp<ngp; gp++)
            sxx += stx[6*mi0*cell + 6*gp];

        sxx /= ngp;

        fprintf(fp," %.8f\n",sxx);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_yy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_yy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double syy=0.0;

        for(int gp=0; gp<ngp; gp++)
            syy += stx[6*mi0*cell + 6*gp+1];

        syy /= ngp;

        fprintf(fp," %.8f\n",syy);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_zz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_zz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double szz=0.0;

        for(int gp=0; gp<ngp; gp++)
            szz += stx[6*mi0*cell + 6*gp + 2];

        szz /= ngp;

        fprintf(fp," %.8f\n",szz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xy=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xy += stx[6*mi0*cell + 6*gp + 3];

        tau_xy /= ngp;

        fprintf(fp," %.8f\n",tau_xy);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xz += stx[6*mi0*cell + 6*gp + 4];

        tau_xz /= ngp;

        fprintf(fp," %.8f\n",tau_xz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_yz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_yz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_yz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_yz += stx[6*mi0*cell + 6*gp + 5];

        tau_yz /= ngp;

        fprintf(fp," %.8f\n",tau_yz);
    }

    fprintf(fp,"        </DataArray>\n");

    // --- Density per cell ---
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    
    for (int cell = 0; cell < ne; cell++) 
    {
        fprintf(fp, "      %.8f\n", rhoPhy[cell]);
    }

    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");


    // Close XML tags
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}

void tecplot_vtu_passive(int nk, int ne,
                         double *co, int *kon, int *ipkon,
                         char *lakon, int mi0,
                         double *v, double *stx, double *rhoPhy,
                         int *passiveIDs, int numPassive, const char *filename)
{

    FILE *fp = fopen(filename, "w");

    if (fp == NULL)
    {
        perror(filename);
        return;
    }

    // --- VTU header ---
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    // Only passive elements are written as cells; keep all nodes
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, numPassive);

    // --- Write displaced nodal coordinates ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n",
                co[3 * node    ] + v[4 * node +1],
                co[3 * node + 1] + v[4 * node +2],
                co[3 * node + 2] + v[4 * node +3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // --- Cells: only passive elements ---
    fprintf(fp, "      <Cells>\n");

    // Connectivity
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;  // convert 1-based element ID to 0-based index
        char *lk = &lakon[8*e];
        int nope = (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        for (int j = 0; j < nope; j++) {  // C3D4 or C3D10
            fprintf(fp, " %d", kon[ipkon[e] + j] - 1);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Offsets (still 4 nodes per cell)
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    int offset = 0;
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        char *lk = &lakon[8*e];
        offset += (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        fprintf(fp, " %d\n", offset);
    }
    fprintf(fp, "        </DataArray>\n");

    // Types (VTK_TETRA = 10 for 4-node tets)
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        char *lk = &lakon[8*e];
        int vtktype = (lk[3]=='1' && lk[4]=='0') ? 24 : 10;
        fprintf(fp, " %d\n", vtktype);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // --- PointData: Displacement at all nodes ---
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        // CalculiX vold-style: [0]=dummy, [1]=ux, [2]=uy, [3]=uz
        fprintf(fp, "        %.8f %.8f %.8f\n",
                v[4 * node + 1],
                v[4 * node + 2],
                v[4 * node + 3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // --- CellData: VonMises, Stress, Density only for passive elements ---
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // Von Mises per passive cell
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        char *lk = &lakon[8*e];
        int ngp = (lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0, syy=0, szz=0;
        double txy=0, txz=0, tyz=0;


        for (int gp = 0; gp < ngp; gp++) {
            sxx += stx[6*mi0*e + 6*gp    ];
            syy += stx[6*mi0*e + 6*gp + 1];
            szz += stx[6*mi0*e + 6*gp + 2];
            txy += stx[6*mi0*e + 6*gp + 3];
            txz += stx[6*mi0*e + 6*gp + 4];
            tyz += stx[6*mi0*e + 6*gp + 5];
        }
        sxx/=ngp; syy/=ngp; szz/=ngp; txy/=ngp; tyz/=ngp; txz/=ngp;
        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
                + 3.0 * (txy*txy + txz*txz + tyz*tyz)
        );
        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");


    // Sigma_xx
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_xx(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0.0;

        for(int gp=0; gp<ngp; gp++)
            sxx += stx[6*mi0*cell + 6*gp];

        sxx /= ngp;

        fprintf(fp," %.8f\n",sxx);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_yy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_yy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double syy=0.0;

        for(int gp=0; gp<ngp; gp++)
            syy += stx[6*mi0*cell + 6*gp+1];

        syy /= ngp;

        fprintf(fp," %.8f\n",syy);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_zz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_zz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double szz=0.0;

        for(int gp=0; gp<ngp; gp++)
            szz += stx[6*mi0*cell + 6*gp + 2];

        szz /= ngp;

        fprintf(fp," %.8f\n",szz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xy=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xy += stx[6*mi0*cell + 6*gp + 3];

        tau_xy /= ngp;

        fprintf(fp," %.8f\n",tau_xy);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xz += stx[6*mi0*cell + 6*gp + 4];

        tau_xz /= ngp;

        fprintf(fp," %.8f\n",tau_xz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_yz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_yz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int k=0; k<numPassive; k++)
    {
        int cell = passiveIDs[k]-1;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_yz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_yz += stx[6*mi0*cell + 6*gp + 5];

        tau_yz /= ngp;

        fprintf(fp," %.8f\n",tau_yz);
    }

    fprintf(fp,"        </DataArray>\n");


    // Density for passive cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int k = 0; k < numPassive; k++) {
        int e = passiveIDs[k] - 1;
        fprintf(fp, "      %.8f\n", rhoPhy[e]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");

    // --- Close ---
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
}


void tecplot_vtu_active(int nk, int ne,
                        double *co, int *kon, int *ipkon,
                        char *lakon, int mi0,
                        double *v, double *stx, double *rhoPhy,
                        int *passiveIDs, int numPassive)
{
    // Build a mask of which elements are passive
    int *isPassive = (int*)calloc(ne, sizeof(int));
    if (!isPassive) {
        perror("Error allocating isPassive");
        return;
    }

    for (int i = 0; i < numPassive; i++) {
        int eID = passiveIDs[i];   // 1-based element ID from skinElementList.nam
        int e   = eID - 1;         // convert to 0-based index
        if (e >= 0 && e < ne) {
            isPassive[e] = 1;
        }
    }

    // Count active elements
    int numActive = 0;
    for (int e = 0; e < ne; e++) {
        if (!isPassive[e]) numActive++;
    }

    if (numActive == 0) {
        fprintf(stderr, "Warning: no active elements found (all are passive?). Skipping VTU write.\n");
        free(isPassive);
        return;
    }

    FILE *fp = fopen("elastic_Field_active.vtu", "w");
    if (fp == NULL) {
        perror("Error opening elastic_Field_active.vtu");
        free(isPassive);
        return;
    }

    // --- VTU header ---
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid>\n");
    // Only active elements as cells; keep all nodes
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nk, numActive);

    // --- write dispaced nodal coordinates ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        fprintf(fp, "        %.5f %.5f %.5f\n",
                co[3 * node    ] + v[4 * node +1],
                co[3 * node + 1] + v[4* node +2],
                co[3 * node + 2] + v[4 * node +3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");

    // --- Cells: only active elements ---
    fprintf(fp, "      <Cells>\n");

    // Connectivity
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        char *lk = &lakon[8*e];
        int nope = (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        for (int j = 0; j < nope; j++) {  // C3D4 or C3D10
            fprintf(fp, " %d", kon[ipkon[e] + j] - 1);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // Offsets: 4 nodes per active cell
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    int offset = 0;
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        char *lk = &lakon[8*e];
        offset += (lk[3]=='1' && lk[4]=='0') ? 10 : 4;
        fprintf(fp, " %d\n", offset);
    }
    fprintf(fp, "        </DataArray>\n");

    // Types: VTK_TET = 10
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        char *lk = &lakon[8*e];
        int vtktype = (lk[3]=='1' && lk[4]=='0') ? 24 : 10;
        fprintf(fp, " %d\n", vtktype);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");

    // --- PointData: Displacement at all nodes ---
    fprintf(fp, "      <PointData Scalars=\"Displacement\">\n");
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    for (int node = 0; node < nk; node++) {
        // CalculiX vold-style: [0]=dummy, [1]=ux, [2]=uy, [3]=uz
        fprintf(fp, "        %.8f %.8f %.8f\n",
                v[4 * node + 1],
                v[4 * node + 2],
                v[4 * node + 3]);
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </PointData>\n");

    // --- CellData: VonMises, Stress, Density only for active elements ---
    fprintf(fp, "      <CellData Scalars=\"VonMises\">\n");

    // Von Mises per active cell
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"VonMises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;

        char *lk = &lakon[8*e];
        int ngp = (lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0, syy=0, szz=0;
        double txy=0, txz=0, tyz=0;

        for (int gp = 0; gp < ngp; gp++) {
            sxx += stx[6*mi0*e + 6*gp    ];
            syy += stx[6*mi0*e + 6*gp + 1];
            szz += stx[6*mi0*e + 6*gp + 2];
            txy += stx[6*mi0*e + 6*gp + 3];
            txz += stx[6*mi0*e + 6*gp + 4];
            tyz += stx[6*mi0*e + 6*gp + 5];
        }

        sxx/=ngp; syy/=ngp; szz/=ngp; txy/=ngp; tyz/=ngp; txz/=ngp;
        double vm = sqrt(
            0.5 * ( (sxx - syy)*(sxx - syy)
                  + (syy - szz)*(syy - szz)
                  + (szz - sxx)*(szz - sxx) )
                + 3.0 * (txy*txy + txz*txz + tyz*tyz)
        );
        fprintf(fp, "      %.8f\n", vm);
    }
    fprintf(fp, "        </DataArray>\n");

    // Sigma_xx
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_xx(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

   for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;

        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double sxx=0.0;

        for(int gp=0; gp<ngp; gp++)
            sxx += stx[6*mi0*cell + 6*gp];

        sxx /= ngp;

        fprintf(fp," %.8f\n",sxx);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_yy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_yy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double syy=0.0;

        for(int gp=0; gp<ngp; gp++)
            syy += stx[6*mi0*cell + 6*gp+1];

        syy /= ngp;

        fprintf(fp," %.8f\n",syy);
    }

    fprintf(fp,"        </DataArray>\n");

    // Sigma_zz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Sigma_zz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double szz=0.0;

        for(int gp=0; gp<ngp; gp++)
            szz += stx[6*mi0*cell + 6*gp + 2];

        szz /= ngp;

        fprintf(fp," %.8f\n",szz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xy
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xy(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;

        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xy=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xy += stx[6*mi0*cell + 6*gp + 3];

        tau_xy /= ngp;

        fprintf(fp," %.8f\n",tau_xy);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_xz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_xz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

      for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_xz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_xz += stx[6*mi0*cell + 6*gp + 4];

        tau_xz /= ngp;

        fprintf(fp," %.8f\n",tau_xz);
    }

    fprintf(fp,"        </DataArray>\n");

    // tau_yz
    fprintf(fp,
    "        <DataArray type=\"Float64\" Name=\"Tau_yz(Pa)\" "
    "NumberOfComponents=\"1\" format=\"ascii\">\n");

    for(int cell=0; cell<ne; cell++)
    {
        if(isPassive[cell]) continue;
        char *lk=&lakon[8*cell];
        int ngp=(lk[3]=='1' && lk[4]=='0') ? 4 : 1;

        double tau_yz=0.0;

        for(int gp=0; gp<ngp; gp++)
            tau_yz += stx[6*mi0*cell + 6*gp + 5];

        tau_yz /= ngp;

        fprintf(fp," %.8f\n",tau_yz);
    }

    fprintf(fp,"        </DataArray>\n");

    // Density for active cells
    fprintf(fp, "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n");
    for (int e = 0; e < ne; e++) {
        if (isPassive[e]) continue;
        fprintf(fp, "      %.8f\n", rhoPhy[e]);
    }
    fprintf(fp, "        </DataArray>\n");

    fprintf(fp, "      </CellData>\n");

    // --- Close ---
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </UnstructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    fclose(fp);
    free(isPassive);
}
