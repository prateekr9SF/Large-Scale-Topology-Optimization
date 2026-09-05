#include "CalculiX.h"


void write_deformed_su2(int nk, double *v)
{
    FILE *fp;
    FILE *out;

    char su2file[512];
    char outfile[512];
    char line[MAXLINE];

    int nnode = 0;
    int node_section = 0;
    int node_counter = 0;

    /* Locate .su2 file before aeroelastic equilibrium */
    find_su2_file(su2file);

    /* Create .su2 file post aeroelastic equilibrium */
    sprintf(outfile, "deformed_%s", su2file);

    /* Try opening the above two files */
    fp = fopen(su2file, "r");
    if (fp == NULL) {
        printf("ERROR: Cannot open %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    out = fopen(outfile, "w");
    if (out == NULL) {
        printf("ERROR: Cannot open %s for writing\n", outfile);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        /* Node count consistency between CalTop and .su2 */
        if (strstr(line, "NPOIN=") != NULL) 
        {
            sscanf(line, "NPOIN= %d", &nnode);
            fprintf(out, "%s", line);

            if (nnode != nk) 
            {
                fprintf(stderr,
                    "\nERROR: Mesh mismatch detected!\n"
                    "       SU2 NPOIN : %d\n"
                    "       CalTop nk : %d\n"
                    "\n"
                    "The SU2 mesh and CalTop displacement vector do not match.\n"
                    "Cannot write deformed SU2 mesh.\n\n",
                nnode, nk);

                fclose(fp);
                fclose(out);

                exit(EXIT_FAILURE);
            }

            node_section = 1;
            node_counter = 0;
            continue;
        }

        if (node_section && node_counter < nnode) 
        {
            double x, y, z;
            int offset = 0;

            if (sscanf(line, "%lf %lf %lf %n",
                       &x, &y, &z, &offset) < 3) 
            {
                printf("ERROR: Could not parse SU2 node line %d\n",
                       node_counter + 1);
                fclose(fp);
                fclose(out);
                exit(EXIT_FAILURE);
            }
            
            /* Update nodal coordinates from aeroelastic solution */
            x += v[4 * node_counter + 1];
            y += v[4 * node_counter + 2];
            z += v[4 * node_counter + 3];

            fprintf(out, "%.15e %.15e %.15e", x, y, z);

            if (offset > 0 && line[offset] != '\0') 
            {
                fprintf(out, " %s", line + offset);
            }
            else 
            {
                fprintf(out, "\n");
            }

            node_counter++;

            if (node_counter == nnode) 
            {
                node_section = 0;
            }

            continue;
        }

        fprintf(out, "%s", line);
    }

    fclose(fp);
    fclose(out);

    printf("Updated SU2 mesh written: %s\n", outfile);
}

double tet_volume(double a[3],
                  double b[3],
                  double c[3],
                  double d[3])
{
    double bx = b[0] - a[0];
    double by = b[1] - a[1];
    double bz = b[2] - a[2];

    double cx = c[0] - a[0];
    double cy = c[1] - a[1];
    double cz = c[2] - a[2];

    double dx = d[0] - a[0];
    double dy = d[1] - a[1];
    double dz = d[2] - a[2];

    double cross_x = by * cz - bz * cy;
    double cross_y = bz * cx - bx * cz;
    double cross_z = bx * cy - by * cx;

    return (cross_x * dx +
            cross_y * dy +
            cross_z * dz) / 6.0;
}


void convert_volume_mesh(char *su2file)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];

    int nnode = 0;
    int nelem = 0;

    printf("Reading SU2 nodes ... ");

    fp = fopen(su2file, "r");

    if (fp == NULL) 
    {
        printf("ERROR: Cannot open SU2 file %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, "NPOIN=") != NULL) 
        {
            sscanf(line, "NPOIN= %d", &nnode);
            break;
        }
    }

    if (nnode <= 0) 
    {
        printf("ERROR: Could not read NPOIN from SU2 file\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    printf("%d nodes\n", nnode);

    double (*xyz)[3];

    xyz = malloc((size_t)nnode * sizeof(*xyz));

    if (xyz == NULL) 
    {
        printf("ERROR: Memory allocation failed for node coordinates\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < nnode; i++) 
    {
        if (fgets(line, MAXLINE, fp) == NULL) 
        {
            printf("ERROR: Unexpected end of file while reading nodes\n");
            free(xyz);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        sscanf(line,
               "%lf %lf %lf",
               &xyz[i][0],
               &xyz[i][1],
               &xyz[i][2]);
    }

    fclose(fp);

    out = fopen("mesh.nam", "w");

    if (out == NULL) 
    {
        printf("ERROR: Cannot open mesh.nam for writing\n");
        free(xyz);
        exit(EXIT_FAILURE);
    }

    fprintf(out, "** Generated from SU2 mesh\n");
    fprintf(out, "*NODE, NSET=NALL\n");

    for (int i = 0; i < nnode; i++) 
    {
        fprintf(out,
                "%d, %.10lf, %.10lf, %.10lf\n",
                i + 1,
                xyz[i][0],
                xyz[i][1],
                xyz[i][2]);
    }

    fp = fopen(su2file, "r");

    if (fp == NULL) 
    {
        printf("ERROR: Cannot reopen SU2 file %s\n", su2file);
        free(xyz);
        fclose(out);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, "NELEM=") != NULL) 
        {
            sscanf(line, "NELEM= %d", &nelem);
            break;
        }
    }

    if (nelem <= 0) 
    {
        printf("ERROR: Could not read NELEM from SU2 file\n");
        free(xyz);
        fclose(fp);
        fclose(out);
        exit(EXIT_FAILURE);
    }

    printf("Reading SU2 elements ... %d elements\n", nelem);

    fprintf(out, "\n*ELEMENT, TYPE=C3D4, ELSET=EALL\n");

    for (int e = 0; e < nelem; e++) 
    {
        int type;
        int n[4];
        int eid;

        if (fgets(line, MAXLINE, fp) == NULL) 
        {
            printf("ERROR: Unexpected end of file while reading elements\n");
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        if (sscanf(line,
                   "%d %d %d %d %d %d",
                   &type,
                   &n[0],
                   &n[1],
                   &n[2],
                   &n[3],
                   &eid) != 6) 
            {
                printf("ERROR: Could not parse element line %d\n", e + 1);
                free(xyz);
                fclose(fp);
                fclose(out);
                exit(EXIT_FAILURE);
            }

        if (type != 10) {
            printf("ERROR: Non-tetrahedral element detected. Element type = %d\n", type);
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        double vol = tet_volume(xyz[n[0]],xyz[n[1]],xyz[n[2]],xyz[n[3]]);

        if (fabs(vol) < 1.0e-30) 
        {
            printf("ERROR: Zero-volume tetrahedron detected at element %d\n", e + 1);
            free(xyz);
            fclose(fp);
            fclose(out);
            exit(EXIT_FAILURE);
        }

        if (vol < 0.0) 
        {
            int tmp = n[1];
            n[1] = n[2];
            n[2] = tmp;
        }

        fprintf(out,
                "%d, %d, %d, %d, %d\n",
                e + 1,
                n[0] + 1,
                n[1] + 1,
                n[2] + 1,
                n[3] + 1);
    }

    fclose(fp);
    fclose(out);

    free(xyz);

    printf("mesh.nam written\n");
}

void extract_marker(char *su2file,char *marker,char *outfile)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];
    char key[128];

    int nnode = 0;
    int active = 0;
    int found_marker = 0;
    int raw_node_count = 0;
    int unique_node_count = 0;

    sprintf(key, "MARKER_TAG= %s", marker);

    /*First pass: read NPOIN so we can allocate a node flag array. */

    fp = fopen(su2file, "r");

    if (fp == NULL) 
    {
        printf("ERROR: Cannot open SU2 file %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, "NPOIN=") != NULL) 
        {
            sscanf(line, "NPOIN= %d", &nnode);
            break;
        }
    }

    fclose(fp);

    if (nnode <= 0) 
    {
        printf("ERROR: Could not read NPOIN from SU2 file\n");
        exit(EXIT_FAILURE);
    }

    int *node_flag = calloc((size_t)nnode, sizeof(int));

    if (node_flag == NULL) 
    {
        printf("ERROR: Memory allocation failed in extract_marker\n");
        exit(EXIT_FAILURE);
    }

    /* Second pass: mark nodes belonging to the requested marker. */

    fp = fopen(su2file, "r");

    if (fp == NULL) 
    {
        printf("ERROR: Cannot open SU2 file %s\n", su2file);
        free(node_flag);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, key) != NULL) 
        {
            active = 1;
            found_marker = 1;
            continue;
        }

        if (active) 
        {
            if (strstr(line, "MARKER_TAG=") != NULL) 
            {
                break;
            }

            if (strstr(line, "MARKER_ELEMS=") != NULL) 
            {
                continue;
            }

            int type;
            int a, b, c;

            if (sscanf(line,
                       "%d %d %d %d",&type,&a,&b,&c) == 4) 
            {
                /*
                   SU2 surface triangle:
                       type a b c

                   Keep SU2 nodes internally as 0-based.
                */

                if (a >= 0 && a < nnode) node_flag[a] = 1;
                if (b >= 0 && b < nnode) node_flag[b] = 1;
                if (c >= 0 && c < nnode) node_flag[c] = 1;

                raw_node_count += 3;
            }
        }
    }

    fclose(fp);

    /*
       Write unique nodes as 1-based CalculiX IDs.
    */

    out = fopen(outfile, "w");

    if (out == NULL) 
    {
        printf("ERROR: Cannot open %s for writing\n", outfile);
        free(node_flag);
        exit(EXIT_FAILURE);
    }

    if (strcmp(marker, "fixed") == 0) 
    {
        fprintf(out, "*NSET, NSET=Nfix1\n");
    } 
    else 
    {
        fprintf(out, "*NSET, NSET=Nsurface\n");
    }

    for (int i = 0; i < nnode; i++) 
    {
        if (node_flag[i]) 
        {
            fprintf(out, "%d\n", i + 1);
            unique_node_count++;
        }
    }

    fclose(out);

    if (!found_marker) 
    {
        printf("WARNING: MARKER_TAG= %s not found. %s may be empty.\n",
               marker,
               outfile);
    } 
    else 
    {
        printf("%s written using MARKER_TAG= %s, raw node entries = %d, unique nodes = %d\n",
               outfile,
               marker,
               raw_node_count,
               unique_node_count);
    }

    free(node_flag);
}



void print_su2_markers(char *su2file)
{
    FILE *fp;
    char line[MAXLINE];
    int count = 0;

    fp = fopen(su2file,"r");

    if(fp == NULL)
    {
        printf("ERROR: Cannot open %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    printf("Detected SU2 boundary markers:\n");
    printf("--------------------------------\n");

    while(fgets(line,MAXLINE,fp))
    {
        if(strstr(line,"MARKER_TAG=") != NULL)
        {
            char tag[256];

            sscanf(line,"MARKER_TAG= %s",tag);

            printf("  %d) %s\n",count+1,tag);
            count++;
        }
    }

    fclose(fp);

    if(count == 0)
    {
        printf("No MARKER_TAG entries detected\n");
    }
    printf("--------------------------------\n\n");
}

char **get_su2_markers(char *su2file, int *marker_count)
{
    FILE *fp;
    char line[MAXLINE];

    char **markers = NULL;
    int count = 0;

    fp = fopen(su2file, "r");

    if (fp == NULL)
    {
        fprintf(stderr,
                "ERROR: Cannot open %s\n",
                su2file);

        exit(EXIT_FAILURE);
    }

    printf("Detected SU2 boundary markers:\n");
    printf("--------------------------------\n");

    while (fgets(line, MAXLINE, fp))
    {
        if (strstr(line, "MARKER_TAG=") != NULL)
        {
            char tag[256];

            if (sscanf(line, "MARKER_TAG= %255s", tag) == 1)
            {
                /*
                 * Increase the marker pointer array by one entry.
                 */
                char **temporary = realloc(
                    markers,
                    (size_t)(count + 1) * sizeof(char *)
                );

                if (temporary == NULL)
                {
                    fprintf(stderr,
                            "ERROR: Memory allocation failed while reading markers\n");

                    fclose(fp);

                    for (int i = 0; i < count; i++)
                    {
                        free(markers[i]);
                    }

                    free(markers);

                    exit(EXIT_FAILURE);
                }

                markers = temporary;

                /*
                 * Allocate memory for the marker name.
                 */
                markers[count] = malloc(strlen(tag) + 1);

                if (markers[count] == NULL)
                {
                    fprintf(stderr,
                            "ERROR: Memory allocation failed for marker name\n");

                    fclose(fp);

                    for (int i = 0; i < count; i++)
                    {
                        free(markers[i]);
                    }

                    free(markers);

                    exit(EXIT_FAILURE);
                }

                strcpy(markers[count], tag);

                printf("  %d) %s\n",
                       count + 1,
                       markers[count]);

                count++;
            }
        }
    }

    fclose(fp);

    if (count == 0)
    {
        printf("No MARKER_TAG entries detected\n");
    }

    printf("--------------------------------\n\n");

    /*
     * Return the number of markers through marker_count.
     */
    *marker_count = count;

    return markers;
}

void remove_existing_nam_files(void)
{
    DIR *dir;
    struct dirent *entry;


    int count = 0;

    dir = opendir(".");

    if(dir == NULL)
    {
        printf("ERROR: Cannot open current directory\n");
        exit(EXIT_FAILURE);
    }

    printf("Checking for old .nam files ...\n");

    while((entry = readdir(dir)) != NULL)
    {
        /*  Check filename extension*/
        if(strstr(entry->d_name,".nam") != NULL)
        {
            printf("Removing old file: %s\n",entry->d_name);

            if(remove(entry->d_name) != 0)
            {
                printf("WARNING: Could not remove %s\n",entry->d_name);
            }
            count++;
        }
    }
    closedir(dir);

    if(count == 0)
    {
        printf("No old .nam files found\n");
    }
    printf("\n");
}

void extract_skin_elements(char *su2file,char *markername,char *outfile)
{
    FILE *fp;
    FILE *out;

    char line[MAXLINE];
    char key[128];

    int nnode = 0;
    int nelem = 0;

    sprintf(key, "MARKER_TAG= %s", markername);

    /* First pass: read NPOIN and NELEM. */

    fp = fopen(su2file, "r");

    if (fp == NULL) 
    {
        printf("ERROR: Cannot open %s\n", su2file);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, "NPOIN=") != NULL) 
        {
            sscanf(line, "NPOIN= %d", &nnode);
        }

        if (strstr(line, "NELEM=") != NULL) 
        {
            sscanf(line, "NELEM= %d", &nelem);
        }
    }

    fclose(fp);

    if (nnode <= 0 || nelem <= 0) 
    {
        printf("ERROR: Could not read NPOIN or NELEM\n");
        exit(EXIT_FAILURE);
    }

    int *skin_nodes = calloc((size_t)nnode, sizeof(int));
    int *skin_elems = calloc((size_t)nelem, sizeof(int));

    if (skin_nodes == NULL || skin_elems == NULL) 
    {
        printf("ERROR: Memory allocation failed in extract_skin_elements\n");
        exit(EXIT_FAILURE);
    }

    /*
       Second pass: collect nodes on the chosen marker.

       SU2 surface triangle format:
           5 n1 n2 n3

       The node IDs are kept 0-based here.
    */

    fp = fopen(su2file, "r");

    int active = 0;
    int found_marker = 0;

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, key) != NULL) 
        {
            active = 1;
            found_marker = 1;
            continue;
        }

        if (active) 
        {
            if (strstr(line, "MARKER_TAG=") != NULL) 
            {
                break;
            }

            if (strstr(line, "MARKER_ELEMS=") != NULL) 
            {
                continue;
            }

            int type;
            int a, b, c;

            if (sscanf(line, "%d %d %d %d", &type, &a, &b, &c) == 4) 
            {
                skin_nodes[a] = 1;
                skin_nodes[b] = 1;
                skin_nodes[c] = 1;
            }
        }
    }

    fclose(fp);

    if (!found_marker) 
    {
        printf("WARNING: MARKER_TAG= %s not found. No skinElementList.nam written.\n",
               markername);

        free(skin_nodes);
        free(skin_elems);

        return;
    }

    /*
       Third pass: scan volume tetrahedra.

       SU2 tetrahedron format:
           10 n1 n2 n3 n4 elementID

       If a tetrahedron contains at least one skin node,
       it is marked as a skin/passive element.
    */

    fp = fopen(su2file, "r");

    int element_section = 0;

    while (fgets(line, MAXLINE, fp)) 
    {
        if (strstr(line, "NELEM=") != NULL) 
        {
            element_section = 1;
            continue;
        }

        if (!element_section) 
        {
            continue;
        }

        int type;
        int n1, n2, n3, n4;
        int eid;

        if (sscanf(line, "%d %d %d %d %d %d", &type, &n1, &n2, &n3, &n4, &eid) == 6) 
        {
            if (type != 10) 
            {
                continue;
            }

            /* if any node lies on current element, assign element to skin */
            if (skin_nodes[n1] || skin_nodes[n2] || skin_nodes[n3] || skin_nodes[n4]) 
            {
                skin_elems[eid] = 1;
            }
        }
    }

    fclose(fp);

    /* Write selected elements as 1-based IDs. */
    out = fopen(outfile, "w");

    if (out == NULL) 
    {
        printf("ERROR: Cannot open %s for writing\n", outfile);
        exit(EXIT_FAILURE);
    }

    int count = 0;

    for (int e = 0; e < nelem; e++) 
    {
        if (skin_elems[e]) 
        {
            fprintf(out, "%d\n", e + 1);
            count++;
        }
    }

    fclose(out);

    printf("%s written. Number of skin/passive elements = %d\n",
           outfile,
           count);

    free(skin_nodes);
    free(skin_elems);
}



void find_su2_file(char *filename)
{
    DIR *dir;
    struct dirent *entry;

    const char *suffix = "Solid.su2";
    size_t suffix_len = strlen(suffix);

    dir = opendir(".");

    if (dir == NULL)
    {
        printf("ERROR: Cannot open current directory\n");
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(dir)) != NULL)
    {
        size_t len = strlen(entry->d_name);

        /* Check whether the filename ends with "Solid.su2" */
        if (len >= suffix_len &&
            strcmp(entry->d_name + len - suffix_len, suffix) == 0)
        {

            printf("Detected SU2 mesh file: %s\n",filename);
            fflush(stdout);
            strcpy(filename, entry->d_name);

            closedir(dir);
            return;
        }
    }

    closedir(dir);

    printf("ERROR: No file ending with \"%s\" found in current directory\n",
           suffix);
    exit(EXIT_FAILURE);
}


void find_AD_su2_file(char *filename)
{
    DIR *dir;
    struct dirent *entry;

    const char *directory = "Adjoint";
    const char *prefix = "deformed_";
    size_t prefix_len = strlen(prefix);

    /* Open Adjoint directory */
    dir = opendir(directory);

    if (dir == NULL)
    {
        printf("ERROR: Cannot open directory \"%s\"\n", directory);
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    /* Search through files in Adjoint directory */
    while ((entry = readdir(dir)) != NULL)
    {
        size_t len = strlen(entry->d_name);

        /* Check whether the filename starts with "deformed_" */
        if (len >= prefix_len &&
            strncmp(entry->d_name, prefix, prefix_len) == 0)
        {
            /* Store complete path to detected file */
            sprintf(filename, "%s/%s", directory, entry->d_name);

            printf("Detected deformed mesh file: %s\n", filename);
            fflush(stdout);

            closedir(dir);
            return;
        }
    }

    /* No matching file was found */
    closedir(dir);

    printf("ERROR: No file starting with \"%s\" found in directory \"%s\"\n",
           prefix, directory);
    fflush(stdout);

    exit(EXIT_FAILURE);
}