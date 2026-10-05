#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <CalFSI.h>

/* ============================================================
   Data structures
   ============================================================ */

typedef struct
{
    double x;
    double y;
    double z;
} HPNode;

typedef struct
{
    int node[4];
    int id;             /* CalculiX/CalTop 1-based element ID */

    double xmin, xmax;
    double ymin, ymax;

    double xc;
    double yc;
    double zc;
} HPTet;


/* ============================================================
   Comparison function for qsort()
   ============================================================ */

static int compare_ints(const void *a, const void *b)
{
    int ia = *(const int *)a;
    int ib = *(const int *)b;

    if (ia < ib) return -1;
    if (ia > ib) return 1;
    return 0;
}


/* ============================================================
   Comparison function for doubles
   ============================================================ */

static int compare_doubles(const void *a, const void *b)
{
    double da = *(const double *)a;
    double db = *(const double *)b;

    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}


/* ============================================================
   Median of array
   ============================================================ */

static double median_double(double *array, int n)
{
    double *tmp;
    double value;

    if (n <= 0)
        return 0.0;

    tmp = (double *)malloc((size_t)n * sizeof(double));

    if (!tmp)
    {
        fprintf(stderr, "ERROR: malloc failed in median_double\n");
        exit(EXIT_FAILURE);
    }

    memcpy(tmp, array, (size_t)n * sizeof(double));

    qsort(tmp, n, sizeof(double), compare_doubles);

    if (n % 2 == 0)
        value = 0.5 * (tmp[n/2 - 1] + tmp[n/2]);
    else
        value = tmp[n/2];

    free(tmp);

    return value;
}


/* ============================================================
   Percentile of sorted/unsorted double array

   p = 0.005 -> 0.5 percentile
   p = 0.995 -> 99.5 percentile
   ============================================================ */

static double percentile_double(double *array, int n, double p)
{
    double *tmp;
    double position;
    int i0, i1;
    double fraction;
    double value;

    if (n <= 0)
        return 0.0;

    tmp = (double *)malloc((size_t)n * sizeof(double));

    if (!tmp)
    {
        fprintf(stderr,
                "ERROR: malloc failed in percentile_double\n");
        exit(EXIT_FAILURE);
    }

    memcpy(tmp, array, (size_t)n * sizeof(double));

    qsort(tmp, n, sizeof(double), compare_doubles);

    position = p * (double)(n - 1);

    i0 = (int)floor(position);
    i1 = (int)ceil(position);

    fraction = position - (double)i0;

    value = tmp[i0] * (1.0 - fraction)
          + tmp[i1] * fraction;

    free(tmp);

    return value;
}


/* ============================================================
   Extract slat hard-point elements from SU2 volume mesh

   Hard points:
       eta = 0.100
             0.275
             0.450
             0.625
             0.800

   Nominal spanwise width:
       2 mm = +/- 1 mm

   Chordwise extent:
       local LE -> 25% local chord

   Neighborhood:
       one-ring node-connected tetrahedra

   ============================================================ */

void extract_slat_hardpoints(const char *su2file,
                             const char *output_file)
{
    FILE *fp;
    FILE *out;

    char line[4096];

    int nelem = 0;
    int npoin = 0;

    HPNode *nodes = NULL;
    HPTet  *tets  = NULL;

    int ntet = 0;

    /* --------------------------------------------------------
       Hard-point parameters
       -------------------------------------------------------- */

    const int nhp = 5;

    const double eta[5] =
    {
        0.100,
        0.275,
        0.450,
        0.625,
        0.800
    };

    const double hardpoint_width = 0.002;
    const double half_width = 0.5 * hardpoint_width;

    const double chord_fraction = 0.25;


    /* --------------------------------------------------------
       Open SU2 mesh
       -------------------------------------------------------- */

    fp = fopen(su2file, "r");

    if (!fp)
    {
        perror("ERROR opening SU2 mesh");
        exit(EXIT_FAILURE);
    }


    /* ========================================================
       PASS 1: Read volume tetrahedra
       ======================================================== */

    while (fgets(line, sizeof(line), fp))
    {
        if (strncmp(line, "NELEM=", 6) == 0 ||
            strncmp(line, "NELEM =", 7) == 0)
        {
            char *eq = strchr(line, '=');

            if (!eq)
                continue;

            nelem = atoi(eq + 1);

            tets = (HPTet *)malloc(
                (size_t)nelem * sizeof(HPTet));

            if (!tets)
            {
                fprintf(stderr,
                        "ERROR allocating tetrahedra\n");
                fclose(fp);
                exit(EXIT_FAILURE);
            }

            for (int e = 0; e < nelem; e++)
            {
                int type;
                int n1, n2, n3, n4;
                int eid;

                if (!fgets(line, sizeof(line), fp))
                    break;

                /*
                   SU2 tetrahedral format:

                   10 n1 n2 n3 n4 element_id

                   Type 10 = tetrahedron
                */

                int nread = sscanf(
                    line,
                    "%d %d %d %d %d %d",
                    &type,
                    &n1,
                    &n2,
                    &n3,
                    &n4,
                    &eid
                );

                if (type == 10 && nread >= 5)
                {
                    tets[ntet].node[0] = n1;
                    tets[ntet].node[1] = n2;
                    tets[ntet].node[2] = n3;
                    tets[ntet].node[3] = n4;

                    /*
                       If SU2 element ID exists, it is normally
                       zero-based.

                       CalculiX/CalTop uses one-based IDs.
                    */

                    if (nread == 6)
                        tets[ntet].id = eid + 1;
                    else
                        tets[ntet].id = e + 1;

                    ntet++;
                }
            }

            break;
        }
    }


    /* ========================================================
       PASS 2: Read nodes
       ======================================================== */

    rewind(fp);

    while (fgets(line, sizeof(line), fp))
    {
        if (strncmp(line, "NPOIN=", 6) == 0 ||
            strncmp(line, "NPOIN =", 7) == 0)
        {
            char *eq = strchr(line, '=');

            if (!eq)
                continue;

            npoin = atoi(eq + 1);

            nodes = (HPNode *)malloc(
                (size_t)npoin * sizeof(HPNode));

            if (!nodes)
            {
                fprintf(stderr,
                        "ERROR allocating mesh nodes\n");

                free(tets);
                fclose(fp);

                exit(EXIT_FAILURE);
            }

            for (int i = 0; i < npoin; i++)
            {
                double x, y, z;
                int node_id;

                if (!fgets(line, sizeof(line), fp))
                    break;

                /*
                   SU2 point format may contain:

                   x y z node_id
                */

                int nread = sscanf(
                    line,
                    "%lf %lf %lf %d",
                    &x,
                    &y,
                    &z,
                    &node_id
                );

                if (nread >= 3)
                {
                    nodes[i].x = x;
                    nodes[i].y = y;
                    nodes[i].z = z;
                }
            }

            break;
        }
    }

    fclose(fp);


    if (!nodes || !tets || ntet == 0)
    {
        fprintf(stderr,
                "ERROR: Failed to read SU2 volume mesh\n");

        free(nodes);
        free(tets);

        exit(EXIT_FAILURE);
    }


    printf("\nGenerating slat hard points\n");
    printf("------------------------------------------\n");
    printf("Nodes        : %d\n", npoin);
    printf("Tetrahedra   : %d\n", ntet);


    /* ========================================================
       Compute tetrahedron geometry
       ======================================================== */

    for (int e = 0; e < ntet; e++)
    {
        int n0 = tets[e].node[0];
        int n1 = tets[e].node[1];
        int n2 = tets[e].node[2];
        int n3 = tets[e].node[3];

        HPNode *p0 = &nodes[n0];
        HPNode *p1 = &nodes[n1];
        HPNode *p2 = &nodes[n2];
        HPNode *p3 = &nodes[n3];

        tets[e].xc =
            0.25 * (p0->x + p1->x + p2->x + p3->x);

        tets[e].yc =
            0.25 * (p0->y + p1->y + p2->y + p3->y);

        tets[e].zc =
            0.25 * (p0->z + p1->z + p2->z + p3->z);


        tets[e].xmin = p0->x;
        tets[e].xmax = p0->x;

        tets[e].ymin = p0->y;
        tets[e].ymax = p0->y;


        HPNode *pv[4] =
        {
            p0, p1, p2, p3
        };

        for (int j = 1; j < 4; j++)
        {
            if (pv[j]->x < tets[e].xmin)
                tets[e].xmin = pv[j]->x;

            if (pv[j]->x > tets[e].xmax)
                tets[e].xmax = pv[j]->x;

            if (pv[j]->y < tets[e].ymin)
                tets[e].ymin = pv[j]->y;

            if (pv[j]->y > tets[e].ymax)
                tets[e].ymax = pv[j]->y;
        }
    }


    /* ========================================================
       Determine total span
       ======================================================== */

    double y_min = nodes[0].y;
    double y_max = nodes[0].y;

    for (int i = 1; i < npoin; i++)
    {
        if (nodes[i].y < y_min)
            y_min = nodes[i].y;

        if (nodes[i].y > y_max)
            y_max = nodes[i].y;
    }

    double span = y_max - y_min;

    printf("Span         : %.10f\n", span);
    printf("y_min        : %.10f\n", y_min);
    printf("y_max        : %.10f\n\n", y_max);


    /* ========================================================
       Build node -> element adjacency

       First count how many tetrahedra touch each node.
       ======================================================== */

    int *node_count = (int *)calloc(
        (size_t)npoin,
        sizeof(int));

    if (!node_count)
    {
        fprintf(stderr,
                "ERROR allocating node_count\n");

        free(nodes);
        free(tets);

        exit(EXIT_FAILURE);
    }


    for (int e = 0; e < ntet; e++)
    {
        for (int j = 0; j < 4; j++)
        {
            int n = tets[e].node[j];

            if (n >= 0 && n < npoin)
                node_count[n]++;
        }
    }


    int **node_elem =
        (int **)malloc((size_t)npoin * sizeof(int *));

    int *node_fill =
        (int *)calloc((size_t)npoin, sizeof(int));


    if (!node_elem || !node_fill)
    {
        fprintf(stderr,
                "ERROR allocating node adjacency\n");

        free(node_count);
        free(node_elem);
        free(node_fill);
        free(nodes);
        free(tets);

        exit(EXIT_FAILURE);
    }


    for (int i = 0; i < npoin; i++)
    {
        if (node_count[i] > 0)
        {
            node_elem[i] = (int *)malloc(
                (size_t)node_count[i] * sizeof(int));
        }
        else
        {
            node_elem[i] = NULL;
        }
    }


    for (int e = 0; e < ntet; e++)
    {
        for (int j = 0; j < 4; j++)
        {
            int n = tets[e].node[j];

            node_elem[n][node_fill[n]++] = e;
        }
    }


    /* ========================================================
       Global selected-element flags
       ======================================================== */

    unsigned char *selected_global =
        (unsigned char *)calloc(
            (size_t)ntet,
            sizeof(unsigned char));


    if (!selected_global)
    {
        fprintf(stderr,
                "ERROR allocating selected_global\n");

        exit(EXIT_FAILURE);
    }


    /* ========================================================
       Loop over five hard-point stations
       ======================================================== */

    for (int hp = 0; hp < nhp; hp++)
    {
        double y_hp =
            y_min + eta[hp] * span;


        /* ----------------------------------------------------
           Determine local chord.

           Same approach used during mesh inspection:

           collect nodes within an adaptive spanwise window.
           ---------------------------------------------------- */

        double dy = fmax(
            0.03,
            0.0025 * span
        );

        int local_count = 0;


        while (1)
        {
            local_count = 0;

            for (int i = 0; i < npoin; i++)
            {
                if (fabs(nodes[i].y - y_hp) <= dy)
                    local_count++;
            }

            if (local_count >= 250)
                break;

            dy *= 1.5;
        }


        double *local_x =
            (double *)malloc(
                (size_t)local_count * sizeof(double));


        if (!local_x)
        {
            fprintf(stderr,
                    "ERROR allocating local_x\n");

            exit(EXIT_FAILURE);
        }


        int k = 0;

        for (int i = 0; i < npoin; i++)
        {
            if (fabs(nodes[i].y - y_hp) <= dy)
            {
                local_x[k++] = nodes[i].x;
            }
        }


        /*
           Robust LE / TE estimates.

           0.5 percentile and 99.5 percentile were used
           rather than absolute min/max to suppress isolated
           mesh-coordinate outliers.
        */

        double x_le =
            percentile_double(
                local_x,
                local_count,
                0.005
            );

        double x_te =
            percentile_double(
                local_x,
                local_count,
                0.995
            );

        double chord = x_te - x_le;

        double x25 =
            x_le + chord_fraction * chord;


        free(local_x);


        /* ----------------------------------------------------
           Find seed elements intersecting nominal 2-mm band.
           ---------------------------------------------------- */

        unsigned char *seed =
            (unsigned char *)calloc(
                (size_t)ntet,
                sizeof(unsigned char));

        unsigned char *neighbor =
            (unsigned char *)calloc(
                (size_t)ntet,
                sizeof(unsigned char));


        if (!seed || !neighbor)
        {
            fprintf(stderr,
                    "ERROR allocating hard-point flags\n");

            exit(EXIT_FAILURE);
        }


        int seed_count = 0;


        for (int e = 0; e < ntet; e++)
        {
            int y_intersects =
                (tets[e].ymax >= y_hp - half_width) &&
                (tets[e].ymin <= y_hp + half_width);

            int x_intersects =
                (tets[e].xmax >= x_le) &&
                (tets[e].xmin <= x25);


            if (y_intersects && x_intersects)
            {
                seed[e] = 1;
                neighbor[e] = 1;

                seed_count++;
            }
        }


        /* ----------------------------------------------------
           Determine median spanwise size of seed elements.
           ---------------------------------------------------- */

        double *seed_size = NULL;

        double local_h = 0.0;


        if (seed_count > 0)
        {
            seed_size = (double *)malloc(
                (size_t)seed_count * sizeof(double));


            if (!seed_size)
            {
                fprintf(stderr,
                        "ERROR allocating seed_size\n");

                exit(EXIT_FAILURE);
            }


            k = 0;

            for (int e = 0; e < ntet; e++)
            {
                if (seed[e])
                {
                    seed_size[k++] =
                        tets[e].ymax -
                        tets[e].ymin;
                }
            }


            local_h =
                median_double(
                    seed_size,
                    seed_count
                );


            free(seed_size);
        }


        /*
           Neighborhood tolerance.

           Same logic used in the tested Python version.
        */

        double tol_y =
            half_width +
            fmax(local_h, 0.002);


        /* ----------------------------------------------------
           One-ring node-connected neighborhood search
           ---------------------------------------------------- */

        for (int e = 0; e < ntet; e++)
        {
            if (!seed[e])
                continue;


            for (int j = 0; j < 4; j++)
            {
                int n = tets[e].node[j];


                for (int q = 0;
                     q < node_count[n];
                     q++)
                {
                    int e2 =
                        node_elem[n][q];

                    neighbor[e2] = 1;
                }
            }
        }


        /* ----------------------------------------------------
           Geometrically restrict neighborhood.

           The centroid must remain:

             LE <= x <= 25%c

           and close to the hard-point station.
           ---------------------------------------------------- */

        int hp_selected = 0;


        for (int e = 0; e < ntet; e++)
        {
            if (!neighbor[e])
                continue;


            if (fabs(tets[e].yc - y_hp) > tol_y)
                continue;


            if (tets[e].xc < x_le)
                continue;


            if (tets[e].xc > x25)
                continue;


            selected_global[e] = 1;

            hp_selected++;
        }


        printf(
            "HP%d: eta=%6.3f  "
            "y=%10.6f  "
            "LE=%10.6f  "
            "x25=%10.6f  "
            "c=%10.6f  "
            "seed=%d  "
            "selected=%d\n",
            hp + 1,
            eta[hp],
            y_hp,
            x_le,
            x25,
            chord,
            seed_count,
            hp_selected
        );


        free(seed);
        free(neighbor);
    }


    /* ========================================================
       Collect selected element IDs
       ======================================================== */

    int total_selected = 0;


    for (int e = 0; e < ntet; e++)
    {
        if (selected_global[e])
            total_selected++;
    }


    int *element_ids =
        (int *)malloc(
            (size_t)total_selected * sizeof(int));


    if (!element_ids)
    {
        fprintf(stderr,
                "ERROR allocating element_ids\n");

        exit(EXIT_FAILURE);
    }


    int count = 0;


    for (int e = 0; e < ntet; e++)
    {
        if (selected_global[e])
        {
            element_ids[count++] =
                tets[e].id;
        }
    }


    /*
       Sort element IDs so slatElementList.nam
       is deterministic and easy to inspect.
    */

    qsort(
        element_ids,
        total_selected,
        sizeof(int),
        compare_ints
    );


    /* ========================================================
       Write slatElementList.nam
       ======================================================== */

    out = fopen(output_file, "w");


    if (!out)
    {
        perror(
            "ERROR opening slatElementList.nam");

        exit(EXIT_FAILURE);
    }


    for (int i = 0;
         i < total_selected;
         i++)
    {
        fprintf(
            out,
            "%d\n",
            element_ids[i]
        );
    }


    fclose(out);


    printf("\n");
    printf("Total slat hard-point elements : %d\n",
           total_selected);

    printf("Written to                    : %s\n\n",
           output_file);


    /* ========================================================
       Cleanup
       ======================================================== */

    for (int i = 0; i < npoin; i++)
    {
        free(node_elem[i]);
    }

    free(node_elem);
    free(node_count);
    free(node_fill);

    free(selected_global);
    free(element_ids);

    free(nodes);
    free(tets);
}