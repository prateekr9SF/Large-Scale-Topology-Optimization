
#include <stdint.h>
#include <limits.h>
#include "CalculiX.h"

/**
 * @file tank_pressure.c
 * @brief Convert SU2 tank boundary triangles into CalculiX pressure loads.
 *
 * This module generates a CalculiX *DLOAD file for applying uniform
 * pressure to the internal surface of a tank represented by an SU2
 * boundary marker.
 *
 * The SU2 mesh contains two relevant types of elements:
 *
 *   1. Volume tetrahedra (SU2 type 10), which become C3D4 elements
 *      in CalculiX.
 *
 *   2. Surface triangles (SU2 type 5), which identify the boundary
 *      of the tank cavity through a marker such as "tank".
 *
 * The objective is to determine which tetrahedral element and local
 * face correspond to each tank surface triangle. CalculiX requires
 * this information to apply pressure using the P1-P4 face labels.
 *
 * ------------------------------------------------------------------
 * ALGORITHM
 * ------------------------------------------------------------------
 *
 * STEP 1: Read the SU2 nodal coordinates
 *
 *   The function searches for NPOIN and reads the coordinates of
 *   every mesh node into the xyz array.
 *
 *   These coordinates are needed to compute the signed volume of
 *   each tetrahedron and reproduce the orientation correction
 *   performed by convert_volume_mesh().
 *
 * STEP 2: Read the volume tetrahedra
 *
 *   The function rewinds the SU2 file, locates NELEM, and reads
 *   each tetrahedral element:
 *
 *       10 n1 n2 n3 n4 elementID
 *
 *   SU2 node indices are zero-based.
 *
 *   The original SU2 elementID is not used for CalculiX numbering.
 *   Instead, CalculiX elements are numbered sequentially according
 *   to their position in the volume-element section:
 *
 *       CalculiX element ID = e + 1
 *
 * STEP 3: Correct tetrahedral orientation
 *
 *   The signed volume of each tetrahedron is computed using
 *   tet_volume().
 *
 *   If the signed volume is negative, local nodes 2 and 3 are
 *   exchanged. This reproduces the connectivity correction in
 *   convert_volume_mesh().
 *
 *   This operation is essential because the CalculiX pressure face
 *   numbers depend on the final local ordering of tetrahedral nodes.
 *
 * STEP 4: Construct the tetrahedral face database
 *
 *   Every C3D4 tetrahedron has four triangular faces:
 *
 *       P1 : nodes 1, 2, 3
 *       P2 : nodes 1, 4, 2
 *       P3 : nodes 2, 4, 3
 *       P4 : nodes 3, 4, 1
 *
 *   For each tetrahedron, four TankFace records are generated.
 *
 *   Each record stores:
 *
 *       a, b, c : Three global node IDs defining the face
 *       element  : Corresponding CalculiX element ID
 *       face     : Corresponding CalculiX pressure face number
 *
 *   The three node IDs are sorted into ascending order using
 *   tank_sort3(). This produces an orientation-independent
 *   representation of each triangular face.
 *
 *   The complete face array therefore contains:
 *
 *       Number of face records = 4 * NELEM
 *
 * STEP 5: Sort the face database
 *
 *   The face records are sorted lexicographically according to
 *   their three node IDs using qsort() and tank_face_compare().
 *
 *   Sorting allows the function to locate matching tetrahedral
 *   faces efficiently using binary search rather than scanning
 *   every volume element for each tank triangle.
 *
 * STEP 6: Read the tank boundary marker
 *
 *   The SU2 file is rewound and searched for:
 *
 *       MARKER_TAG= tank
 *
 *   The associated MARKER_ELEMS section contains the triangular
 *   faces forming the tank cavity boundary:
 *
 *       5 n1 n2 n3
 *
 *   Each triangle is represented by its three zero-based SU2
 *   node IDs.
 *
 * STEP 7: Match each tank triangle to a tetrahedral face
 *
 *   The three node IDs of the tank triangle are sorted using
 *   tank_sort3().
 *
 *   A binary search (bsearch) is performed against the sorted
 *   TankFace array.
 *
 *   A successful match identifies:
 *
 *       - The C3D4 element adjacent to the tank surface.
 *       - The local face number (P1, P2, P3, or P4).
 *
 *   The function checks that the triangle matches exactly one
 *   tetrahedral face.
 *
 *   If multiple tetrahedra share the same face, the match is
 *   considered ambiguous because the function does not determine
 *   which side should receive the pressure load.
 *
 *   A separate used array prevents duplicate pressure entries
 *   for the same element face.
 *
 * STEP 8: Write the CalculiX pressure loads
 *
 *   For every successfully matched tank triangle, the function
 *   writes one distributed pressure load:
 *
 *       elementID, Pn, pressure
 *
 *   The output file begins with:
 *
 *       *DLOAD
 *
 *   For example, a uniform pressure of 4 bar produces:
 *
 *       *DLOAD
 *       100, P1, 400000
 *       105, P3, 400000
 *       110, P2, 400000
 *
 *   The generated file can be included in the CalculiX input
 *   deck inside the appropriate analysis step.
 *
 * STEP 9: Validate and release memory
 *
 *   The function terminates if it encounters:
 *
 *       - Invalid nodal or tetrahedral data.
 *       - Zero-volume tetrahedra.
 *       - Missing or empty tank markers.
 *       - Tank triangles without matching tetrahedral faces.
 *       - Multiple tetrahedral faces matching one tank triangle.
 *       - Duplicate tank triangles.
 *       - Memory allocation failures.
 *
 *   After successful extraction, the function closes the files,
 *   releases the allocated arrays, and returns the number of
 *   pressure-loaded faces.
 *
 * ------------------------------------------------------------------
 * COMPUTATIONAL CHARACTERISTICS
 * ------------------------------------------------------------------
 *
 * The algorithm stores four TankFace records per tetrahedron.
 *
 * Sorting the face database requires approximately:
 *
 *       O(NELEM log NELEM)
 *
 * Matching NTANK boundary triangles using binary search requires:
 *
 *       O(NTANK log NELEM)
 *
 * The face database requires O(NELEM) memory.
 *
 * ------------------------------------------------------------------
 * ASSUMPTIONS
 * ------------------------------------------------------------------
 *
 *   - The volume mesh consists of linear C3D4 tetrahedra.
 *   - The tank marker contains triangular SU2 boundary elements.
 *   - Each tank triangle corresponds to exactly one solid face.
 *   - The pressure is uniform over the selected tank surface.
 *   - Element numbering and orientation match convert_volume_mesh().
 *   - The pressure sign follows the CalculiX P1-P4 convention.
 *
 * @note
 * This module generates the pressure load definitions only.
 * The actual pressure forces are assembled by CalculiX during
 * the finite-element analysis.
 *
 * @note
 * The tank surface is obtained directly from the SU2 boundary
 * marker. The passive-element list in tankElementList.nam is
 * not used for pressure-face identification.
 */



/*double tet_volume(double a[3],
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
}*/


typedef struct 
{
    int a, b, c;        /* sorted, zero-based node IDs */
    int element, face;  /* CalculiX element number and P face */
} TankFace;

static void tank_sort3(int *a, int *b, int *c)
{
    int t;
    if (*a > *b) { t=*a; *a=*b; *b=t; }
    if (*b > *c) { t=*b; *b=*c; *c=t; }
    if (*a > *b) { t=*a; *a=*b; *b=t; }
}

static int tank_face_compare(const void *pa, const void *pb)
{
    const TankFace *a = pa, *b = pb;
    if (a->a != b->a) return (a->a > b->a) - (a->a < b->a);
    if (a->b != b->b) return (a->b > b->b) - (a->b < b->b);
    return (a->c > b->c) - (a->c < b->c);
}

/* Returns number of unique loaded faces; exits on invalid or ambiguous matches.
 * Requires the same mesh connectivity and numbering conventions as convert_volume_mesh().
 * Output is a CalculiX *DLOAD block for inclusion INSIDE a *STEP.
 */
size_t extract_pressure_faces(const char *su2file, const char *marker,
                              const char *outfile, double pressure)
{
    FILE *fp = fopen(su2file, "r");
    char line[MAXLINE], tag[256];
    int nnode=0, nelem=0;
    double (*xyz)[3] = NULL;
    TankFace *faces = NULL;
    unsigned char *used = NULL;
    size_t nf, loaded=0;
    FILE *out = NULL;

    /* CalculiX C3D4 pressure faces, zero-based local node indices */
    const int local_face[4][3] = 
    {
        {0,1,2}, {0,3,1}, {1,3,2}, {2,3,0}
    };

    if (!fp || pressure <= 0.0) 
    {
        fprintf(stderr, "ERROR: Cannot open SU2 mesh or invalid pressure.\n");
        exit(EXIT_FAILURE);
    }

    /* First pass: read coordinates; SU2 commonly lists NELEM before NPOIN. */
    while (fgets(line, sizeof line, fp)) 
    {
        if (sscanf(line, "NPOIN= %d", &nnode) == 1) 
        {
            if (nnode <= 0 || xyz) goto bad_mesh;
            xyz = malloc((size_t)nnode * sizeof *xyz);
            if (!xyz) goto out_of_memory;
            for (int i=0; i<nnode; ++i) {
                if (!fgets(line,sizeof line,fp) ||
                    sscanf(line,"%lf %lf %lf",&xyz[i][0],&xyz[i][1],&xyz[i][2]) != 3)
                    goto bad_mesh;
            }
        }
    }

    rewind(fp);

    /* Second pass: read volume connectivity, now that coordinates are available. */
    while (fgets(line, sizeof line, fp)) 
    {
        if (sscanf(line, "NELEM= %d", &nelem) == 1) {
            if (nelem <= 0 || !xyz || (size_t)nelem > SIZE_MAX/(4*sizeof *faces))
                goto bad_mesh;
            nf = (size_t)nelem*4;
            faces = malloc(nf * sizeof *faces);
            if (!faces) goto out_of_memory;
            for (int e=0; e<nelem; ++e) {
                int type, n[4], eid;
                if (!fgets(line,sizeof line,fp) ||
                    sscanf(line,"%d %d %d %d %d %d",&type,&n[0],&n[1],&n[2],&n[3],&eid)!=6 || type!=10)
                    goto bad_mesh;
                for (int k=0;k<4;k++) if (n[k]<0 || n[k]>=nnode) goto bad_mesh;
                double vol=tet_volume(xyz[n[0]],xyz[n[1]],xyz[n[2]],xyz[n[3]]);
                if (fabs(vol)<1.0e-30) goto bad_mesh;
                if (vol<0) { int t=n[1];n[1]=n[2];n[2]=t; }
                for (int f=0; f<4; ++f) {
                    TankFace *p=&faces[4*(size_t)e+f];
                    p->a=n[local_face[f][0]];
                    p->b=n[local_face[f][1]];
                    p->c=n[local_face[f][2]];
                    tank_sort3(&p->a,&p->b,&p->c);
                    p->element=e+1;
                    p->face=f+1;
                }
            }
        }
    }
    if (!faces || !xyz) goto bad_mesh;
    free(xyz); xyz=NULL;
    nf=(size_t)nelem*4;
    qsort(faces,nf,sizeof *faces,tank_face_compare);
    used=calloc(nf,sizeof *used);
    if (!used) goto out_of_memory;
    rewind(fp);
    out=fopen(outfile,"w");
    if (!out) { perror(outfile); exit(EXIT_FAILURE); }
    fprintf(out,"*DLOAD\n");
    int in_marker=0, found=0, remaining=0;
    while (fgets(line,sizeof line,fp)) {
        if (sscanf(line,"MARKER_TAG= %255s",tag)==1) {
            in_marker=(strcmp(tag,marker)==0);
            if (in_marker) { found=1; remaining=-1; }
            continue;
        }
        if (!in_marker) continue;
        if (sscanf(line,"MARKER_ELEMS= %d",&remaining)==1) continue;
        if (remaining<=0) continue;
        int type,a,b,c;
        if (sscanf(line,"%d %d %d %d",&type,&a,&b,&c)!=4 || type!=5 ||
            a<0 || b<0 || c<0 || a>=nnode || b>=nnode || c>=nnode)
            goto bad_mesh;
        --remaining;
        tank_sort3(&a,&b,&c);
        TankFace key={.a=a,.b=b,.c=c,.element=0,.face=0};
        TankFace *match=bsearch(&key,faces,nf,sizeof *faces,tank_face_compare);
        if (!match) {
            fprintf(stderr,"ERROR: tank triangle (%d,%d,%d) has no C3D4 face.\n",a,b,c);
            exit(EXIT_FAILURE);
        }
        /* A face shared by multiple tetrahedra is ambiguous without a side designation. */
        size_t idx=(size_t)(match-faces);
        size_t first=idx, last=idx;
        while (first>0 && tank_face_compare(&faces[first-1],&key)==0) --first;
        while (last+1<nf && tank_face_compare(&faces[last+1],&key)==0) ++last;
        if (first!=last) {
            fprintf(stderr,"ERROR: tank triangle (%d,%d,%d) belongs to %zu tetrahedral faces; cavity side ambiguous.\n",a,b,c,last-first+1);
            exit(EXIT_FAILURE);
        }
        if (used[idx]) {
            fprintf(stderr,"ERROR: duplicate tank triangle (%d,%d,%d).\n",a,b,c);
            exit(EXIT_FAILURE);
        }
        used[idx]=1;
        fprintf(out,"%d, P%d, %.12g\n",match->element,match->face,pressure);
        ++loaded;
    }
    if (!found || loaded==0 || remaining>0) {
        fprintf(stderr,"ERROR: marker '%s' missing, empty, or truncated.\n",marker);
        exit(EXIT_FAILURE);
    }
    fclose(out); fclose(fp); free(faces); free(used);
    printf("Tank pressure: wrote %zu faces at %.9g Pa to %s\n",loaded,pressure,outfile);
    return loaded;

bad_mesh:
    fprintf(stderr,"ERROR: malformed/unsupported SU2 mesh while extracting tank pressure.\n");
    exit(EXIT_FAILURE);
out_of_memory:
    fprintf(stderr,"ERROR: out of memory extracting tank pressure.\n");
    exit(EXIT_FAILURE);
}
