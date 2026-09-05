#ifndef SU2_PREPROCESS_H
#define SU2_PREPROCESS_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>

#define MAXLINE 1024


/**
 * @brief Finds the SU2 solid mesh file in the current working directory.
 *
 * This function scans the current directory for a file ending with
 * the `.su2` extension. The discovered filename is copied into the
 * user-provided character buffer.
 *
 * @param filename Character array where the SU2 filename is stored.
 *
 * @note
 * The function exits with EXIT_FAILURE if no SU2 mesh is found.
 */
void find_su2_file(char *filename);

/**
 * @brief Finds the deformed SU2 solid mesh file in the current 
 * working directory.
 *
 * This function scans the current directory for a file ending with
 * the `.su2` extension. The discovered filename is copied into the
 * user-provided character buffer.
 *
 * @param filename Character array where the SU2 filename is stored.
 *
 * @note
 * The function exits with EXIT_FAILURE if no SU2 mesh is found.
 */
void find_AD_su2_file(char *filename);



/**
 * @brief Removes existing CalculiX/CalTop NAM files.
 *
 * Searches the current working directory and removes all files
 * with the `.nam` extension.
 *
 * This prevents stale mesh preprocessing files from being reused.
 */
void remove_existing_nam_files(void);



/**
 * @brief Prints all boundary markers present in an SU2 mesh.
 *
 * Searches for all MARKER_TAG entries in the SU2 file and prints
 * the marker names.
 *
 * @param su2file Input SU2 mesh filename.
 */
void print_su2_markers(char *su2file);

/**
 * @brief Prints all boundary markers present in an SU2 mesh and returns a list of markers.
 *
 * Searches for all MARKER_TAG entries in the SU2 file, prints
 * and returns the marker names.
 *
 * @param su2file Input SU2 mesh filename.
 */
char **get_su2_markers(char *su2file, int *marker_count);

/**
 * @brief Computes signed tetrahedral volume.
 *
 * Computes the determinant-based volume of a tetrahedral element.
 * The sign is used to detect incorrect element orientation.
 *
 * @param a Coordinates of node 1.
 * @param b Coordinates of node 2.
 * @param c Coordinates of node 3.
 * @param d Coordinates of node 4.
 *
 * @return Signed tetrahedron volume.
 */
double tet_volume(double a[3],
                  double b[3],
                  double c[3],
                  double d[3]);


/**
 * @brief Converts SU2 tetrahedral mesh to CalculiX format.
 *
 * Reads:
 *  - NPOIN node coordinates
 *  - NELEM tetrahedral connectivity
 *
 * Generates:
 *  - mesh.nam
 *
 * @param su2file Input SU2 mesh filename.
 */  
void convert_volume_mesh(char *su2file);


/**
 * @brief Extracts boundary nodes from an SU2 marker.
 *
 * Reads a given MARKER_TAG and extracts all associated
 * surface nodes.
 *
 * Generates:
 *  - Nfix1.nam for fixed markers
 *  - NSurface.nam for surface markers
 *
 * @param su2file Input SU2 mesh.
 * @param marker Marker name to extract.
 * @param outfile Output node set file.
 */
void extract_marker(char *su2file,
                    char *marker,
                    char *outfile);


/**
 * @brief Extracts passive skin elements.
 *
 * Finds all volume tetrahedral elements touching a specified
 * SU2 surface marker.
 *
 * Generates:
 *  - skinElementList.nam
 *
 * @param su2file Input SU2 mesh.
 * @param markername SU2 surface marker.
 * @param outfile Output element list filename.
 */
void extract_skin_elements(char *su2file,
                           char *markername,
                           char *outfile);


                           /**
 * @brief Writes a displaced SU2 mesh.
 *
 * Reads the SU2 mesh from CWD and updates nodal coordinates using
 * the structural displacement vector from CalTop.
 *
 * Updated coordinates:
 *
 *     x_new = x + ux
 *     y_new = y + uy
 *     z_new = z + uz
 *
 * where:
 *
 *     ux = v[4*i+1]
 *     uy = v[4*i+2]
 *     uz = v[4*i+3]
 *
 * The remainder of the SU2 file is preserved.
 *
 * @param nk Number of nodes.
 * @param v CalculiX displacement vector.
 */
void write_deformed_su2(int nk, double *v);


#endif