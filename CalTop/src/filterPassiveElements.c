#include <stdio.h>

/**
 * Sets the design variable to 1.0 for passive elements.
 *
 * @param design        Pointer to the filtered/un-filtered density array (length = ne).
 * @param ne            Number of elements in the mesh.
 * @param passiveIDs    Array of passive element IDs (1-based).
 * @param numPassive    Number of passive element IDs.
 */
void filterOutPassiveElems_density(double *design, int ne, int *passiveIDs, int numPassive) 
{
    for (int i = 0; i < numPassive; i++) 
    {
        int eid = passiveIDs[i];  // 1-based ID

        if (eid >= 1 && eid <= ne) 
        {
            design[eid - 1] = 0.0;  // Convert to 0-based index
        } 
        else 
        {
            printf("Warning: Passive element ID %d is out of bounds [1, %d]\n", eid, ne);
        }
    }
}

/**
 * Sets the compliance sensitivity to 0.0 for passive elements.
 *
 * @param sens          Pointer to the filtered sens array (length = ne).
 * @param ne            Number of elements in the mesh.
 * @param passiveIDs    Array of passive element IDs (1-based).
 * @param numPassive    Number of passive element IDs.
 */
void filterOutPassiveElems_sens(double *sens, int ne, int *passiveIDs, int numPassive) 
{
    for (int i = 0; i < numPassive; i++) 
    {
        int eid = passiveIDs[i];  // 1-based ID

        if (eid >= 1 && eid <= ne) 
        {
            sens[eid - 1] = 0.0;  // Convert to 0-based index
        } 
        else 
        {
            printf("Warning: Passive element ID %d is out of bounds [1, %d]\n", eid, ne);
        }
    }
}
