/* Include the POSIX thread libraries for threading support.*/
#include <pthread.h>


/* Define constants for different operating systems to use in conditional compilation. */
#define Linux 1


// Define the FORTRAN macro differently based on the operating system.
#if ARCH == Linux
#define FORTRAN(A, B) A##_ B // For Linux, concatenate A and B with an underscore.
#endif
// Define the CEE macro similarly to the FORTRAN macro, for potentially different usage contexts.
#if ARCH == Linux
#define CEE(A, B) A##_ B // Same as FORTRAN macro on Linux.
#endif


// Memory management macros to simplify dynamic memory operations with debugging information.
#define NNEW(a, b, c) a = (b *)u_calloc((c), sizeof(b), __FILE__, __LINE__, #a) // Allocate memory and zero-initialize.
#define RENEW(a, b, c) a = (b *)u_realloc((b *)(a), (c) * sizeof(b), __FILE__, __LINE__, #a) // Reallocate memory.
#define SFREE(a) u_free(a, __FILE__, __LINE__, #a) // Free memory.

// Simplified version of memset for specific data manipulation within arrays.
#define DMEMSET(a, b, c, d)    \
	for (im = b; im < c; im++) \
	a[im] = d // Set elements from index b to c of array a to value d.

// Conditional compilation for supporting different integer types.
#ifdef LONGLONG
#define ITG long long // Use long long for integer type if LONGLONG is defined.
#define ITGFORMAT "lld" // Format specifier for long long.
#else
#define ITG int // Use int for integer type by default.
#define ITGFORMAT "d" // Format specifier for int.
#endif

/* Function declaration for Topology Optimization */

/* Function to compute element density*/
void rho(double *design, int ne);

/* Treat passive element */
void elementPassiveTreatment(ITG nset, ITG *ialset, char *set, ITG *istartset, 
                             ITG *iendset, double *vector, double setToValue); // Function to handle passive elements in a mesh.
char *getSubstring2(char *string, int position, int length); // Function to extract a substring.

void dVoldRhoPhys(double *eleVol, double *gradVol, ITG ne); // Function to compute volume gradient with respect to physical density.


#if ARCH == Linux
#define FORTRAN(A, B) A##_ B // For Linux, concatenate A and B with an underscore.
#endif

ITG strcmp1(const char *s1, const char *s2);

ITG strcmp2(const char *s1, const char *s2, ITG length);

ITG strcpy1(char *s1, const char *s2, ITG length);

void FORTRAN(stop, ());