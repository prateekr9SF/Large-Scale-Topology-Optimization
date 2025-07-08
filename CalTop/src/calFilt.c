
/* Insert header later*/


#ifdef __WIN32
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

#ifdef CALCULIX_MPI
#include <spoolesMPI.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"
#include <time.h>
#include <unistd.h>
#include <sys/stat.h> 



#ifdef CALCULIX_MPI
ITG myid = 0, nproc = 0;
#endif