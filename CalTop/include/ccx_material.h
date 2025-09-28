#pragma once
#include "CalculiX.h"
#include "stdio.h"

/* Fetch isotropic linear E, nu for an element's Gauss point.
   Returns 0 on success and sets *E, *nu; returns -1 if the material
   at this GP is not isotropic linear (mattyp!=1). */
int ccx_fetch_Enu_for_elem(
    ITG imat,                 /* material number (ielmat(1,elem)) */
    ITG iorien,               /* orientation number (0 if none) */
    const double pgauss[3],   /* GP coords (used for cylindrical systems) */
    double t0l, double t1l,   /* reference and current temperatures at GP */
    /* material DB (same arrays you already have) */
    const double *elcon, const ITG *nelcon,
    const double *rhcon,  const ITG *nrhcon,
    const double *alcon,  const ITG *nalcon, const double *alzero,
    ITG ntmat_, ITG ncmat_, ITG npmat_,
    const double *plicon, const double *plkcon,
    /* misc */
    ITG mi1, double dtime, ITG jj,
    /* outputs */
    double *E, double *nu
);
