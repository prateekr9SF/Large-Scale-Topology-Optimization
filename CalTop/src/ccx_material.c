#include <string.h>
#include "CalculiX.h"
#include "ccx_material.h"

/* Fortran subroutine prototype (material property fetch) */
extern void materialdata_me_(
    const double *elcon, const ITG *nelcon,
    const double *rhcon, const ITG *nrhcon,
    const double *alcon, const ITG *nalcon,
    ITG *imat, char *amat,
    ITG *iorien, const double *pgauss, const double *orab,
    ITG *ntmat_, double *elas, double *rho,
    ITG *i, ITG *ithermal, const double *alzero, ITG *mattyp,
    double *t0l, double *t1l, ITG *ihyper, ITG *istiff,
    double *elconloc, double *eth, ITG *kode,
    const double *plicon, const ITG *nplicon,
    const double *plkcon, const ITG *nplkcon,
    ITG *npmat_, double *plconloc, ITG *mi1,
    double *dtime, ITG *jj, double *xstiff, ITG *ncmat_
);

/* Fortran interface is available through FORTRAN(...) macro from CalculiX. */
int ccx_fetch_Enu_for_elem(
    ITG imat, ITG iorien, const double pgauss[3],
    double t0l, double t1l,
    const double *elcon, const ITG *nelcon,
    const double *rhcon,  const ITG *nrhcon,
    const double *alcon,  const ITG *nalcon, const double *alzero,
    ITG ntmat_, ITG ncmat_, ITG npmat_,
    const double *plicon, const double *plkcon,
    ITG mi1, double dtime, ITG jj,
    double *E, double *nu)
{
    /* Scratch matching materialdata_me arguments */
    double elas[21];           /* Fortran 1..21 → C [0..20] */
    double eth[6];
    double elconloc[21];
    double plconloc[802];
    double orab_dummy[7] = {0.0};
    char amat_dummy[80]; memset(amat_dummy, ' ', sizeof(amat_dummy));

    ITG istiff = 0, ihyper = 0, mattyp = 0;
    ITG kode   = 0;                 /* small-strain, linear path */
    ITG ielem  = 1;                 /* dummy element ID for the call */
    ITG ithermal_local[2] = {0,0};  /* purely mechanical in this helper */
    double rho_dummy = 0.0;
    double xstiff_dummy[27];        /* not used, but required by interface */

    /* Hyperelastic flag in CalculiX: nelcon(1,imat) < 0  (Fortran) */
    /* nelcon is (2,ntmat_). We need its first entry for imat. */
    if (nelcon[0 + (imat-1)*2] < 0) ihyper = 1;

    FORTRAN(materialdata_me,(
        elcon, (ITG*)nelcon, rhcon, (ITG*)nrhcon, alcon, (ITG*)nalcon,
        &imat, amat_dummy, &iorien, (double*)pgauss, orab_dummy, &ntmat_,
        elas, &rho_dummy, &ielem, ithermal_local, (double*)alzero, &mattyp,
        &t0l, &t1l, &ihyper, &istiff, elconloc, eth,
        &kode, (double*)plicon, &npmat_, (double*)plkcon, &npmat_, plconloc,
        &mi1, &dtime, &jj, xstiff_dummy, &ncmat_
    ));

    /* Isotropic linear elasticity → mattyp==1, elas(1)=E, elas(2)=nu */
    if (mattyp != 1) return -1;

    *E  = elas[0];
    *nu = elas[1];
    return 0;
}
