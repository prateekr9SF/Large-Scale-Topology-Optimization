/* resultsmech_c3d4_linear.c
   Minimal C path for linear-elastic C3D4 stresses (one GP).
   - Uses materialdata_me via ccx_fetch_Enu_for_elem() to get E,nu
   - Writes stx(6,1,i) in CalculiX storage (Fortran layout)
*/

#include <math.h>
#include <string.h>
#include "CalculiX.h"
#include "ccx_material.h"    /* ccx_fetch_Enu_for_elem(E,nu) wrapper */

/* ---- small utilities -------------------------------------------------- */

static inline const char* lakon_i(const char *lakon, ITG i)
/* Return pointer to the 8-char lakon record for element i (1-based). */
{
    return lakon + 8*(i-1);
}

static inline int is_C3D4(const char *lk)
/* True if lakon starts with "C3D4" (ignores trailing chars/spaces). */
{
    return (lk[0]=='C' && lk[1]=='3' && lk[2]=='D' && lk[3]=='4');
}

/* Indexing helpers for Fortran-shaped arrays seen from C
   Fortran: co(3, nk)           → C: co[ (k-1) + (node-1)*3 ], k=1..3
   Fortran: v(0:mi2, nk) (mt=mi2+1) → C: v[ k + (node-1)*mt ], k=1..3
   Fortran: stx(6, mi1, ne)     → C: stx[ (k-1) + (jj-1)*6 + (elem-1)*6*mi1 ]
*/

static inline double co_at(const double *co, ITG node, int k)
/* k=1..3 */
{
    return co[(k-1) + (node-1)*3];
}

static inline double v_at(const double *v, ITG node, int k, ITG mt)
/* k=1..3; mt = mi[1]+1 */
{
    return v[k + (node-1)*mt];
}

static inline double* stx_ptr(double *stx, ITG elem, ITG jj, ITG mi1)
/* Returns base pointer to stx(:,jj,elem) (Fortran layout) */
{
    return stx + (jj-1)*6 + (elem-1)*6*mi1;
}

/* ---- Tet4 geometry: constant gradients and volume --------------------- */

typedef struct {
    double dNdx[4][3];  /* ∂N_a/∂x_k */
    double V;           /* element volume */
} Tet4Geom;

static inline double det3(const double A[3][3])
{
    return  A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
          - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
          + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}

static int tet4_geom(const double x[4][3], Tet4Geom *tg)
/* x[a][k], a=0..3, k=0..2 */
{
    double J[3][3] = {
        { x[1][0]-x[0][0], x[2][0]-x[0][0], x[3][0]-x[0][0] },
        { x[1][1]-x[0][1], x[2][1]-x[0][1], x[3][1]-x[0][1] },
        { x[1][2]-x[0][2], x[2][2]-x[0][2], x[3][2]-x[0][2] }
    };
    double detJ = det3(J);
    double V = detJ / 6.0;
    if (V <= 0.0) return -1; /* inverted/degenerate */

    double invJ[3][3];
    /* inverse(J) = adj(J)/detJ ; we need invJ^T for gradients */
    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1]) / detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1]) / detJ;
    invJ[0][2] =  (J[0][1]*J[1][2]-J[0][2]*J[1][1]) / detJ;

    invJ[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0]) / detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0]) / detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0]) / detJ;

    invJ[2][0] =  (J[1][0]*J[2][1]-J[1][1]*J[2][0]) / detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0]) / detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0]) / detJ;

    /* Reference ∇_ξ N: N0=[-1,-1,-1], N1=[1,0,0], N2=[0,1,0], N3=[0,0,1] */
    const double dNr[4][3] = {
        {-1.0, -1.0, -1.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0}
    };
    for (int a=0; a<4; ++a)
        for (int k=0; k<3; ++k)
            tg->dNdx[a][k] = invJ[0][k]*dNr[a][0] + invJ[1][k]*dNr[a][1] + invJ[2][k]*dNr[a][2];

    tg->V = V;
    return 0;
}

/* ---- Build B, iso D, and compute stress ------------------------------- */

static void build_B(const Tet4Geom *tg, double B[6][12])
/* Voigt order: (xx,yy,zz,xy,xz,yz); tensor shear (ε_xy = 1/2(du_x/dy+du_y/dx)) */
{
    for (int i=0;i<6;++i) for (int j=0;j<12;++j) B[i][j]=0.0;

    for (int a=0;a<4;++a) {
        const double dNx = tg->dNdx[a][0];
        const double dNy = tg->dNdx[a][1];
        const double dNz = tg->dNdx[a][2];
        const int c = 3*a;

        /* normal */
        B[0][c+0] = dNx;
        B[1][c+1] = dNy;
        B[2][c+2] = dNz;
        /* tensor shear */
        B[3][c+0] = 0.5*dNy; B[3][c+1] = 0.5*dNx; /* ε_xy */
        B[4][c+0] = 0.5*dNz; B[4][c+2] = 0.5*dNx; /* ε_xz */
        B[5][c+1] = 0.5*dNz; B[5][c+2] = 0.5*dNy; /* ε_yz */
    }
}

static void iso_D(double E, double nu, double D[6][6])
/* Isotropic 3D elasticity (Voigt: xx,yy,zz,xy,xz,yz), tensor shear */
{
    const double lam = (nu*E)/((1.0+nu)*(1.0-2.0*nu));
    const double G   = E/(2.0*(1.0+nu));
    for (int i=0;i<6;++i) for (int j=0;j<6;++j) D[i][j]=0.0;
    D[0][0]=lam+2*G; D[0][1]=lam;     D[0][2]=lam;
    D[1][0]=lam;     D[1][1]=lam+2*G; D[1][2]=lam;
    D[2][0]=lam;     D[2][1]=lam;     D[2][2]=lam;
    D[3][3]=G; D[4][4]=G; D[5][5]=G;
}

/* y = A(6x12) * x(12) */
static inline void matvec_6x12(const double A[6][12], const double x[12], double y[6])
{
    for (int i=0;i<6;++i) {
        double s=0.0;
        for (int j=0;j<12;++j) s += A[i][j]*x[j];
        y[i]=s;
    }
}

/* y = A(6x6) * x(6) */
static inline void matvec_6x6(const double A[6][6], const double x[6], double y[6])
{
    for (int i=0;i<6;++i) {
        double s=0.0;
        for (int j=0;j<6;++j) s += A[i][j]*x[j];
        y[i]=s;
    }
}

/* ---- Public entry ----------------------------------------------------- */

void resultsmech_c3d4_linear(
    const double *co,       /* co(3,*) */
    const ITG    *kon,      /* global connectivity (1-based node IDs) */
    const ITG    *ipkon,    /* ipkon(i): 1-based start index in kon for elem i */
    const char   *lakon,    /* 8-char type per element */
    const ITG    *ielmat,   /* ielmat(layer,elem) in C: ielmat[(elem-1)*mi3 + (layer-1)] */
    const ITG    *ielorien, /* may be NULL if norien==0 */
    ITG norien,
    const double *v,        /* v(0:mi(2),*) nodal displacements */
    double *stx,            /* stx(6,mi(1),*) output stresses */
    /* material DB */
    const double *elcon, const ITG *nelcon,
    const double *rhcon,  const ITG *nrhcon,
    const double *alcon,  const ITG *nalcon, const double *alzero,
    ITG ntmat_, ITG ncmat_, ITG npmat_,
    const double *plicon, const double *plkcon,
    /* sizes */
    const ITG *mi,         /* mi(1)=mi1 (#GP lines), mi(2)=#dof-1 per node, mi(3)=max layers */
    /* element range to process (1-based, inclusive) */
    ITG nea, ITG neb
)
{
    const ITG mi1 = mi[0];          /* #integration point lines in stx dimension */
    const ITG mi2 = mi[1];          /* dof-1 per node → mt = mi2+1 */
    const ITG mi3 = mi[2];          /* layer slots */
    const ITG mt  = mi2 + 1;

    /* Only jj=1 exists for C3D4 (single GP) */
    const ITG jj  = 1;

    for (ITG i = nea; i <= neb; ++i)
    {
        if (ipkon[i-1] < 0) continue;                  /* inactive */
        const char *lk = lakon_i(lakon, i);
        if (lk[0]=='F') continue;                      /* fluid */
        if (!is_C3D4(lk)) continue;                    /* only C3D4 */

        /* Get element connectivity (4 nodes) */
        ITG start = ipkon[i-1];    /* Fortran 1-based pointer into kon */
        /* Convert to C 0-based array index */
        ITG k0 = start - 1;
        ITG n1 = kon[k0+0];
        ITG n2 = kon[k0+1];
        ITG n3 = kon[k0+2];
        ITG n4 = kon[k0+3];

        /* Coordinates */
        double x[4][3] = {
            { co_at(co, n1,1), co_at(co, n1,2), co_at(co, n1,3) },
            { co_at(co, n2,1), co_at(co, n2,2), co_at(co, n2,3) },
            { co_at(co, n3,1), co_at(co, n3,2), co_at(co, n3,3) },
            { co_at(co, n4,1), co_at(co, n4,2), co_at(co, n4,3) }
        };

        Tet4Geom tg;
        if (tet4_geom(x, &tg) != 0) {
            /* degenerate/inverted → write zeros */
            double *st = stx_ptr(stx, i, jj, mi1);
            for (int k=0;k<6;++k) st[k] = 0.0;
            continue;
        }

        /* Displacements (ux,uy,uz) for each node */
        double ue[12] = {
            v_at(v,n1,1,mt), v_at(v,n1,2,mt), v_at(v,n1,3,mt),
            v_at(v,n2,1,mt), v_at(v,n2,2,mt), v_at(v,n2,3,mt),
            v_at(v,n3,1,mt), v_at(v,n3,2,mt), v_at(v,n3,3,mt),
            v_at(v,n4,1,mt), v_at(v,n4,2,mt), v_at(v,n4,3,mt)
        };

        /* Gauss point position (centroid) — used by materialdata for some systems */
        double xgp[3] = {
            0.25*(x[0][0]+x[1][0]+x[2][0]+x[3][0]),
            0.25*(x[0][1]+x[1][1]+x[2][1]+x[3][1]),
            0.25*(x[0][2]+x[1][2]+x[2][2]+x[3][2])
        };

        /* Material selection for this element (layer 1) */
        ITG imat   = ielmat ? ielmat[(i-1)*mi3 + 0] : 1;
        ITG iorien = (norien>0 && ielorien) ? ielorien[(i-1)*mi3 + 0] : 0;

        /* Fetch E and nu through CalculiX material path */
        double E=0.0, nu=0.0;
        if (ccx_fetch_Enu_for_elem(imat, iorien, xgp,
                /* t0l,t1l: set to 0 here; For thermo-mech, compute at GP */
                0.0, 0.0,
                elcon, nelcon, rhcon, nrhcon, alcon, nalcon, alzero,
                ntmat_, ncmat_, npmat_, plicon, plkcon,
                mi[0], /* dtime */ 0.0, /* jj */ 1,
                &E, &nu) != 0)
        {
            /* Not isotropic linear → write zeros and continue */
            double *st = stx_ptr(stx, i, jj, mi1);
            for (int k=0;k<6;++k) st[k] = 0.0;
            continue;
        }

        /* Operators */
        double B[6][12]; build_B(&tg, B);
        double D[6][6];  iso_D(E, nu, D);

        /* Strain and Stress (tensor shear conventions) */
        double eps[6];   matvec_6x12(B, ue, eps);
        double sig[6];   matvec_6x6(D, eps, sig);

        /* Store in stx(:,1,i) */
        double *st = stx_ptr(stx, i, jj, mi1);
        /* Voigt order: 1..6 = (xx,yy,zz,xy,xz,yz) */
        st[0] = sig[0];
        st[1] = sig[1];
        st[2] = sig[2];
        st[3] = sig[3];
        st[4] = sig[4];
        st[5] = sig[5];
    }
}
