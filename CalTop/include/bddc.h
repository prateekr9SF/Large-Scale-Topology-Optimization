/* CALTOP_SOLVER=bddc: PETSc BDDC solver (src/bddc.c; stubs without -DCALTOP_BDDC) */
#ifndef CALTOP_BDDC_H
#define CALTOP_BDDC_H
#ifndef ITG
#define ITG int
#endif
void bddc_launch(int *argc, char ***argv);
int  bddc_active(void);
void bddc_capture_begin(ITG ne, ITG *ipkon, char *lakon);
void bddc_elem(ITG e, ITG nope, ITG ndof, double *sm);
void bddc_capture_end(ITG *kon, ITG *ipkon, ITG ne, ITG *nactdof, ITG *mi, ITG nk, double *co, ITG neq);
int  bddc_solve(double *b, ITG neq, double tol, int adjoint, double *ad, double *au, ITG *jq, ITG *irow);
/* preconditioner rebuild: the caller assembles the matrix with the stresses
   set to zero inside bddc_capture_begin/end and calls bddc_rebuild(), (i) at
   the next Newton iteration when bddc_want_rebuild() says a solve
   needed more than CALTOP_BDDC_REBUILD CG iterations, (ii) once at the
   converged configuration before the adjoint solves (if bddc_active()).
   bddc_rebuild() releases the old preconditioner before the new set-up. */
int  bddc_want_rebuild(void);
void bddc_rebuild(void);

/* distributed element work: assembly (mafillsm) and internal forces
   (resultsmech) computed by the BDDC ranks, each for the elements of its own
   subdomain. The caller fills the arguments it has; bddc_dist_asm returns 1
   when the element stiffness was computed and captured on the ranks (the
   caller then only assembles the right-hand side), bddc_dist_res returns 1
   when fn, qa and nal (slot 0) and the element arrays hold the result. */
typedef struct {
  double *co; ITG *nk, *kon, *ipkon; char *lakon; ITG *ne, *nactdof, *mi;
  double *elcon; ITG *nelcon; double *rhcon; ITG *nrhcon; double *alcon; ITG *nalcon; double *alzero;
  ITG *ielmat, *ielorien, *norien; double *orab; ITG *ntmat_;
  double *plicon; ITG *nplicon; double *plkcon; ITG *nplkcon; ITG *npmat_; char *matname; ITG *ncmat_;
  double *thicke; ITG *ielprop; double *prop; double *design, *penal, *matdens, *physcon;
  ITG *ithermal, *iprestr, *nstate_, *mortar, *nasym, *ne0, *nbody, *nload;
  ITG *nmethod, *iperturb, *iexpl, *mass, *stiffness, *buckling, *intscheme, *istep, *iinc, *coriolis, *kscale, *network, *neq, *nzs;
  double *dtime, *ttime, *time, *reltime;
  double *vold, *sti;                                     /* assembly */
  double *v, *stx, *eme, *xstiff, *stiini, *ener, *eei, *enerini, *emeini, *rhomin;   /* internal forces */
  ITG *iout, *ielas, *icmd, calcul_fn, calcul_qa, calcul_cauchy, *nener, ikin;
  double *fn, *qa; ITG *nal;
} bddc_eargs;
int  bddc_dist_asm(bddc_eargs *a);
int  bddc_dist_on(void);   /* 1: element work is distributed (no per-thread copies needed) */
int  bddc_dist_res(bddc_eargs *a);
#endif
