/*
 * ----------------------------------------------------------------------
 * CALTOP_SOLVER=bddc: PETSc BDDC solver for the NLGEOM runs (nonlingeo)
 * ----------------------------------------------------------------------
 * Process model: calTop.exe relaunches itself under MPI ($CALTOP_MPIEXEC,
 * default "mpiexec -n $CALTOP_BDDC_NP --bind-to core"; the FADO command does
 * not change). Rank 0 runs CalTop (OpenMP) and owns subdomain 0; ranks
 * 1..P-1 wait in worker_loop for commands. One METIS subdomain per rank.
 *
 * Per run
 *  - first assembly: rank 0 captures every element matrix (bddc_elem, via
 *    ctelemcapture_ from mafillsm.f), partitions the mesh (METIS, contiguous)
 *    and sends each rank its elements (bddc_capture_end, CMD_TOPO).
 *  - from the second assembly on: element stiffness and internal forces are
 *    computed on the ranks, each for its own elements (bddc_dist_asm,
 *    bddc_dist_res; "distributed element work" below). Element matrices
 *    stay on their rank; CalTop's global K and its structure are never built.
 *  - preconditioner: built from the first (stress-free) matrix, PSD in
 *    every subdomain; CG on the current K_T for all Newton solves and the
 *    adjoints. Rebuilt from the stress-free matrix at the current
 *    configuration (i) at the next Newton iteration after a solve
 *    needs more than CALTOP_BDDC_REBUILD iterations and (ii) once at the
 *    converged configuration before the adjoint solves (nonlingeo.c). A
 *    rebuild frees the old preconditioner first (new PC, pc_customize).
 *  - PETSc: KSPCG + PCBDDC (faces, corner selection, deluxe scaling, adaptive
 *    coarse space, threshold 5), MUMPS local solves; CALTOP_BDDC_OPTS adds
 *    options. CG stops after 500 iterations; a solve that stops there is
 *    accepted if its true residual is below max(tol, 1e-5).
 *  - elements: C3D4, C3D10 or a mix (12x12 or 30x30 matrix per element).
 *  - no fallback: anything BDDC cannot do stops the run with *ERROR.
 *
 * Environment: CALTOP_SOLVER=bddc selects this solver (nr and nk are the
 * PARDISO solvers, pardiso_select_solver), CALTOP_BDDC_NP (24), CALTOP_MPIEXEC,
 * CALTOP_BDDC_TOL (Newton, 1e-3), CALTOP_BDDC_TOL_ADJ (adjoints, 1e-5),
 * CALTOP_BDDC_REBUILD (150, 0 = never during Newton), CALTOP_BDDC_OPTS,
 * CALTOP_BDDC_PROFILE=1 (coarse size and per-rank setup table).
 *
 * Build: make BDDC=1 links with -Wl,--exclude-libs,ALL. CalculiX bundles
 * reference BLAS/LAPACK (dgesv.f, dgmres.f); without the flag MUMPS and PETSc
 * call those instead of MKL and the setup is about 3x slower.
 *
 * Sections: launch | capture | collective pieces (topology, matrix, solve)
 * | distributed element work | worker loop | rank 0 (partition, ship,
 * rebuild, solve) | stubs (build without PETSc)
 * ----------------------------------------------------------------------
 */
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <execinfo.h>
#include <malloc.h>
#include <omp.h>
#include "CalculiX.h"
#include "bddc.h"
/* TAU regions (profiling build): see the category list in nonlingeo.c.
   Every rank writes its own profile (node id = rank); rank 0's profile is
   the run's breakdown. */
#ifdef PROFILING_ON
#include <TAU.h>
#define TAU_REGION_START(n) TAU_START(n)
#define TAU_REGION_STOP(n) do{TAU_TRACK_MEMORY_FOOTPRINT_HERE();TAU_STOP(n);}while(0)
#else
#define TAU_REGION_START(n) ((void)0)
#define TAU_REGION_STOP(n) ((void)0)
#endif

#ifdef CALTOP_BDDC
#include <petscksp.h>
#include <petsc/private/pcbddcimpl.h>
#if defined(PETSC_HAVE_METIS)
#include <metis.h>
#endif
#include "mkl_service.h"

enum { CMD_TOPO = 1, CMD_MATRIX = 2, CMD_SOLVE = 3, CMD_PMAT = 4, CMD_DMODEL = 5, CMD_DASM = 6, CMD_DRES = 7, CMD_LMATRIX = 8,
       CMD_QUIT = 9, CMD_LPMAT = 10 };
#define TAG_CMD 7101
#define TAG_DATA 7102

static int  bd_on = 0;            /* this process runs BDDC (rank 0) */
static cpu_set_t bd_core, bd_all; /* rank 0: own core (solves), all cores (CalTop) */
static void pin_own(void)  { sched_getaffinity(0, sizeof(bd_all), &bd_all); sched_setaffinity(0, sizeof(bd_core), &bd_core); }
static void pin_back(void) { sched_setaffinity(0, sizeof(bd_all), &bd_all); }
/* after a solve: CalTop's thread counts back, every OpenMP thread free again */
static void unpin_all(int nt, int nm)
{
  pin_back();
  omp_set_num_threads(nt);
  mkl_set_num_threads(nm);
#pragma omp parallel num_threads(nt)
  sched_setaffinity(0, sizeof(bd_all), &bd_all);
}
static int  bd_rank = 0, bd_size = 1, bd_disabled = 0;

/* ---- rank 0: captured element matrices ---- */
static double *bd_K = NULL;       /* bd_s2 per element (row-major bd_s x bd_s) */
static int     bd_nn = 10, bd_s = 30, bd_s2 = 900;   /* nodes, unknowns, matrix entries per element slot */
static unsigned char *bd_enope = NULL;               /* nodes of each element (4 or 10) */
static char   *bd_have = NULL;
static ITG     bd_ne = 0;
static int     bd_wcap0 = 0;     /* this rank computes its own element matrices (distributed assembly) */
static int     bd_dist;          /* rank 0: 1 = element work distributed (bddc_dist_asm) */
static int     bd_capturing = 0, bd_dirty = 0, bd_topo = 0;
static int     bd_nel = 0, *bd_el = NULL;  /* captured elements (index) */
static int    *bd_nodeid = NULL;            /* bs=3: PETSc node number of each CalTop node (-1: fixed) */
static int    *bd_edof = NULL;    /* per captured element, bd_s equations (-1: fixed or padding) */
static int    *bd_part = NULL;    /* subdomain of each captured element */
static MPI_Datatype *bd_type = NULL; /* per rank: its elements in bd_K */
static int    *bd_cnt = NULL, *bd_dsp = NULL; /* vector scatter layout */
static int    *bd_own0 = NULL, bd_n0 = 0;     /* rank 0's elements */
static int    *bd_g = NULL;                   /* CalTop equation -> PETSc index */
static double *bd_vp = NULL;                  /* vector in PETSc order (rank 0) */
static double  bd_tsetup = 0., bd_tsolve = 0., bd_tship = 0.;
static long    bd_nsolve = 0, bd_nit = 0, bd_nfail = 0, bd_nmat = 0;
static int     bd_rebuild_its = -1;            /* CALTOP_BDDC_REBUILD (0: never) */
static int     bd_want = 0, bd_nrebuild = 0;   /* rank 0: rebuild requested / done */
static double  bd_trebuild = 0.;
#ifndef BDDC_REBUILD_DEFAULT
#define BDDC_REBUILD_DEFAULT 150
#endif

/* ---- every rank: its subdomain ---- */
static PetscInt n_g = 0, m_own = 0, r_0 = 0;  /* global size, owned rows */
static PetscInt bs_g = 1;                     /* 3: unknowns grouped per node */
static int      nown = 0, *ledof = NULL;      /* own elements, local dofs */
static PetscInt nloc = 0, *lglob = NULL;      /* local dof -> equation */
static double  *lK = NULL;                    /* own element matrices */
static Mat      A = NULL, B = NULL;
static KSP      ksp = NULL;
static Vec      vx = NULL, vb = NULL, vr = NULL;
static int      pc_ready = 0;
static Vec      bd_coords = NULL;             /* owned node coordinates (bs 3) */
/* BDDC customization, kept so a rebuild can start from a new PC (pc_customize) */
static PetscInt  *bd_xa = NULL, *bd_aj = NULL; /* local adjacency graph (CSR) */
static IS         bd_fld[3] = {NULL, NULL, NULL}; /* x / y / z unknowns */
static PetscReal *bd_cr = NULL;                /* coordinates of the owned equations */

static double getd(const char *name, double def)
{
  char *e = getenv(name);
  return (e != NULL && e[0]) ? atof(e) : def;
}

typedef struct { double v; int i; } slab_t;
static int cmpslab(const void *x, const void *y)
{
  double d = ((const slab_t *)x)->v - ((const slab_t *)y)->v;
  return (d > 0) - (d < 0);
}

static int cmpint(const void *a, const void *b)
{
  PetscInt x = *(const PetscInt *)a, y = *(const PetscInt *)b;
  return (x > y) - (x < y);
}

static int cmp_ll2(const void *a, const void *b)
{
  long long x = *(const long long *)a, y = *(const long long *)b;
  return (x > y) - (x < y);
}

/* ---------------- launch ---------------- */

static int under_mpi(void)
{
  return getenv("CALTOP_BDDC_CHILD") != NULL || getenv("OMPI_COMM_WORLD_SIZE") != NULL || getenv("PMI_SIZE") != NULL ||
         getenv("PMIX_RANK") != NULL;
}

static void bddc_finish(void);
static void worker_loop(void);

static void crash_handler(int sig)
{
  void       *bt[64];
  int         n = backtrace(bt, 64);
  const char *m = "\n*** calTop rank 0 crashed (signal ";
  char        num[16];
  int         k = snprintf(num, sizeof(num), "%d", sig);
  fflush(stdout);
  if (write(1, m, strlen(m)) < 0 || write(1, num, k) < 0 || write(1, "), backtrace:\n", 14) < 0) {}
  backtrace_symbols_fd(bt, n, 1);
  if (write(1, "*** end of backtrace\n", 21) < 0) {}
  signal(sig, SIG_DFL);
  raise(sig);
}

void bddc_launch(int *argc, char ***argv)
{
  char *s = getenv("CALTOP_SOLVER");
  if (s == NULL || strcmp(s, "bddc") != 0) return;
  if (!under_mpi()) {
    /* relaunch this executable under MPI with the same arguments */
    char  launcher[2048], exe[4096], *tok[256], *nargv[512];
    int   nt = 0, i, np = (int)getd("CALTOP_BDDC_NP", 24.);
    char *l = getenv("CALTOP_MPIEXEC"), *p;
    ssize_t k = readlink("/proc/self/exe", exe, sizeof(exe) - 1);
    if (k <= 0) { strncpy(exe, (*argv)[0], sizeof(exe) - 1); k = strlen(exe); }
    exe[k] = '\0';
    if (l != NULL && l[0]) snprintf(launcher, sizeof(launcher), "%s", l);
    else snprintf(launcher, sizeof(launcher), "mpiexec -n %d --bind-to core", np);
    for (p = strtok(launcher, " \t"); p != NULL && nt < 250; p = strtok(NULL, " \t")) tok[nt++] = p;
    for (i = 0; i < nt; i++) nargv[i] = tok[i];
    nargv[nt] = exe;
    for (i = 1; i < *argc && nt + i < 510; i++) nargv[nt + i] = (*argv)[i];
    nargv[nt + i] = NULL;
    setenv("CALTOP_BDDC_CHILD", "1", 1);
    /* the solver ranks start with OMP_NUM_THREADS / MKL_NUM_THREADS unset:
       OpenMP then uses the CPUs of the rank's
       binding (2 hyperthreads of its core). With both set to 1 the same
       setup was 3-4x slower (wing: 407 s vs 136 s). Rank 0 restores the
       user's values for CalTop. */
    {
      char *o = getenv("OMP_NUM_THREADS"), *m = getenv("MKL_NUM_THREADS");
      setenv("CALTOP_BDDC_OMP0", (o && o[0]) ? o : "", 1);
      setenv("CALTOP_BDDC_MKL0", (m && m[0]) ? m : "", 1);
      unsetenv("OMP_NUM_THREADS");
      unsetenv("MKL_NUM_THREADS");
    }
    printf(" BDDC: relaunching under MPI:");
    for (i = 0; nargv[i] != NULL; i++) printf(" %s", nargv[i]);
    printf("\n");
    fflush(stdout);
    execvp(nargv[0], nargv);
    printf(" *ERROR: BDDC: relaunch failed (%s not found?); check CALTOP_MPIEXEC\n", nargv[0]);
    fflush(stdout);
    exit(1);
  }
  if (PetscInitializeNoArguments() != PETSC_SUCCESS) {
    printf(" *ERROR: BDDC: PETSc initialization failed\n");
    fflush(stdout);
    exit(1);
  }
  MPI_Comm_rank(PETSC_COMM_WORLD, &bd_rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &bd_size);
#ifdef PROFILING_ON
  TAU_PROFILE_SET_NODE(bd_rank);   /* profile.<rank>.0.<thread>: ranks do not overwrite each other */
#endif
  if (getd("CALTOP_BDDC_PROFILE", 0.) > 0.) PetscLogDefaultBegin();   /* per-rank setup report */
  if (bd_size == 1 && getd("CALTOP_BDDC_NP", 24.) > 1.) {
    /* every process thinks it is alone: the launcher belongs to a different
       MPI than the one calTop.exe was linked with (e.g. Intel MPI's mpiexec
       first in PATH, Open MPI in the build). Stop all copies at once. */
    printf(" *ERROR: BDDC: this process sees 1 MPI rank, CALTOP_BDDC_NP=%s. The launcher (CALTOP_MPIEXEC) is from a\n"
           "         different MPI than calTop.exe was built with; give its full path, e.g.\n"
           "         CALTOP_MPIEXEC=\"/usr/bin/mpiexec -n 24 --bind-to core\" (check: ldd calTop.exe | grep mpi)\n",
           getenv("CALTOP_BDDC_NP"));
    fflush(stdout);
    PetscFinalize();
    exit(1);
  }
  if (bd_rank > 0) {
    /* solver worker: never runs CalTop; thread counts left at their defaults */
    if (bd_rank == 1) {
      cpu_set_t cs;
      sched_getaffinity(0, sizeof(cs), &cs);
      char ver[256];
      mkl_get_version_string(ver, sizeof(ver));
      printf(" BDDC: solver ranks run %d OpenMP / %d MKL thread(s) on %d core(s) each\n BDDC: %s\n", omp_get_max_threads(),
             mkl_get_max_threads(), CPU_COUNT(&cs), ver);
      fflush(stdout);
    }
    worker_loop();
    PetscFinalize();
    exit(0);
  }
  /* rank 0 runs CalTop: give back the thread counts the user asked for */
  {
    char *o = getenv("CALTOP_BDDC_OMP0"), *m = getenv("CALTOP_BDDC_MKL0");
    if (o != NULL) {
      if (o[0]) { setenv("OMP_NUM_THREADS", o, 1); omp_set_num_threads(atoi(o)); }
      else unsetenv("OMP_NUM_THREADS");
    }
    if (m != NULL) {
      if (m[0]) { setenv("MKL_NUM_THREADS", m, 1); mkl_set_num_threads(atoi(m)); }
      else unsetenv("MKL_NUM_THREADS");
    }
  }
  /* rank 0 runs CalTop's OpenMP parts: free it from the launcher's core
     binding (the solver ranks stay pinned, one core each) */
  {
    cpu_set_t all;
    long      nc = sysconf(_SC_NPROCESSORS_ONLN), c;
    CPU_ZERO(&all);
    for (c = 0; c < nc && c < CPU_SETSIZE; c++) CPU_SET(c, &all);
    sched_getaffinity(0, sizeof(bd_core), &bd_core);   /* the core the launcher gave us */
    {
      int cc = -1, nn = CPU_COUNT(&bd_core);
      for (c = 0; c < CPU_SETSIZE; c++) if (CPU_ISSET(c, &bd_core)) { cc = (int)c; break; }
      printf(" BDDC: rank 0 solves on core %d%s; CalTop's own work uses all cores\n", cc,
             nn > 1 ? " (launcher gave it several: use --bind-to core)" : "");
    }
    sched_setaffinity(0, sizeof(all), &all);
  }
  setvbuf(stdout, NULL, _IOLBF, 0);   /* log stays complete if a later step crashes */
  {
    /* crash report with a backtrace (replaces PETSc's handler on rank 0):
       addresses can be resolved with addr2line -f -e calTop.exe <offset> */
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = crash_handler;
    sigaction(SIGSEGV, &sa, NULL);
    sigaction(SIGBUS, &sa, NULL);
    sigaction(SIGABRT, &sa, NULL);
    sigaction(SIGFPE, &sa, NULL);
  }
  bd_on = 1;
  atexit(bddc_finish);
  printf(" BDDC: %d MPI ranks (rank 0 runs CalTop and subdomain 0)\n", bd_size);
}

int bddc_active(void) { return bd_on && !bd_disabled; }

static void bd_fatal(const char *msg)
{
  printf(" *ERROR: BDDC: %s\n", msg);
  fflush(stdout);
  MPI_Abort(PETSC_COMM_WORLD, 1);
}

/* ---------------- capture (rank 0, called from mafillsm threads) ---------------- */

/* mafillsm.f calls ctelemcapture for every element matrix */
void ctelemcapture_(ITG *iel, ITG *nope, ITG *ndof, double *sm) { bddc_elem(*iel - 1, *nope, *ndof, sm); }

void bddc_capture_begin(ITG ne, ITG *ipkon, char *lakon)
{
  if (!bddc_active()) return;
  if (bd_dist == 1) {                /* element matrices are computed on the ranks (bddc_dist_asm) */
    bd_capturing = 1;
    return;
  }
  if (bd_K == NULL || bd_ne != ne) {
    /* element types: C3D4 and C3D10 (any mix). The slot size is the largest
       element: 12 unknowns for an all-C3D4 mesh, 30 as soon as there is a
       C3D10 (C3D4 elements then use the first 12 rows/columns, the rest is
       padding that is never assembled) */
    ITG  e, n4 = 0, n10 = 0;
    char msg[160];
    free(bd_K); free(bd_have); free(bd_enope);
    bd_enope = calloc(ne, 1);
    if (bd_enope == NULL) bd_fatal("no memory for the element list");
    for (e = 0; e < ne; e++) {
      const char *t = lakon + 8 * (size_t)e;
      if (ipkon[e] < 0) continue;
      if (strncmp(t, "C3D10", 5) == 0) { bd_enope[e] = 10; n10++; }
      else if (strncmp(t, "C3D4", 4) == 0 && (t[4] == ' ' || t[4] == '\0')) { bd_enope[e] = 4; n4++; }
      else {
        snprintf(msg, sizeof(msg), "element %d is %.8s; BDDC supports C3D4 and C3D10 only", (int)e + 1, t);
        bd_fatal(msg);
      }
    }
    if (n4 + n10 == 0) bd_fatal("no C3D4 or C3D10 elements in the model");
    bd_nn = (n10 > 0) ? 10 : 4;
    bd_s = 3 * bd_nn; bd_s2 = bd_s * bd_s;
    bd_K    = malloc(sizeof(double) * bd_s2 * (size_t)ne);
    bd_have = calloc(ne, 1);
    if (bd_K == NULL || bd_have == NULL) {
      snprintf(msg, sizeof(msg), "no memory for the element matrices (%.1f GB)", 8. * bd_s2 * (double)ne / 1e9);
      bd_fatal(msg);
    }
    bd_ne = ne;
    printf(" BDDC: elements: %d C3D10, %d C3D4; %dx%d matrix per element\n", (int)n10, (int)n4, bd_s, bd_s);
  } else memset(bd_have, 0, ne);
  bd_capturing = 1;
}

/* mafillsm / mafillsmmain: 1 = do not assemble CalTop's global K in this pass.
   A BDDC run never needs it: BDDC works from the element matrices, and there is
   no fallback solver. */
ITG ctnokmat_(void) { return (bddc_active() && bd_capturing) || bd_wcap0; }

/* s(60,60) column-major, upper triangle filled (mafillsm) */
void bddc_elem(ITG e, ITG nope, ITG ndof, double *sm)
{
  ITG a, b;
  double *K;
  ITG m;
  /* bd_K == NULL: distributed element work, rank 0's right-hand-side pass (no matrices kept) */
  if (!bd_capturing || bd_K == NULL || e < 0 || e >= bd_ne || ndof != 3 || nope != bd_enope[e]) return;
  m = 3 * nope;
  K = bd_K + (size_t)e * bd_s2;
  if (m < bd_s) memset(K, 0, sizeof(double) * bd_s2);   /* C3D4 in a mixed mesh: padding */
  for (a = 0; a < m; a++)
    for (b = a; b < m; b++) {
      double v = sm[(size_t)b * 60 + a];
      K[a * bd_s + b] = v; K[b * bd_s + a] = v;
    }
  bd_have[e] = 1;
}

/* ---------------- collective pieces ---------------- */

static void send_cmd(double c, double a1, double a2, double a3)
{
  double msg[4] = {c, a1, a2, a3};
  int    r;
  for (r = 1; r < bd_size; r++) MPI_Send(msg, 4, MPI_DOUBLE, r, TAG_CMD, PETSC_COMM_WORLD);
}

/* topology: rank 0 sends each rank its elements' equations; all build
   the local numbering, the MATIS matrices and the KSP shell */
/* pass the BDDC customization to a new PC: at the first set-up and again at
   every rebuild (pmat_collective). The graph stays CalTop's
   (PETSC_USE_POINTER), so it is passed again without a copy. */
static PetscErrorCode pc_customize(PC pc)
{
  IS isd;
  PetscFunctionBeginUser;
  PetscCall(PCBDDCSetLocalAdjacencyGraph(pc, nloc, bd_xa, bd_aj, PETSC_USE_POINTER));
  PetscCall(PCBDDCSetDofsSplittingLocal(pc, 3, bd_fld));
  PetscCall(PCSetCoordinates(pc, 3, m_own, bd_cr));
  /* Empty Dirichlet set (CalTop has already removed the fixed unknowns).
     Only telling BDDC that a Dirichlet set exists switches PETSc's
     adaptive Schur code to the path that reuses the MUMPS factor to get
     the inverse Schur complements (no MatFactorInvertSchurComplement);
     the standalone test always had it and set up 2.5x faster. */
  PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, 0, NULL, PETSC_COPY_VALUES, &isd));
  PetscCall(PCBDDCSetDirichletBoundaries(pc, isd));
  PetscCall(ISDestroy(&isd));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode topo_collective(int *myedof, int mynel, PetscInt n, double *coeq, PetscInt bs)
{
  PetscInt *tmp, cnt = 0, i, k;
  ISLocalToGlobalMapping map;
  PetscFunctionBeginUser;
  n_g = n;
  bs_g = bs;
  r_0 = bs * ((n / bs) * bd_rank / bd_size);
  m_own = bs * ((n / bs) * (bd_rank + 1) / bd_size) - r_0;
  nown = mynel; ledof = myedof;
  /* local dofs = equations of own elements */
  PetscCall(PetscMalloc1((size_t)bd_s * nown + 1, &tmp));
  for (i = 0; i < (PetscInt)bd_s * nown; i++) if (ledof[i] >= 0) tmp[cnt++] = ledof[i];
  qsort(tmp, cnt, sizeof(PetscInt), cmpint);
  nloc = 0;
  for (i = 0; i < cnt; i++) if (nloc == 0 || tmp[i] != tmp[nloc - 1]) tmp[nloc++] = tmp[i];
  lglob = tmp;
  /* element dofs -> local indices */
  for (i = 0; i < (PetscInt)bd_s * nown; i++)
    if (ledof[i] >= 0) {
      PetscInt g = ledof[i], *p = bsearch(&g, lglob, nloc, sizeof(PetscInt), cmpint);
      ledof[i] = (int)(p - lglob);
    }
  if (bs == 3) {
    /* complete nodes: block indices; local index 3i+c <-> global 3*block_i+c */
    PetscInt *blk;
    PetscCall(PetscMalloc1(nloc / 3 + 1, &blk));
    for (i = 0; i < nloc / 3; i++) blk[i] = lglob[3 * i] / 3;
    PetscCall(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 3, nloc / 3, blk, PETSC_OWN_POINTER, &map));
  } else PetscCall(ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, nloc, lglob, PETSC_COPY_VALUES, &map));
  PetscCall(MatCreateIS(PETSC_COMM_WORLD, bs, m_own, m_own, n, n, map, map, &A));
  PetscCall(ISLocalToGlobalMappingDestroy(&map));
  /* local dof graph: preallocation and BDDC adjacency */
  {
    PetscInt *c2, *ptr, *lst, *mark, *nnz, *xa, *aj, a, c;
    PetscCall(PetscCalloc1(nloc + 1, &c2));
    for (k = 0; k < nown; k++)
      for (i = 0; i < bd_s; i++) if (ledof[bd_s * k + i] >= 0) c2[ledof[bd_s * k + i] + 1]++;
    PetscCall(PetscMalloc1(nloc + 1, &ptr));
    ptr[0] = 0;
    for (a = 0; a < nloc; a++) ptr[a + 1] = ptr[a] + c2[a + 1];
    PetscCall(PetscMalloc1(ptr[nloc] + 1, &lst));
    for (a = 0; a < nloc; a++) c2[a] = ptr[a];
    for (k = 0; k < nown; k++)
      for (i = 0; i < bd_s; i++) if (ledof[bd_s * k + i] >= 0) lst[c2[ledof[bd_s * k + i]]++] = k;
    PetscCall(PetscMalloc1(nloc, &mark));
    PetscCall(PetscMalloc1(nloc + 1, &nnz));
    PetscCall(PetscMalloc1(nloc + 1, &xa));
    for (a = 0; a < nloc; a++) mark[a] = -1;
    xa[0] = 0;
    for (a = 0; a < nloc; a++) {
      PetscInt d = 0;
      for (c = ptr[a]; c < ptr[a + 1]; c++)
        for (i = 0; i < bd_s; i++) {
          PetscInt b2 = ledof[bd_s * lst[c] + i];
          if (b2 >= 0 && mark[b2] != a) { mark[b2] = a; d++; }
        }
      nnz[a] = d; xa[a + 1] = xa[a] + d;
    }
    PetscCall(PetscMalloc1(xa[nloc] + 1, &aj));
    for (a = 0; a < nloc; a++) mark[a] = -1;
    for (a = 0; a < nloc; a++) {
      PetscInt d = xa[a];
      for (c = ptr[a]; c < ptr[a + 1]; c++)
        for (i = 0; i < bd_s; i++) {
          PetscInt b2 = ledof[bd_s * lst[c] + i];
          if (b2 >= 0 && mark[b2] != a) { mark[b2] = a; aj[d++] = b2; }
        }
      qsort(aj + xa[a], xa[a + 1] - xa[a], sizeof(PetscInt), cmpint);
    }
    {
      Mat Al;
      PetscCall(MatISGetLocalMat(A, &Al));
      PetscCall(MatSeqAIJSetPreallocation(Al, 0, nnz));
      if (bs == 3) {
        PetscInt *bnz;
        PetscCall(PetscMalloc1(nloc / 3 + 1, &bnz));
        for (a = 0; a < nloc / 3; a++) bnz[a] = nnz[3 * a] / 3;
        PetscCall(MatSeqBAIJSetPreallocation(Al, 3, 0, bnz));
        PetscCall(PetscFree(bnz));
      }
      PetscCall(MatSetOption(Al, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
      PetscCall(MatISRestoreLocalMat(A, &Al));
    }
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, KSPCG));
    {
      PC        pc;
      PetscInt *fi[3], fn[3] = {0, 0, 0}, f;
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCBDDC));
      bd_xa = xa; bd_aj = aj;
      /* x / y / z unknowns as separate fields (constraints per component) */
      for (f = 0; f < 3; f++) PetscCall(PetscMalloc1(nloc + 1, &fi[f]));
      for (a = 0; a < nloc; a++) { f = (PetscInt)lrint(coeq[4 * lglob[a] + 3]); fi[f][fn[f]++] = a; }
      for (f = 0; f < 3; f++) PetscCall(ISCreateGeneral(PETSC_COMM_SELF, fn[f], fi[f], PETSC_OWN_POINTER, &bd_fld[f]));
      /* coordinates of the owned equations (corner selection) */
      PetscCall(PetscMalloc1(3 * m_own + 1, &bd_cr));
      for (i = 0; i < m_own; i++)
        for (k = 0; k < 3; k++) bd_cr[3 * i + k] = coeq[4 * (r_0 + i) + k];
      PetscCall(pc_customize(pc));
      /* rigid-body modes from the coordinates (as the standalone test) */
      if (bs == 3) {
        PetscScalar *xc;
        PetscCall(VecCreateMPI(PETSC_COMM_WORLD, m_own, n, &bd_coords));
        PetscCall(VecSetBlockSize(bd_coords, 3));
        PetscCall(VecGetArray(bd_coords, &xc));
        for (i = 0; i < m_own; i++) xc[i] = coeq[4 * (r_0 + i) + i % 3];
        PetscCall(VecRestoreArray(bd_coords, &xc));
      }
    }
    PetscCall(PetscFree(c2)); PetscCall(PetscFree(ptr)); PetscCall(PetscFree(lst));
    PetscCall(PetscFree(mark)); PetscCall(PetscFree(nnz));
  }
  PetscCall(MatCreateVecs(A, &vx, &vb));
  PetscCall(VecDuplicate(vb, &vr));
  PetscCall(PetscMalloc1((size_t)bd_s2 * nown + 1, &lK));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static void set_default(const char *name, const char *val)
{
  PetscBool has;
  PetscOptionsHasName(NULL, NULL, name, &has);
  if (!has) PetscOptionsSetValue(NULL, name, val);
}

/* resident memory over the ranks, printed by rank 0 (collective) */
static PetscErrorCode report_memory(const char *when)
{
  PetscLogDouble m, msum = 0., mmax = 0.;
  PetscFunctionBeginUser;
  PetscCall(PetscMemoryGetCurrentUsage(&m));
  PetscCallMPI(MPI_Reduce(&m, &msum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&m, &mmax, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  if (bd_rank == 0)
    printf(" BDDC: memory after %s: %.1f GB over all ranks, %.1f GB largest rank, rank 0 (CalTop + subdomain 0) %.1f GB\n", when, msum / 1e9,
           mmax / 1e9, m / 1e9);
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* new element matrices (in lK): assemble A; the first one also becomes
   the preconditioner matrix B and BDDC is set up from it */
static PetscErrorCode matrix_collective(double *tsetup)
{
  PetscInt li[30], c, j;
  PetscFunctionBeginUser;
  PetscCall(MatZeroEntries(A));
  for (c = 0; c < nown; c++) {
    for (j = 0; j < bd_s; j++) li[j] = ledof[bd_s * c + j];
    PetscCall(MatSetValuesLocal(A, bd_s, li, bd_s, li, lK + (size_t)bd_s2 * c, ADD_VALUES));
  }
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  if (!pc_ready) {
    TAU_REGION_START("Forward analysis: Preconditioner construction");
    double t0;
    char   thr[32];
    /* flags exactly as the standalone test: A marked SPD before it is
       copied, so the subdomain matrices of B are known to be SPD too (MUMPS
       then uses its positive definite factorization, and the Schur
       complements need no explicit inverse) */
    PetscCall(MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE));
    PetscCall(MatSetOption(A, MAT_SPD, PETSC_TRUE));
    PetscCall(MatDuplicate(A, MAT_COPY_VALUES, &B));
    PetscCall(MatSetOption(B, MAT_SYMMETRIC, PETSC_TRUE));
    PetscCall(MatSetOption(B, MAT_SPD, PETSC_TRUE));
    {
      Mat Bl;
      PetscCall(MatISGetLocalMat(B, &Bl));
      PetscCall(MatSetOption(Bl, MAT_SYMMETRIC, PETSC_TRUE));
      PetscCall(MatSetOption(Bl, MAT_SPD, PETSC_TRUE));
      PetscCall(MatISRestoreLocalMat(B, &Bl));
    }
    if (bd_coords) {
      MatNullSpace nsp;
      PetscCall(MatNullSpaceCreateRigidBody(bd_coords, &nsp));
      PetscCall(MatSetNearNullSpace(A, nsp));
      PetscCall(MatSetNearNullSpace(B, nsp));
      PetscCall(MatNullSpaceDestroy(&nsp));
    }
    snprintf(thr, sizeof(thr), "%g", 5.);
    set_default("-pc_bddc_use_faces", "");
    set_default("-pc_bddc_corner_selection", "");
    set_default("-pc_bddc_use_deluxe_scaling", "");
    set_default("-pc_bddc_adaptive_threshold", thr);
#if defined(PETSC_HAVE_MUMPS)
    set_default("-pc_bddc_dirichlet_pc_factor_mat_solver_type", "mumps");
    set_default("-pc_bddc_neumann_pc_factor_mat_solver_type", "mumps");
#endif
    {
      char *o = getenv("CALTOP_BDDC_OPTS");
      if (o != NULL && o[0]) PetscCall(PetscOptionsInsertString(NULL, o));
    }
    PetscCall(KSPSetOperators(ksp, A, B));
    PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    t0 = MPI_Wtime();
    PetscCall(KSPSetUp(ksp));
    {
      PC pc;
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetUp(pc));
    }
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (tsetup) *tsetup = MPI_Wtime() - t0;
    PetscCall(KSPSetReusePreconditioner(ksp, PETSC_TRUE));
    if (getd("CALTOP_BDDC_PROFILE", 0.) > 0.) {
      /* per rank: subdomain size and the time of the main setup steps */
      const char *ev[4] = {"MatCholFctrNum", "MatMatSolve", "MatFctrInvSchur", "PCBDDCSchr"};
      double      mine[8], *all = NULL;
      int         r, e;
      cpu_set_t   cs;
      mine[0] = (double)nloc;
      {
        PC pcx;
        PetscCall(KSPGetPC(ksp, &pcx));
        mine[6] = (double)((PC_IS *)pcx->data)->n_B;
        mine[7] = (double)((PC_BDDC *)pcx->data)->local_primal_size;
        if (bd_rank == 0) printf(" BDDC: coarse size %d\n", (int)((PC_BDDC *)pcx->data)->coarse_size);
      }
      for (e = 0; e < 4; e++) {
        PetscLogEvent      id;
        PetscEventPerfInfo info;
        mine[1 + e] = -1.;
        if (PetscLogEventGetId(ev[e], &id) == PETSC_SUCCESS && PetscLogEventGetPerfInfo(0, id, &info) == PETSC_SUCCESS) mine[1 + e] = info.time;
      }
      sched_getaffinity(0, sizeof(cs), &cs);
      mine[5] = -1.;
      for (r = 0; r < CPU_SETSIZE; r++) if (CPU_ISSET(r, &cs)) { mine[5] = r; break; }
      if (bd_rank == 0) all = malloc(sizeof(double) * 8 * bd_size);
      PetscCallMPI(MPI_Gather(mine, 8, MPI_DOUBLE, all, 8, MPI_DOUBLE, 0, PETSC_COMM_WORLD));
      if (bd_rank == 0) {
        printf(" BDDC profile per rank (seconds):\n  rank  core  local unknowns  interface  primal  Cholesky  MatMatSolve  InvSchur  SchurTotal\n");
        for (r = 0; r < bd_size; r++)
          printf("  %4d  %4d  %14.0f  %9.0f  %6.0f  %8.1f  %11.1f  %8.1f  %10.1f\n", r, (int)all[8 * r + 5], all[8 * r], all[8 * r + 6], all[8 * r + 7],
                 all[8 * r + 1], all[8 * r + 2], all[8 * r + 3], all[8 * r + 4]);
        free(all);
      }
    }
    PetscCall(report_memory("setup"));
    pc_ready = 1;
    TAU_REGION_STOP("Forward analysis: Preconditioner construction");
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* preconditioner rebuild: element matrices of the zero-stress assembly at
   the current configuration (in lK) become B; BDDC is set up again on the
   same subdomains (same topology, new values) */
static PetscErrorCode pmat_collective(double *tsetup)
{
  PetscInt li[30], c, j;
  Mat      Bl;
  PC       pc;
  double   t0;
  PetscFunctionBeginUser;
  TAU_REGION_START("Forward analysis: Preconditioner construction");
  PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
  t0 = MPI_Wtime();
  /* free the old preconditioner first, so the old and the new one are never
     in memory together: the KSP gets a new, empty PC, which releases the
     old one. The new PC gets the same customization and options; its
     PCSetUp analyses the subdomain topology again (same result, small
     extra cost). */
  {
    PC pcn;
    PetscCall(PCCreate(PETSC_COMM_WORLD, &pcn));
    PetscCall(KSPSetPC(ksp, pcn));  /* the old PC is destroyed here */
    PetscCall(PCDestroy(&pcn));     /* the KSP keeps its reference */
  }
#ifdef __GLIBC__
  malloc_trim(0);                   /* return the freed pages to the system */
#endif
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCBDDC));
  PetscCall(pc_customize(pc));
  PetscCall(MatZeroEntries(B));
  for (c = 0; c < nown; c++) {
    for (j = 0; j < bd_s; j++) li[j] = ledof[bd_s * c + j];
    PetscCall(MatSetValuesLocal(B, bd_s, li, bd_s, li, lK + (size_t)bd_s2 * c, ADD_VALUES));
  }
  PetscCall(MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY));
  PetscCall(MatSetOption(B, MAT_SYMMETRIC, PETSC_TRUE));
  PetscCall(MatSetOption(B, MAT_SPD, PETSC_TRUE));
  PetscCall(MatISGetLocalMat(B, &Bl));
  PetscCall(MatSetOption(Bl, MAT_SYMMETRIC, PETSC_TRUE));
  PetscCall(MatSetOption(Bl, MAT_SPD, PETSC_TRUE));
  PetscCall(MatISRestoreLocalMat(B, &Bl));
  PetscCall(KSPSetOperators(ksp, A, B));
  PetscCall(PCSetFromOptions(pc));
  PetscCall(KSPSetReusePreconditioner(ksp, PETSC_FALSE));
  PetscCall(PCSetUp(pc));
  PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
  if (tsetup) *tsetup = MPI_Wtime() - t0;
  PetscCall(KSPSetReusePreconditioner(ksp, PETSC_TRUE));
  PetscCall(report_memory("rebuild"));
  TAU_REGION_STOP("Forward analysis: Preconditioner construction");
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* rhs owned slices are in vb; result in vx */
/* CG from a zero start; the true relative residual is returned in res */
static PetscErrorCode solve_collective(double tol, int maxit, PetscInt *its, double *res, int *reason)
{
  KSPConvergedReason r;
  PetscReal          rn, bn;
  PetscFunctionBeginUser;
  PetscCall(KSPSetTolerances(ksp, tol, 0., 1.e30, maxit));
  PetscCall(VecSet(vx, 0.));
  PetscCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));
  PetscCall(KSPSolve(ksp, vb, vx));
  PetscCall(KSPGetIterationNumber(ksp, its));
  PetscCall(KSPGetConvergedReason(ksp, &r));
  PetscCall(MatMult(A, vx, vr));
  PetscCall(VecAYPX(vr, -1., vb));
  PetscCall(VecNorm(vr, NORM_2, &rn));
  PetscCall(VecNorm(vb, NORM_2, &bn));

  *res = bn > 0. ? rn / bn : 0.;
  *reason = (int)r;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ---------------- distributed element work (assembly, internal forces) ----------------
   Every rank computes the element stiffness and the stresses/internal forces of
   the elements of its own subdomain. Rank 0 keeps all of CalTop's arrays and
   logic: for each call it sends each rank the displacements of its nodes and the
   element state of its elements, and takes back what was computed. Between calls
   a rank keeps only its read-only local model (its elements, their nodes with
   local numbers, the materials). Element and node numbers in CalTop, the files
   and the optimisation are not changed. */

#if defined(LONGLONG)
#define MPI_ITG MPI_LONG_LONG
#else
#define MPI_ITG MPI_INT
#endif
#define TAG_D 7200

typedef struct {
  ITG     ne, nk, nkon, mi[3], mt, nmat, ntmat_, ncmat_, npmat_, norien;
  ITG    *kon, *ipkon, *ielmat, *ielorien, *ielprop, *nactdof, *nelcon, *nrhcon, *nalcon, *nplicon, *nplkcon;
  unsigned char *enope;
  char   *have;
  char   *lakon, *matname;
  double *co, *design, *thicke, *elcon, *rhcon, *alcon, *alzero, *plicon, *plkcon, *orab, physcon[13], prop[2];
  double  penal, matdens;
  double *v, *vold, *zn, *fn, *sti, *stx, *xstiff, *eei, *eme, *emeini, *stiini, *ener, *enerini;
} dmodel_t;
static dmodel_t dm;
static int      bd_dist = 0;          /* rank 0: 1 element work distributed, -1 stays on rank 0, 0 not decided */
static int    **bd_elst = NULL, *bd_ecnt = NULL;   /* rank 0: elements of each rank (global, sorted) */
static int    **bd_nlst = NULL, *bd_ncnt = NULL;   /* rank 0: nodes of each rank (global, sorted) */
static MPI_Datatype *bd_tnode = NULL, *bd_te6 = NULL, *bd_te27 = NULL, *bd_te1 = NULL;
static double   bd_tdasm = 0., bd_tdres = 0.;
static int      bd_ndasm = 0, bd_ndres = 0;
static ITG      d_idum[512];
static double   d_ddum[512];
static char     d_cdum[1024];

/* sizes of the material arrays (same on every rank) */
#define D_NELCON (2 * dm.nmat)
#define D_ELCON  ((size_t)(dm.ncmat_ + 1) * dm.ntmat_ * dm.nmat)
#define D_RHCON  ((size_t)2 * dm.ntmat_ * dm.nmat)
#define D_ALCON  ((size_t)7 * dm.ntmat_ * dm.nmat)
#define D_PLICON ((size_t)(2 * dm.npmat_ + 1) * dm.ntmat_ * dm.nmat)
#define D_NPLICON ((size_t)(dm.ntmat_ + 1) * dm.nmat)

static void d_alloc_work(void)
{
  size_t nn = (size_t)dm.mt * dm.nk, e6 = (size_t)6 * dm.mi[0] * dm.ne;
  dm.v = calloc(nn + 1, sizeof(double)); dm.vold = calloc(nn + 1, sizeof(double));
  dm.zn = calloc(nn + 1, sizeof(double)); dm.fn = calloc(nn + 1, sizeof(double));
  dm.sti = calloc(e6 + 1, sizeof(double)); dm.stx = calloc(e6 + 1, sizeof(double));
  dm.eei = calloc(e6 + 1, sizeof(double)); dm.eme = calloc(e6 + 1, sizeof(double));
  dm.emeini = calloc(e6 + 1, sizeof(double)); dm.stiini = calloc(e6 + 1, sizeof(double));
  dm.xstiff = calloc((size_t)27 * dm.mi[0] * dm.ne + 1, sizeof(double));
  dm.ener = calloc((size_t)2 * dm.mi[0] * dm.ne + 1, sizeof(double));
  dm.enerini = calloc((size_t)2 * dm.mi[0] * dm.ne + 1, sizeof(double));
  dm.have = calloc(dm.ne + 1, 1);
  if (!dm.v || !dm.vold || !dm.zn || !dm.fn || !dm.sti || !dm.stx || !dm.eei || !dm.eme || !dm.emeini || !dm.stiini || !dm.xstiff ||
      !dm.ener || !dm.enerini || !dm.have)
    bd_fatal("no memory for the local element model");
}

/* model of rank r (rank 0): elements bd_elst[r], their nodes renumbered locally */
static void d_build_local(bddc_eargs *a, int r, int nmat, dmodel_t *o)
{
  ITG  mt = a->mi[1] + 1, m3 = a->mi[2], k, j, c, q;
  int *el = bd_elst[r], ne_l = bd_ecnt[r], *nl = bd_nlst[r], nk_l = bd_ncnt[r];
  int *g2l = malloc(sizeof(int) * (size_t)(*a->nk));
  for (k = 0; k < *a->nk; k++) g2l[k] = -1;
  for (k = 0; k < nk_l; k++) g2l[nl[k]] = (int)k;
  memset(o, 0, sizeof(*o));
  o->ne = ne_l; o->nk = nk_l; o->mi[0] = a->mi[0]; o->mi[1] = a->mi[1]; o->mi[2] = a->mi[2]; o->mt = mt;
  o->nmat = nmat; o->ntmat_ = *a->ntmat_; o->ncmat_ = *a->ncmat_; o->npmat_ = *a->npmat_; o->norien = *a->norien;
  o->nkon = 0;
  for (k = 0; k < ne_l; k++) o->nkon += bd_enope[el[k]];
  o->kon = malloc(sizeof(ITG) * (o->nkon + 1)); o->ipkon = malloc(sizeof(ITG) * (ne_l + 1));
  o->lakon = malloc(8 * (size_t)ne_l + 1); o->ielmat = malloc(sizeof(ITG) * (size_t)m3 * ne_l + 1);
  o->ielorien = malloc(sizeof(ITG) * (size_t)m3 * ne_l + 1); o->ielprop = malloc(sizeof(ITG) * (ne_l + 1));
  o->design = malloc(sizeof(double) * (ne_l + 1)); o->thicke = calloc((size_t)m3 * o->nkon + 1, sizeof(double));
  o->co = malloc(sizeof(double) * 3 * (size_t)nk_l + 1); o->nactdof = malloc(sizeof(ITG) * (size_t)mt * nk_l + 1);
  for (k = 0, q = 0; k < ne_l; k++) {
    int e = el[k], n = bd_enope[e];
    o->ipkon[k] = q;
    for (j = 0; j < n; j++) {
      o->kon[q + j] = g2l[a->kon[a->ipkon[e] + j] - 1] + 1;
      if (a->thicke != NULL)
        for (c = 0; c < m3; c++) o->thicke[(size_t)m3 * (q + j) + c] = a->thicke[(size_t)m3 * (a->ipkon[e] + j) + c];
    }
    q += n;
    memcpy(o->lakon + 8 * (size_t)k, a->lakon + 8 * (size_t)e, 8);
    for (c = 0; c < m3; c++) {
      o->ielmat[(size_t)m3 * k + c] = a->ielmat[(size_t)m3 * e + c];
      /* norien == 0: CalTop has freed ielorien (SFREE does not reset the pointer) */
      o->ielorien[(size_t)m3 * k + c] = (*a->norien > 0 && a->ielorien != NULL) ? a->ielorien[(size_t)m3 * e + c] : 0;
    }
    o->ielprop[k] = -1;
    o->design[k] = a->design[e];
  }
  for (k = 0; k < nk_l; k++) {
    for (c = 0; c < 3; c++) o->co[3 * (size_t)k + c] = a->co[3 * (size_t)nl[k] + c];
    for (c = 0; c < mt; c++) o->nactdof[(size_t)mt * k + c] = a->nactdof[(size_t)mt * nl[k] + c];
  }
  free(g2l);
}

static void d_send_local(dmodel_t *o, int r)
{
  ITG hdr[12] = {o->ne, o->nk, o->nkon, o->mi[0], o->mi[1], o->mi[2], o->nmat, o->ntmat_, o->ncmat_, o->npmat_, o->norien, o->mt};
  MPI_Send(hdr, 12, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->kon, o->nkon, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->ipkon, o->ne, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->lakon, 8 * o->ne, MPI_CHAR, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->ielmat, o->mi[2] * o->ne, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->ielorien, o->mi[2] * o->ne, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->design, o->ne, MPI_DOUBLE, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->thicke, o->mi[2] * o->nkon, MPI_DOUBLE, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->co, 3 * o->nk, MPI_DOUBLE, r, TAG_D, PETSC_COMM_WORLD);
  MPI_Send(o->nactdof, o->mt * o->nk, MPI_ITG, r, TAG_D, PETSC_COMM_WORLD);
}

static void d_recv_local(void)
{
  ITG        hdr[12];
  MPI_Status st;
  memset(&dm, 0, sizeof(dm));
  MPI_Recv(hdr, 12, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  dm.ne = hdr[0]; dm.nk = hdr[1]; dm.nkon = hdr[2]; dm.mi[0] = hdr[3]; dm.mi[1] = hdr[4]; dm.mi[2] = hdr[5];
  dm.nmat = hdr[6]; dm.ntmat_ = hdr[7]; dm.ncmat_ = hdr[8]; dm.npmat_ = hdr[9]; dm.norien = hdr[10]; dm.mt = hdr[11];
  dm.kon = malloc(sizeof(ITG) * (dm.nkon + 1)); dm.ipkon = malloc(sizeof(ITG) * (dm.ne + 1));
  dm.lakon = malloc(8 * (size_t)dm.ne + 1); dm.ielmat = malloc(sizeof(ITG) * (size_t)dm.mi[2] * dm.ne + 1);
  dm.ielorien = malloc(sizeof(ITG) * (size_t)dm.mi[2] * dm.ne + 1); dm.ielprop = malloc(sizeof(ITG) * (dm.ne + 1));
  dm.design = malloc(sizeof(double) * (dm.ne + 1)); dm.thicke = malloc(sizeof(double) * ((size_t)dm.mi[2] * dm.nkon + 1));
  dm.co = malloc(sizeof(double) * (3 * (size_t)dm.nk + 1)); dm.nactdof = malloc(sizeof(ITG) * ((size_t)dm.mt * dm.nk + 1));
  MPI_Recv(dm.kon, dm.nkon, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.ipkon, dm.ne, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.lakon, 8 * dm.ne, MPI_CHAR, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.ielmat, dm.mi[2] * dm.ne, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.ielorien, dm.mi[2] * dm.ne, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.design, dm.ne, MPI_DOUBLE, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.thicke, dm.mi[2] * dm.nkon, MPI_DOUBLE, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.co, 3 * dm.nk, MPI_DOUBLE, 0, TAG_D, PETSC_COMM_WORLD, &st);
  MPI_Recv(dm.nactdof, dm.mt * dm.nk, MPI_ITG, 0, TAG_D, PETSC_COMM_WORLD, &st);
  { ITG k; for (k = 0; k < dm.ne; k++) dm.ielprop[k] = -1; }
}

/* materials and scalars of the model: broadcast (rank 0 fills dm before) */
static void d_bcast_materials(void)
{
  int me = bd_rank;
  if (me != 0) {
    dm.nelcon = malloc(sizeof(ITG) * (D_NELCON + 1)); dm.elcon = malloc(sizeof(double) * (D_ELCON + 1));
    dm.rhcon = malloc(sizeof(double) * (D_RHCON + 1)); dm.nrhcon = malloc(sizeof(ITG) * (dm.nmat + 1));
    dm.alcon = malloc(sizeof(double) * (D_ALCON + 1)); dm.nalcon = malloc(sizeof(ITG) * (D_NELCON + 1));
    dm.alzero = malloc(sizeof(double) * (dm.nmat + 1)); dm.plicon = malloc(sizeof(double) * (D_PLICON + 1));
    dm.nplicon = malloc(sizeof(ITG) * (D_NPLICON + 1)); dm.plkcon = malloc(sizeof(double) * (D_PLICON + 1));
    dm.nplkcon = malloc(sizeof(ITG) * (D_NPLICON + 1)); dm.matname = malloc(80 * (size_t)dm.nmat + 1);
    dm.orab = malloc(sizeof(double) * (7 * (size_t)dm.norien + 1));
  }
  MPI_Bcast(dm.nelcon, D_NELCON, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.elcon, (int)D_ELCON, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.rhcon, (int)D_RHCON, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.nrhcon, dm.nmat, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.alcon, (int)D_ALCON, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.nalcon, D_NELCON, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.alzero, dm.nmat, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.plicon, (int)D_PLICON, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.nplicon, (int)D_NPLICON, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.plkcon, (int)D_PLICON, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.nplkcon, (int)D_NPLICON, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.matname, 80 * dm.nmat, MPI_CHAR, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.orab, 7 * dm.norien, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dm.physcon, 13, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&dm.penal, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&dm.matdens, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  dm.prop[0] = dm.prop[1] = 0.;
  {
    ITG k;
    dm.enope = malloc(dm.ne + 1);
    for (k = 0; k < dm.ne; k++) dm.enope[k] = (strncmp(dm.lakon + 8 * (size_t)k, "C3D10", 5) == 0) ? 10 : 4;
  }
  if (lK == NULL || (ITG)nown != dm.ne) bd_fatal("local element model does not match the subdomain");
  d_alloc_work();
}

/* can the element work be distributed? (rank 0) */
static int d_eligible(bddc_eargs *a, char *why, size_t nw, int *nmat)
{
  ITG e, c, m3 = a->mi[2];
  *nmat = 0;
  if (bd_size < 2) { snprintf(why, nw, "one rank"); return 0; }
  if (a->ithermal[0] != 0) { snprintf(why, nw, "thermal analysis"); return 0; }
  if (*a->iprestr != 0) { snprintf(why, nw, "initial stresses or strains"); return 0; }
  if (*a->nstate_ != 0) { snprintf(why, nw, "material state variables"); return 0; }
  if (*a->mortar != 0 || *a->nasym != 0 || *a->ne0 != *a->ne) { snprintf(why, nw, "contact elements"); return 0; }
  if (*a->nmethod != 1 || a->mass[0] != 0 || a->mass[1] != 0 || *a->buckling != 0) { snprintf(why, nw, "not a static analysis"); return 0; }
  for (e = 0; e < *a->ne; e++) {
    if (a->ipkon[e] < 0) continue;
    if (a->ielprop != NULL && a->ielprop[e] >= 0) { snprintf(why, nw, "element properties"); return 0; }
    for (c = 0; c < m3; c++) {
      ITG im = a->ielmat[(size_t)m3 * e + c];
      if (im < 0) { snprintf(why, nw, "deactivated elements"); return 0; }
      if (im > *nmat) *nmat = (int)im;
    }
  }
  for (c = 0; c < *nmat; c++) {
    if (a->nelcon[2 * c] <= 0) { snprintf(why, nw, "material %d is not linear elastic", (int)c + 1); return 0; }
    if (*a->npmat_ > 0 && a->nplicon != NULL && a->nplkcon != NULL && (a->nplicon[(size_t)(*a->ntmat_ + 1) * c] != 0 || a->nplkcon[(size_t)(*a->ntmat_ + 1) * c] != 0)) {
      snprintf(why, nw, "material %d is plastic", (int)c + 1); return 0;
    }
  }
  return 1;
}

/* rank 0: build and send the local models, datatypes (first distributed assembly) */
static int d_setup(bddc_eargs *a)
{
  char why[160];
  int  nmat, r, k, j;
  ITG  mt = a->mi[1] + 1;
  double t0 = MPI_Wtime();
  if (!d_eligible(a, why, sizeof(why), &nmat)) {
    printf(" BDDC: element work stays on rank 0 (%s)\n", why);
    bd_dist = -1;
    return 0;
  }
  /* nodes of each rank */
  bd_nlst = malloc(sizeof(int *) * bd_size); bd_ncnt = calloc(bd_size, sizeof(int));
  {
    int *mark = malloc(sizeof(int) * (size_t)(*a->nk));
    for (k = 0; k < *a->nk; k++) mark[k] = -1;
    for (r = 0; r < bd_size; r++) {
      int n = 0, *lst;
      for (k = 0; k < bd_ecnt[r]; k++) {
        int e = bd_elst[r][k];
        for (j = 0; j < bd_enope[e]; j++) {
          int nd = a->kon[a->ipkon[e] + j] - 1;
          if (mark[nd] != r) { mark[nd] = r; n++; }
        }
      }
      lst = malloc(sizeof(int) * (n + 1)); n = 0;
      for (k = 0; k < *a->nk; k++) if (mark[k] == r) lst[n++] = (int)k;   /* sorted */
      bd_nlst[r] = lst; bd_ncnt[r] = n;
    }
    free(mark);
  }
  /* datatypes: node blocks of mt doubles, element blocks of 6, 27 and 1 x mi[0] doubles */
  bd_tnode = malloc(sizeof(MPI_Datatype) * bd_size); bd_te6 = malloc(sizeof(MPI_Datatype) * bd_size);
  bd_te27 = malloc(sizeof(MPI_Datatype) * bd_size); bd_te1 = malloc(sizeof(MPI_Datatype) * bd_size);
  {
    MPI_Datatype bn, b6, b27, b1;
    MPI_Type_contiguous((int)mt, MPI_DOUBLE, &bn); MPI_Type_contiguous(6 * a->mi[0], MPI_DOUBLE, &b6);
    MPI_Type_contiguous(27 * a->mi[0], MPI_DOUBLE, &b27); MPI_Type_contiguous(a->mi[0], MPI_DOUBLE, &b1);
    for (r = 0; r < bd_size; r++) {
      MPI_Type_create_indexed_block(bd_ncnt[r], 1, bd_nlst[r], bn, &bd_tnode[r]); MPI_Type_commit(&bd_tnode[r]);
      MPI_Type_create_indexed_block(bd_ecnt[r], 1, bd_elst[r], b6, &bd_te6[r]); MPI_Type_commit(&bd_te6[r]);
      MPI_Type_create_indexed_block(bd_ecnt[r], 1, bd_elst[r], b27, &bd_te27[r]); MPI_Type_commit(&bd_te27[r]);
      MPI_Type_create_indexed_block(bd_ecnt[r], 1, bd_elst[r], b1, &bd_te1[r]); MPI_Type_commit(&bd_te1[r]);
    }
    MPI_Type_free(&bn); MPI_Type_free(&b6); MPI_Type_free(&b27); MPI_Type_free(&b1);
  }
  /* local models */
  send_cmd(CMD_DMODEL, 0, 0, 0);
  for (r = 1; r < bd_size; r++) {
    dmodel_t o;
    d_build_local(a, r, nmat, &o);
    d_send_local(&o, r);
    free(o.kon); free(o.ipkon); free(o.lakon); free(o.ielmat); free(o.ielorien); free(o.ielprop);
    free(o.design); free(o.thicke); free(o.co); free(o.nactdof);
  }
  d_build_local(a, 0, nmat, &dm);
  /* rank 0 broadcasts its own material arrays (zeros for any it does not have) */
#define D_OR0(p, n, T) ((p) != NULL ? (p) : (T *)calloc((size_t)(n) + 1, sizeof(T)))
  dm.nelcon = a->nelcon; dm.elcon = a->elcon;
  dm.rhcon = D_OR0(a->rhcon, D_RHCON, double); dm.nrhcon = D_OR0(a->nrhcon, dm.nmat, ITG);
  dm.alcon = D_OR0(a->alcon, D_ALCON, double); dm.nalcon = D_OR0(a->nalcon, D_NELCON, ITG);
  dm.alzero = D_OR0(a->alzero, dm.nmat, double); dm.plicon = D_OR0(a->plicon, D_PLICON, double);
  dm.nplicon = D_OR0(a->nplicon, D_NPLICON, ITG); dm.plkcon = D_OR0(a->plkcon, D_PLICON, double);
  dm.nplkcon = D_OR0(a->nplkcon, D_NPLICON, ITG); dm.matname = a->matname;
  dm.orab = (dm.norien > 0) ? D_OR0(a->orab, 7 * dm.norien, double) : calloc(1, sizeof(double));   /* freed when norien == 0 */
#undef D_OR0
  memcpy(dm.physcon, a->physcon, sizeof(dm.physcon));
  dm.penal = *a->penal; dm.matdens = *a->matdens;
  d_bcast_materials();
  /* the element matrices now stay on the ranks: rank 0's global store is not needed */
  free(bd_K); bd_K = NULL;
  bd_dist = 1;
  {
    long long nn = 0; for (r = 0; r < bd_size; r++) nn += bd_ncnt[r];
    printf(" BDDC: element work distributed: assembly and internal forces on %d ranks (%lld local nodes for %d nodes), %.1f s\n",
           bd_size, nn, (int)*a->nk, MPI_Wtime() - t0);
  }
  return 1;
}

/* scatter rank 0's global node/element arrays to the local ones (all ranks) */
static void d_scatter(double *g, MPI_Datatype *t, double *l, int nl, int tag)
{
  MPI_Request *rq = malloc(sizeof(MPI_Request) * (bd_size + 1));
  int          n = 0, r;
  TAU_REGION_START("Forward analysis: Data exchange");
  MPI_Irecv(l, nl, MPI_DOUBLE, 0, tag, PETSC_COMM_WORLD, &rq[n++]);
  if (bd_rank == 0)
    for (r = 0; r < bd_size; r++) MPI_Isend(g, 1, t[r], r, tag, PETSC_COMM_WORLD, &rq[n++]);
  MPI_Waitall(n, rq, MPI_STATUSES_IGNORE);
  free(rq);
  TAU_REGION_STOP("Forward analysis: Data exchange");
}
/* gather the local element arrays into rank 0's global array (disjoint elements) */
static void d_gather(double *g, MPI_Datatype *t, double *l, int nl, int tag)
{
  MPI_Request *rq = malloc(sizeof(MPI_Request) * (bd_size + 1));
  int          n = 0, r;
  TAU_REGION_START("Forward analysis: Data exchange");
  if (bd_rank == 0)
    for (r = 0; r < bd_size; r++) MPI_Irecv(g, 1, t[r], r, tag, PETSC_COMM_WORLD, &rq[n++]);
  MPI_Isend(l, nl, MPI_DOUBLE, 0, tag, PETSC_COMM_WORLD, &rq[n++]);
  MPI_Waitall(n, rq, MPI_STATUSES_IGNORE);
  free(rq);
  TAU_REGION_STOP("Forward analysis: Data exchange");
}

/* element matrices of this rank's elements into lK (all ranks) */
static void d_asm_local(ITG *ih, double *dh)
{
  ITG     nmeth = ih[0], iperturb[2] = {ih[1], ih[2]}, iexpl = ih[3], mass[2] = {0, 0}, stiff = 1, buck = 0, rhs0 = 0,
          intscheme = ih[4], istep = ih[5], iinc = ih[6], coriolis = ih[7], kscale = ih[8], network = ih[9],
          neq[2] = {ih[10], ih[11]}, nzs[3] = {ih[12], ih[13], ih[14]}, zero = 0, one = 1, ithermal[2] = {0, 0}, nzl = 0, nstate0 = 0,
          nasym = 0, mortar = 0, nelem = dm.ne, k;
  double  dtime = dh[0], ttime = dh[1], time = dh[2], reltime = dh[3];
  double *sK = bd_K;
  unsigned char *sE = bd_enope;
  char   *sH = bd_have;
  ITG     sN = bd_ne;
  int     sC = bd_capturing;
  memset(dm.have, 0, dm.ne);
  bd_K = lK; bd_enope = dm.enope; bd_have = dm.have; bd_ne = dm.ne; bd_capturing = 1; bd_wcap0 = 1;
  FORTRAN(mafillsm, (dm.co, &dm.nk, dm.kon, dm.ipkon, dm.lakon, &nelem, d_idum, d_idum, d_ddum, &zero, d_idum, d_idum, d_ddum, &zero,
                     d_idum, d_idum, d_ddum, &zero, d_idum, d_cdum, d_ddum, &zero, d_ddum, d_idum, &zero, d_ddum, d_ddum, d_ddum, d_ddum,
                     dm.nactdof, d_idum, d_idum, d_idum, neq, &nzl, &nmeth, d_idum, d_idum, d_idum, d_idum, dm.elcon, dm.nelcon, dm.rhcon,
                     dm.nrhcon, dm.alcon, dm.nalcon, dm.alzero, dm.ielmat, dm.ielorien, &dm.norien, dm.orab, &dm.ntmat_, dm.zn, dm.zn,
                     ithermal, d_ddum, &zero, dm.vold, iperturb, dm.sti, nzs, dm.sti, d_ddum, d_ddum, &iexpl, dm.plicon, dm.nplicon,
                     dm.plkcon, dm.nplkcon, dm.xstiff, &dm.npmat_, &dtime, dm.matname, dm.mi, &dm.ncmat_, mass, &stiff, &buck, &rhs0,
                     &intscheme, dm.physcon, d_ddum, d_idum, d_ddum, d_idum, &ttime, &time, &istep, &iinc, &coriolis, d_idum, d_ddum,
                     &reltime, dm.zn, d_ddum, &nstate0, d_ddum, d_ddum, dm.thicke, d_idum, d_ddum, d_cdum, d_idum, d_idum, d_idum, &zero,
                     &nasym, d_ddum, d_ddum, &mortar, d_ddum, dm.ielprop, dm.prop, &nelem, dm.zn, &one, &nelem, &kscale, d_idum, d_idum,
                     &network, dm.design, &dm.penal, &dm.matdens));
  bd_K = sK; bd_enope = sE; bd_have = sH; bd_ne = sN; bd_capturing = sC; bd_wcap0 = 0;
  for (k = 0; k < dm.ne; k++)
    if (!dm.have[k]) { char m[96]; snprintf(m, sizeof(m), "rank %d: element %d of its subdomain handed over no matrix", bd_rank, (int)k + 1); bd_fatal(m); }
  ih[0] = nmeth;
}

/* collective part of a distributed assembly (all ranks) */
static void d_asm_collective(bddc_eargs *a)
{
  ITG    ih[16];
  double dh[4];
  int    nmin, nme;
  if (bd_rank == 0) {
    ih[0] = *a->nmethod; ih[1] = a->iperturb[0]; ih[2] = a->iperturb[1]; ih[3] = *a->iexpl; ih[4] = *a->intscheme;
    ih[5] = *a->istep; ih[6] = *a->iinc; ih[7] = *a->coriolis; ih[8] = *a->kscale; ih[9] = *a->network;
    ih[10] = a->neq[0]; ih[11] = a->neq[1]; ih[12] = a->nzs[0]; ih[13] = a->nzs[1]; ih[14] = a->nzs[2]; ih[15] = 0;
    dh[0] = *a->dtime; dh[1] = *a->ttime; dh[2] = *a->time; dh[3] = *a->reltime;
  }
  MPI_Bcast(ih, 16, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dh, 4, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  d_scatter(bd_rank == 0 ? a->vold : NULL, bd_tnode, dm.vold, dm.mt * dm.nk, TAG_D + 1);
  d_scatter(bd_rank == 0 ? a->sti : NULL, bd_te6, dm.sti, 6 * dm.mi[0] * dm.ne, TAG_D + 2);
  /* the material tangent at the integration points (from results): e_c3d takes the
     elastic matrix from it (materialdata_me, istiff=1) */
  d_scatter(bd_rank == 0 ? a->xstiff : NULL, bd_te27, dm.xstiff, 27 * dm.mi[0] * dm.ne, TAG_D + 6);
  d_asm_local(ih, dh);
  nme = (int)ih[0];
  MPI_Reduce(&nme, &nmin, 1, MPI_INT, MPI_MIN, 0, PETSC_COMM_WORLD);
  if (bd_rank == 0 && nmin < *a->nmethod) *a->nmethod = nmin;   /* nonpositive Jacobian somewhere: nmethod=0 */
}

/* stresses, tangent and internal forces of this rank's elements (all ranks) */
static void d_res_local(ITG *ih, double *dh, double *qa, ITG *nal)
{
  ITG    iout = ih[0], nmeth = ih[1], iperturb[2] = {ih[2], ih[3]}, ielas = ih[4], icmd = ih[5], istep = ih[6], iinc = ih[7],
         cfn = ih[8], cqa = ih[9], cca = ih[10], nener = ih[11], ikin = ih[12], kscale = ih[13], zero = 0, one = 1,
         ithermal[2] = {0, 0}, nstate0 = 0, mortar = 0, nelem = dm.ne;
  double dtime = dh[0], time = dh[1], ttime = dh[2], reltime = dh[3], rhomin = dh[5];
  memset(dm.fn, 0, sizeof(double) * (size_t)dm.mt * dm.nk);
  qa[0] = qa[1] = qa[2] = qa[3] = 0.;
  *nal = 0;
  FORTRAN(resultsmech, (dm.co, dm.kon, dm.ipkon, dm.lakon, &nelem, dm.v, dm.stx, dm.elcon, dm.nelcon, dm.rhcon, dm.nrhcon, dm.alcon,
                        dm.nalcon, dm.alzero, dm.ielmat, dm.ielorien, &dm.norien, dm.orab, &dm.ntmat_, dm.zn, dm.zn, ithermal, d_ddum,
                        &zero, dm.eme, iperturb, dm.fn, &iout, qa, dm.vold, &nmeth, dm.zn, &dtime, &time, &ttime, dm.plicon,
                        dm.nplicon, dm.plkcon, dm.nplkcon, d_ddum, dm.xstiff, d_ddum, &dm.npmat_, dm.matname, dm.mi, &ielas, &icmd,
                        &dm.ncmat_, &nstate0, dm.stiini, dm.zn, dm.ener, dm.eei, dm.enerini, &istep, &iinc, d_ddum, &reltime, &cfn,
                        &cqa, &cca, &nener, &ikin, nal, &nelem, dm.thicke, dm.emeini, d_ddum, d_ddum, &mortar, d_ddum, &one, &nelem,
                        dm.ielprop, dm.prop, &kscale, &zero, d_idum, dm.design, &dm.penal, &rhomin));
}

/* element arrays of the internal-force call: bit k of ih[14] says that array k exists on rank 0 */
enum { DR_STX, DR_XSTIFF, DR_EEI, DR_EME, DR_EMEINI, DR_STIINI, DR_ENER, DR_ENERINI, DR_N };

static void d_res_collective(bddc_eargs *a)
{
  ITG     ih[16], nal = 0;
  double  dh[6], qa[4], *gq = NULL;
  int     me = bd_rank, have, w[DR_N] = {6, 27, 6, 6, 6, 6, 1, 1}, k, r;
  double *loc[DR_N] = {dm.stx, dm.xstiff, dm.eei, dm.eme, dm.emeini, dm.stiini, dm.ener, dm.enerini};
  double *glo[DR_N] = {NULL};
  int     out[DR_N] = {1, 1, 1, 1, 1, 0, 1, 0};
  if (me == 0) {
    double *g[DR_N] = {a->stx, a->xstiff, a->eei, a->eme, a->emeini, a->stiini, a->ener, a->enerini};
    ih[0] = *a->iout; ih[1] = *a->nmethod; ih[2] = a->iperturb[0]; ih[3] = a->iperturb[1]; ih[4] = *a->ielas; ih[5] = *a->icmd;
    ih[6] = *a->istep; ih[7] = *a->iinc; ih[8] = a->calcul_fn; ih[9] = a->calcul_qa; ih[10] = a->calcul_cauchy; ih[11] = *a->nener;
    ih[12] = a->ikin; ih[13] = *a->kscale; ih[14] = 0; ih[15] = 0;
    for (k = 0; k < DR_N; k++) { glo[k] = g[k]; if (g[k] != NULL) ih[14] |= (1 << k); }
    if (*a->nener != 1) ih[14] &= ~((1 << DR_ENER) | (1 << DR_ENERINI));
    dh[0] = *a->dtime; dh[1] = *a->time; dh[2] = *a->ttime; dh[3] = *a->reltime; dh[4] = *a->penal; dh[5] = *a->rhomin;
  }
  MPI_Bcast(ih, 16, MPI_ITG, 0, PETSC_COMM_WORLD);
  MPI_Bcast(dh, 6, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  have = (int)ih[14];
  d_scatter(me == 0 ? a->v : NULL, bd_tnode, dm.v, dm.mt * dm.nk, TAG_D + 3);
  d_scatter(me == 0 ? a->vold : NULL, bd_tnode, dm.vold, dm.mt * dm.nk, TAG_D + 4);
  for (k = 0; k < DR_N; k++) {
    MPI_Datatype *t = (w[k] == 27) ? bd_te27 : (w[k] == 6) ? bd_te6 : bd_te1;
    size_t        nl = (size_t)w[k] * dm.mi[0] * dm.ne;
    if (have & (1 << k)) d_scatter(glo[k], t, loc[k], (int)nl, TAG_D + 10 + k);
    else memset(loc[k], 0, sizeof(double) * nl);
  }
  d_res_local(ih, dh, qa, &nal);
  /* element arrays back */
  for (k = 0; k < DR_N; k++) {
    MPI_Datatype *t = (w[k] == 27) ? bd_te27 : (w[k] == 6) ? bd_te6 : bd_te1;
    if (out[k] && (have & (1 << k))) d_gather(glo[k], t, loc[k], w[k] * dm.mi[0] * (int)dm.ne, TAG_D + 20 + k);
  }
  /* internal forces: added at the global nodes, rank by rank (fixed order) */
  if (me == 0) {
    ITG    mt = dm.mt;
    double *tmp = NULL;
    int    mx = 0, c;
    for (r = 1; r < bd_size; r++) if (bd_ncnt[r] > mx) mx = bd_ncnt[r];
    tmp = malloc(sizeof(double) * ((size_t)mt * mx + 1));
    for (r = 0; r < bd_size; r++) {
      double *src = dm.fn;
      if (r > 0) { MPI_Recv(tmp, (int)mt * bd_ncnt[r], MPI_DOUBLE, r, TAG_D + 5, PETSC_COMM_WORLD, MPI_STATUS_IGNORE); src = tmp; }
      for (k = 0; k < bd_ncnt[r]; k++)
        for (c = 0; c < mt; c++) a->fn[(size_t)mt * bd_nlst[r][k] + c] += src[(size_t)mt * k + c];
    }
    free(tmp);
  } else MPI_Send(dm.fn, (int)(dm.mt * dm.nk), MPI_DOUBLE, 0, TAG_D + 5, PETSC_COMM_WORLD);
  /* qa and nal: the rules of results.c for its threads, over the ranks */
  {
    double mine[5] = {qa[0], qa[1], qa[2], qa[3], (double)nal};
    if (me == 0) gq = malloc(sizeof(double) * 5 * bd_size);
    MPI_Gather(mine, 5, MPI_DOUBLE, gq, 5, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    if (me == 0) {
      double q0 = gq[0], q1 = gq[1], q2 = gq[2], q3 = gq[3], na = gq[4];
      for (r = 1; r < bd_size; r++) {
        double *g = gq + 5 * r;
        q0 += g[0]; q1 += g[1]; na += g[4];
        if (g[2] > 0.) { if (q2 < 0.) q2 = g[2]; else if (g[2] < q2) q2 = g[2]; }
        if (g[3] > 0.) { if (q3 < 0.) q3 = g[3]; else if (g[3] > q3) q3 = g[3]; }
      }
      a->qa[0] = q0; a->qa[1] = q1; a->qa[2] = q2; a->qa[3] = q3; *a->nal = (ITG)na;
      free(gq);
    }
  }
}

/* rank 0 entry points (called from mafillsmmain and results) */
int bddc_dist_on(void) { return bddc_active() && bd_dist == 1; }

int bddc_dist_asm(bddc_eargs *a)
{
  int    nt, nm;
  double t0;
  if (!bddc_active() || !bd_topo || bd_dist < 0 || *a->stiffness != 1) return 0;
  if (bd_dist == 0) {
    int ok;
    TAU_REGION_START("Domain decomposition");   /* local models sent to the ranks (once) */
    ok = d_setup(a);
    TAU_REGION_STOP("Domain decomposition");
    if (!ok) return 0;
  }
  t0 = MPI_Wtime();
  nt = omp_get_max_threads(); nm = mkl_get_max_threads();
  pin_own(); omp_set_num_threads(1);
  send_cmd(CMD_DASM, 0, 0, 0);
  d_asm_collective(a);
  unpin_all(nt, nm);
  bd_tdasm += MPI_Wtime() - t0; bd_ndasm++;
  return 1;
}

int bddc_dist_res(bddc_eargs *a)
{
  int    nt, nm;
  double t0;
  if (!bddc_active() || bd_dist != 1) return 0;
  t0 = MPI_Wtime();
  nt = omp_get_max_threads(); nm = mkl_get_max_threads();
  pin_own(); omp_set_num_threads(1);
  send_cmd(CMD_DRES, 0, 0, 0);
  d_res_collective(a);
  unpin_all(nt, nm);
  bd_tdres += MPI_Wtime() - t0; bd_ndres++;
  return 1;
}

/* ---------------- worker loop ---------------- */

static void worker_loop(void)
{
  double          msg[4];
  int             flag, mynel = 0, *myedof = NULL;
  MPI_Status      st;
  struct timespec ts = {0, 1000000};       /* 1 ms: sleep, do not spin */
  for (;;) {
    MPI_Iprobe(0, TAG_CMD, PETSC_COMM_WORLD, &flag, &st);
    if (!flag) { nanosleep(&ts, NULL); continue; }
    MPI_Recv(msg, 4, MPI_DOUBLE, 0, TAG_CMD, PETSC_COMM_WORLD, &st);
    if ((int)msg[0] == CMD_QUIT) return;
    if ((int)msg[0] == CMD_TOPO) {
      PetscInt n = (PetscInt)msg[1];
      double  *coeq = malloc(sizeof(double) * 4 * (size_t)n);
      bd_s = (int)msg[3]; bd_s2 = bd_s * bd_s; bd_nn = bd_s / 3;
      MPI_Recv(&mynel, 1, MPI_INT, 0, TAG_DATA, PETSC_COMM_WORLD, &st);
      myedof = malloc(sizeof(int) * bd_s * (size_t)(mynel + 1));
      MPI_Recv(myedof, bd_s * mynel, MPI_INT, 0, TAG_DATA, PETSC_COMM_WORLD, &st);
      MPI_Bcast(coeq, 4 * (int)n, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
      if (topo_collective(myedof, mynel, n, coeq, (PetscInt)msg[2])) MPI_Abort(PETSC_COMM_WORLD, 1);
      free(coeq);
    } else if ((int)msg[0] == CMD_MATRIX) {
      MPI_Recv(lK, bd_s2 * nown, MPI_DOUBLE, 0, TAG_DATA, PETSC_COMM_WORLD, &st);
      if (matrix_collective(NULL)) MPI_Abort(PETSC_COMM_WORLD, 1);
    } else if ((int)msg[0] == CMD_PMAT) {
      MPI_Recv(lK, bd_s2 * nown, MPI_DOUBLE, 0, TAG_DATA, PETSC_COMM_WORLD, &st);
      if (pmat_collective(NULL)) MPI_Abort(PETSC_COMM_WORLD, 1);
    } else if ((int)msg[0] == CMD_LMATRIX) {
      if (matrix_collective(NULL)) MPI_Abort(PETSC_COMM_WORLD, 1);
    } else if ((int)msg[0] == CMD_LPMAT) {
      if (pmat_collective(NULL)) MPI_Abort(PETSC_COMM_WORLD, 1);
    } else if ((int)msg[0] == CMD_DMODEL) {
      TAU_REGION_START("Domain decomposition");
      d_recv_local();
      d_bcast_materials();
      TAU_REGION_STOP("Domain decomposition");
    } else if ((int)msg[0] == CMD_DASM) {
      d_asm_collective(NULL);
    } else if ((int)msg[0] == CMD_DRES) {
      d_res_collective(NULL);
    } else if ((int)msg[0] == CMD_SOLVE) {
      PetscInt     its;
      double       res;
      int          reason;
      PetscScalar *a;
      VecGetArray(vb, &a);
      MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, a, (int)m_own, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
      VecRestoreArray(vb, &a);
      if (solve_collective(msg[1], (int)msg[2], &its, &res, &reason)) MPI_Abort(PETSC_COMM_WORLD, 1);
      VecGetArray(vx, &a);
      MPI_Gatherv(a, (int)m_own, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
      VecRestoreArray(vx, &a);
    }
  }
}

/* ---------------- rank 0 ---------------- */


/* after an assembly: first time, build the topology and the partition */
void bddc_capture_end(ITG *kon, ITG *ipkon, ITG ne, ITG *nactdof, ITG *mi, ITG nk, double *co, ITG neq)
{
  ITG     mt = mi[1] + 1, e, j, c, k;
  double  t0 = omp_get_wtime(), *coeq;
  int     r, *cnt, nact = 0;
  if (!bddc_active() || !bd_capturing) return;
  bd_capturing = 0;
  bd_dirty = 1;
  if (bd_topo) return;
  TAU_REGION_START("Domain decomposition");
  /* every element in use must have handed over its matrix */
  bd_nel = 0;
  bd_el = malloc(sizeof(int) * (size_t)ne);
  for (e = 0; e < ne; e++) {
    if (ipkon[e] < 0) continue;
    nact++;
    if (bd_have[e]) bd_el[bd_nel++] = (int)e;
  }
  if (bd_nel == 0 || bd_nel != nact) {
    char msg[160];
    snprintf(msg, sizeof(msg), "only %d of %d elements handed over a stiffness matrix (C3D4/C3D10 with 3 dofs per node)", bd_nel, nact);
    bd_fatal(msg);
  }
  /* PETSc numbering: if every node has all 3 directions free or all 3
     fixed (no symmetry/roller supports), the unknowns are grouped per node
     (block size 3: faster BDDC setup and solves); otherwise CalTop's
     equation order, one unknown at a time */
  {
    ITG nb = 0, partial = 0;
    bd_g = malloc(sizeof(int) * (size_t)neq);
    for (k = 0; k < nk; k++) {
      int na = 0;
      for (c = 0; c < 3; c++) if (nactdof[mt * k + c + 1] > 0) na++;
      if (na > 0 && na < 3) partial++;
    }
    bs_g = (partial == 0) ? 3 : 1;
    if (bs_g == 3) {
      /* nodes numbered in the order they first appear in the element list
         (fixed nodes left out); the METIS partition below uses the same
         numbering (fast setup: compact subdomain numbering, block size 3) */
      bd_nodeid = malloc(sizeof(int) * (size_t)nk);
      for (k = 0; k < nk; k++) bd_nodeid[k] = -1;
      for (e = 0; e < bd_nel; e++)
        for (j = 0; j < bd_enope[bd_el[e]]; j++) {
          k = kon[ipkon[bd_el[e]] + j] - 1;
          if (bd_nodeid[k] < 0 && nactdof[mt * k + 1] > 0) bd_nodeid[k] = (int)nb++;
        }
      for (k = 0; k < nk; k++)
        if (bd_nodeid[k] >= 0)
          for (c = 0; c < 3; c++) {
            ITG q = nactdof[mt * k + c + 1];
            if (q > 0 && q <= neq) bd_g[q - 1] = 3 * bd_nodeid[k] + (int)c;
          }
    }
    if (bs_g == 3 && 3 * nb != neq) { bs_g = 1; free(bd_nodeid); bd_nodeid = NULL; } /* safety: identity */
    if (bs_g == 1)
      for (k = 0; k < neq; k++) bd_g[k] = (int)k;
    printf(" BDDC: unknowns %s (%d nodes with partly fixed directions)\n",
           bs_g == 3 ? "grouped per node (x,y,z), nodes numbered in element order" : "one at a time", (int)partial);
  }
  bd_edof = malloc(sizeof(int) * bd_s * (size_t)bd_nel);
  for (k = 0; k < bd_nel; k++) {
    e = bd_el[k];
    for (j = 0; j < bd_nn; j++)
      for (c = 0; c < 3; c++) {
        ITG q = (j < bd_enope[e]) ? nactdof[mt * (kon[ipkon[e] + j] - 1) + c + 1] : 0;   /* 0: padding */
        bd_edof[bd_s * k + 3 * j + c] = (q > 0) ? bd_g[q - 1] : -1;
      }
  }
  /* coordinates and component of every unknown (PETSc order) */
  coeq = calloc(4 * (size_t)neq, sizeof(double));
  for (k = 0; k < nk; k++)
    for (c = 0; c < 3; c++) {
      ITG q = nactdof[mt * k + c + 1];
      if (q > 0 && q <= neq) {
        int g = bd_g[q - 1];
        for (j = 0; j < 3; j++) coeq[4 * (size_t)g + j] = co[3 * k + j];
        coeq[4 * (size_t)g + 3] = (double)c;
      }
    }
  /* partition of the elements */
  bd_part = calloc(bd_nel, sizeof(int));
  if (bd_size > 1) {
#if defined(PETSC_HAVE_METIS)
    idx_t ne_ = bd_nel, nn_ = bd_nodeid ? (idx_t)(neq / 3) : (idx_t)nk, ncommon = 3, np_ = bd_size, obj, *eptr, *eind, *ep, *npt, opt[METIS_NOPTIONS], q = 0;
    eptr = malloc(sizeof(idx_t) * (bd_nel + 1)); eind = malloc(sizeof(idx_t) * bd_nn * (size_t)bd_nel);
    ep = malloc(sizeof(idx_t) * bd_nel); npt = malloc(sizeof(idx_t) * (nn_ + 1));
    eptr[0] = 0;
    for (k = 0; k < bd_nel; k++) {
      for (j = 0; j < bd_enope[bd_el[k]]; j++) {
        ITG nd = kon[ipkon[bd_el[k]] + j] - 1;
        if (!bd_nodeid) eind[q++] = nd;
        else if (bd_nodeid[nd] >= 0) eind[q++] = bd_nodeid[nd];   /* fixed nodes left out */
      }
      eptr[k + 1] = q;
    }
    METIS_SetDefaultOptions(opt);
    opt[METIS_OPTION_NUMBERING] = 0;
    /* contiguous subdomains: a subdomain in several pieces has a floating
       piece that the BDDC vertices do not hold, and its Schur complement is
       singular (64 ranks on Delta: zero pivot in PCBDDCSubSchursSetUp) */
    opt[METIS_OPTION_CONTIG] = 1;
    if (METIS_PartMeshDual(&ne_, &nn_, eptr, eind, NULL, NULL, &ncommon, &np_, NULL, opt, &obj, ep, npt) != METIS_OK) {
      printf(" BDDC: METIS cannot make contiguous subdomains (mesh graph not connected?), partitioning without\n");
      opt[METIS_OPTION_CONTIG] = 0;
      if (METIS_PartMeshDual(&ne_, &nn_, eptr, eind, NULL, NULL, &ncommon, &np_, NULL, opt, &obj, ep, npt) != METIS_OK) {
        bd_fatal("METIS could not partition the mesh");
      }
    }
    for (k = 0; k < bd_nel; k++) bd_part[k] = (int)ep[k];
    free(eptr); free(eind); free(ep); free(npt);
#else
    /* slabs along the longest axis of the element centroids */
    {
      double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
      int    ax = 0;
      slab_t *sl = malloc(sizeof(slab_t) * bd_nel);
      for (k = 0; k < bd_nel; k++)
        for (c = 0; c < 3; c++) {
          double v = co[3 * (kon[ipkon[bd_el[k]]] - 1) + c];
          if (v < lo[c]) lo[c] = v;
          if (v > hi[c]) hi[c] = v;
        }
      for (c = 1; c < 3; c++) if (hi[c] - lo[c] > hi[ax] - lo[ax]) ax = c;
      for (k = 0; k < bd_nel; k++) {
        sl[k].v = 0.; sl[k].i = (int)k;
        for (j = 0; j < 4; j++) sl[k].v += co[3 * (kon[ipkon[bd_el[k]] + j] - 1) + ax];
      }
      qsort(sl, bd_nel, sizeof(slab_t), cmpslab);
      for (k = 0; k < bd_nel; k++) bd_part[sl[k].i] = (int)((long long)k * bd_size / bd_nel);
      free(sl);
    }
#endif
  }
  {                                     /* pieces per subdomain (elements sharing a face = connected) */
    int   *up = malloc(sizeof(int) * (size_t)bd_nel), *pieces = calloc(bd_size, sizeof(int)), worst = 0, nbad = 0;
    size_t nf = 4 * (size_t)bd_nel, f;
    long long *key = malloc(sizeof(long long) * nf * 2);
    static const int fc[4][3] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    for (k = 0; k < bd_nel; k++) up[k] = (int)k;
    for (k = 0; k < bd_nel; k++)          /* faces as sorted corner triples, with their element */
      for (j = 0; j < 4; j++) {
        long long a[3], t;
        for (c = 0; c < 3; c++) a[c] = kon[ipkon[bd_el[k]] + fc[j][c]];
        if (a[0] > a[1]) { t = a[0]; a[0] = a[1]; a[1] = t; }
        if (a[1] > a[2]) { t = a[1]; a[1] = a[2]; a[2] = t; }
        if (a[0] > a[1]) { t = a[0]; a[0] = a[1]; a[1] = t; }
        key[2 * (4 * (size_t)k + j)] = (a[0] * 2000003LL + a[1]) * 2000003LL + a[2];
        key[2 * (4 * (size_t)k + j) + 1] = k;
      }
    qsort(key, nf, 2 * sizeof(long long), cmp_ll2);
    for (f = 1; f < nf; f++)
      if (key[2 * f] == key[2 * f - 2]) {
        int x = (int)key[2 * f + 1], y = (int)key[2 * f - 1];
        if (bd_part[x] != bd_part[y]) continue;
        while (up[x] != x) x = up[x] = up[up[x]];
        while (up[y] != y) y = up[y] = up[up[y]];
        if (x != y) up[x] = y;
      }
    for (k = 0; k < bd_nel; k++) if (up[k] == (int)k) pieces[bd_part[k]]++;
    for (r = 0; r < bd_size; r++) { if (pieces[r] > 1) nbad++; if (pieces[r] > worst) worst = pieces[r]; }
    printf(" BDDC: subdomains in more than one piece: %d of %d (most pieces in one subdomain: %d)\n", nbad, bd_size, worst);
    free(up); free(pieces); free(key);
  }
  /* per-rank element lists, datatypes for the element matrices */
  cnt = calloc(bd_size, sizeof(int));
  for (k = 0; k < bd_nel; k++) cnt[bd_part[k]]++;
  bd_type = malloc(sizeof(MPI_Datatype) * bd_size);
  {
    MPI_Datatype blk;
    int        **lst = malloc(sizeof(int *) * bd_size), *fill = calloc(bd_size, sizeof(int));
    MPI_Type_contiguous(bd_s2, MPI_DOUBLE, &blk);
    for (r = 0; r < bd_size; r++) lst[r] = malloc(sizeof(int) * (cnt[r] + 1));
    for (k = 0; k < bd_nel; k++) { r = bd_part[k]; lst[r][fill[r]++] = bd_el[k]; }
    for (r = 0; r < bd_size; r++) {
      MPI_Type_create_indexed_block(cnt[r], 1, lst[r], blk, &bd_type[r]);
      MPI_Type_commit(&bd_type[r]);
    }
    /* topology to the workers */
    send_cmd(CMD_TOPO, (double)neq, (double)bs_g, (double)bd_s);
    for (r = 1; r < bd_size; r++) {
      int *ed = malloc(sizeof(int) * bd_s * (size_t)(cnt[r] + 1)), q = 0;
      for (k = 0; k < bd_nel; k++)
        if (bd_part[k] == r) { memcpy(ed + (size_t)bd_s * q, bd_edof + (size_t)bd_s * k, sizeof(int) * bd_s); q++; }
      MPI_Send(&cnt[r], 1, MPI_INT, r, TAG_DATA, PETSC_COMM_WORLD);
      MPI_Send(ed, bd_s * cnt[r], MPI_INT, r, TAG_DATA, PETSC_COMM_WORLD);
      free(ed);
    }
    MPI_Bcast(coeq, 4 * (int)neq, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    {
      int *ed0 = malloc(sizeof(int) * bd_s * (size_t)(cnt[0] + 1)), q = 0;
      for (k = 0; k < bd_nel; k++)
        if (bd_part[k] == 0) { memcpy(ed0 + (size_t)bd_s * q, bd_edof + (size_t)bd_s * k, sizeof(int) * bd_s); q++; }
      if (topo_collective(ed0, cnt[0], (PetscInt)neq, coeq, bs_g)) MPI_Abort(PETSC_COMM_WORLD, 1);
      /* rank 0's own element list, for copying its matrices */
      bd_n0 = cnt[0];
      bd_own0 = malloc(sizeof(int) * (cnt[0] + 1));
      memcpy(bd_own0, lst[0], sizeof(int) * cnt[0]);   /* ed0 is kept: local numbering */
    }
    bd_elst = lst;                     /* kept: per-rank element lists (distributed element work) */
    bd_ecnt = malloc(sizeof(int) * bd_size);
    memcpy(bd_ecnt, cnt, sizeof(int) * bd_size);
    free(fill);
    MPI_Type_free(&blk);
  }
  free(coeq);
  /* vector layout (owned rows per rank) */
  bd_cnt = malloc(sizeof(int) * bd_size);
  bd_dsp = malloc(sizeof(int) * bd_size);
  for (r = 0; r < bd_size; r++) {
    PetscInt nb = (PetscInt)neq / bs_g, a0 = bs_g * (nb * r / bd_size), a1 = bs_g * (nb * (r + 1) / bd_size);
    bd_cnt[r] = (int)(a1 - a0); bd_dsp[r] = (int)a0;
  }
  bd_vp = malloc(sizeof(double) * (size_t)neq);
  bd_topo = 1;
  printf(" BDDC: %d elements in %d subdomains (%s), %d equations, %.1f s\n", bd_nel, bd_size,
#if defined(PETSC_HAVE_METIS)
         "METIS",
#else
         "slabs",
#endif
         (int)neq, omp_get_wtime() - t0);
  free(cnt);
  TAU_REGION_STOP("Domain decomposition");
}

static int ship_matrix(void)
{
  double t0 = omp_get_wtime(), ts = 0.;
  int    r, k;
  if (bd_dist == 1) send_cmd(CMD_LMATRIX, 0, 0, 0);   /* element matrices already on their ranks */
  else {
    send_cmd(CMD_MATRIX, 0, 0, 0);
    for (r = 1; r < bd_size; r++) MPI_Send(bd_K, 1, bd_type[r], r, TAG_DATA, PETSC_COMM_WORLD);
    for (k = 0; k < bd_n0; k++) memcpy(lK + (size_t)bd_s2 * k, bd_K + (size_t)bd_s2 * bd_own0[k], sizeof(double) * bd_s2);
  }
  bd_tship += omp_get_wtime() - t0;
  if (matrix_collective(&ts)) MPI_Abort(PETSC_COMM_WORLD, 1);
  bd_nmat++;
  if (ts > 0.) {
    bd_tsetup = ts;
    if (bd_rebuild_its < 0) bd_rebuild_its = (int)getd("CALTOP_BDDC_REBUILD", (double)BDDC_REBUILD_DEFAULT);
    printf(" BDDC: preconditioner set up from the stress-free matrix in %.1f s\n", ts);
    if (bd_rebuild_its > 0)
      printf(" BDDC: rebuilt at the next Newton iteration after a solve needs more than %d CG iterations (CALTOP_BDDC_REBUILD)\n", bd_rebuild_its);
    else
      printf(" BDDC: not rebuilt during the Newton iterations (CALTOP_BDDC_REBUILD=0)\n");
    printf(" BDDC: rebuilt once at the converged configuration if adjoint solves follow\n");
  }
  bd_dirty = 0;
  return 1;
}

int bddc_want_rebuild(void)
{
  return bddc_active() && pc_ready && bd_want && (bd_K != NULL || bd_dist == 1);
}

/* the matrix just captured (zero stresses, current displacements: the
   material stiffness in the current configuration, positive semi-definite
   element by element) becomes the preconditioner matrix */
void bddc_rebuild(void)
{
  double t0 = omp_get_wtime(), ts = 0.;
  int    r, k, nt = omp_get_max_threads(), nm = mkl_get_max_threads();
  if (!bddc_active() || !pc_ready || (bd_K == NULL && bd_dist != 1)) return;
  bd_want = 0;
  pin_own();
  omp_set_num_threads(1);
  if (bd_dist == 1) send_cmd(CMD_LPMAT, 0, 0, 0);    /* zero-stress element matrices already on their ranks */
  else {
    send_cmd(CMD_PMAT, 0, 0, 0);
    for (r = 1; r < bd_size; r++) MPI_Send(bd_K, 1, bd_type[r], r, TAG_DATA, PETSC_COMM_WORLD);
    for (k = 0; k < bd_n0; k++) memcpy(lK + (size_t)bd_s2 * k, bd_K + (size_t)bd_s2 * bd_own0[k], sizeof(double) * bd_s2);
  }
  bd_tship += omp_get_wtime() - t0;
  if (pmat_collective(&ts)) MPI_Abort(PETSC_COMM_WORLD, 1);
  unpin_all(nt, nm);
  bd_nrebuild++;
  bd_trebuild += ts;
  bd_dirty = 1;                       /* the Newton matrix (K_T) follows */
  printf(" BDDC: preconditioner rebuilt from the stress-free matrix at the current configuration in %.1f s (rebuild %d)\n", ts, bd_nrebuild);
  fflush(stdout);
}

/* solve K x = b (b overwritten by x); a solve that does not converge stops the run */
int bddc_solve(double *b, ITG neq, double tol, int adjoint, double *ad, double *au, ITG *jq, ITG *irow)
{
  double       t0, res;
  ITG          i;
  PetscInt     its;
  int          reason, maxit = 500, nt = omp_get_max_threads(), nm = mkl_get_max_threads(), ok;
  double       accept = 1.e-5;
  PetscScalar *a;
  if (!bddc_active() || bd_topo == 0 || (PetscInt)neq != n_g) {
    printf(" *ERROR: BDDC: solve requested but BDDC is not set up for this system (%d equations, %d expected)\n", (int)neq, (int)n_g);
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  pin_own();          /* rank 0 back on its own core: no sharing with a solver rank */
  {                   /* and its threads as a solver rank's: one per CPU of that core */
    cpu_set_t cs;
    int       nc = 1;
    if (sched_getaffinity(0, sizeof(cs), &cs) == 0 && CPU_COUNT(&cs) > 0) nc = CPU_COUNT(&cs);
    if (nc > 4) nc = 1;                     /* not pinned to one core: stay single-threaded */
    omp_set_num_threads(nc);
    mkl_set_num_threads(1);
#pragma omp parallel num_threads(nc)
    sched_setaffinity(0, sizeof(bd_core), &bd_core);   /* CalTop's pool threads onto this core too */
  }
  if (bd_dirty) {
    if (!ship_matrix()) {
      printf(" *ERROR: BDDC: the element matrices could not be handed to the solver ranks; stopping\n");
      fflush(stdout);
      MPI_Abort(PETSC_COMM_WORLD, 1);
    }
  }
  t0 = MPI_Wtime();
  send_cmd(CMD_SOLVE, tol, (double)maxit, 0.);
  for (i = 0; i < neq; i++) bd_vp[bd_g[i]] = b[i];
  VecGetArray(vb, &a);
  MPI_Scatterv(bd_vp, bd_cnt, bd_dsp, MPI_DOUBLE, a, (int)m_own, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
  VecRestoreArray(vb, &a);
  if (solve_collective(tol, maxit, &its, &res, &reason)) MPI_Abort(PETSC_COMM_WORLD, 1);
  unpin_all(nt, nm);
  ok = (reason > 0) || (res <= (tol > accept ? tol : accept));
  {
    double *x = malloc(sizeof(double) * neq);   /* gathered either way (collective) */
    VecGetArray(vx, &a);
    MPI_Gatherv(a, (int)m_own, MPI_DOUBLE, bd_vp, bd_cnt, bd_dsp, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    VecRestoreArray(vx, &a);
    if (ok) for (i = 0; i < neq; i++) b[i] = bd_vp[bd_g[i]];
    free(x);
  }
  t0 = MPI_Wtime() - t0;
  bd_tsolve += t0;
  if (ok) { bd_nsolve++; bd_nit += its; } else bd_nfail++;
  if (adjoint && bd_K != NULL) {
    /* the converged K_T is on the solver ranks: the element matrices on
       rank 0 (900 doubles per element) are not needed for the adjoints */
    free(bd_K); free(bd_have);
    bd_K = NULL; bd_have = NULL; bd_ne = 0;
  }
  printf(" BDDC: %s CG %d it, |b-Kx|/|b| = %.1e (asked %.0e), %.2f s\n", adjoint ? "adjoint" : "Newton", (int)its, res, tol, t0);
  if (!ok) {
    printf(" *ERROR: BDDC: %s CG did not converge (true residual %.1e, PETSc reason %d); stopping.\n"
           "         Try CALTOP_BDDC_REBUILD lower, or more CG iterations via CALTOP_BDDC_OPTS\n", adjoint ? "adjoint" : "Newton", res, reason);
    fflush(stdout);
    MPI_Abort(PETSC_COMM_WORLD, 1);
  }
  if (bd_rebuild_its < 0) bd_rebuild_its = (int)getd("CALTOP_BDDC_REBUILD", (double)BDDC_REBUILD_DEFAULT);
  if (!adjoint && bd_rebuild_its > 0 && its > bd_rebuild_its && !bd_want) {
    bd_want = 1;
    printf(" BDDC: more than %d CG iterations (CALTOP_BDDC_REBUILD): the preconditioner is rebuilt at the next Newton iteration\n", bd_rebuild_its);
  }
  return ok;
}

static void bddc_finish(void)
{
  if (!bd_on) return;
  if (bd_nmat > 0)
    printf(" BDDC: matrices=%ld solves=%ld (avg %.1f it) failed=%ld rebuilds=%d; time: setup %.1f s, rebuilds %.1f s, matrix transfer %.1f s, solves %.1f s\n", bd_nmat,
           bd_nsolve, bd_nsolve ? (double)bd_nit / bd_nsolve : 0., bd_nfail, bd_nrebuild, bd_tsetup, bd_trebuild, bd_tship, bd_tsolve);
  if (bd_dist == 1)
    printf(" BDDC: distributed element work: %d assemblies %.1f s, %d internal-force calls %.1f s\n", bd_ndasm, bd_tdasm, bd_ndres, bd_tdres);
  if (!bd_disabled) send_cmd(CMD_QUIT, 0, 0, 0);
  bd_on = 0;
  fflush(stdout);
  PetscFinalize();
}

#else /* no PETSc: stubs */

void bddc_launch(int *argc, char ***argv)
{
  char *s = getenv("CALTOP_SOLVER");
  (void)argc; (void)argv;
  if (s != NULL && strcmp(s, "bddc") == 0) {
    printf(" *ERROR: CALTOP_SOLVER=bddc, but this calTop.exe was built without PETSc (make BDDC=1)\n");
    fflush(stdout);
    exit(1);
  }
}
int  bddc_active(void) { return 0; }
int  bddc_want_rebuild(void) { return 0; }
ITG  ctnokmat_(void) { return 0; }
void bddc_rebuild(void) {}
int  bddc_dist_asm(bddc_eargs *a) { (void)a; return 0; }
int  bddc_dist_on(void) { return 0; }
int  bddc_dist_res(bddc_eargs *a) { (void)a; return 0; }
void bddc_capture_begin(ITG ne, ITG *ipkon, char *lakon) { (void)ne; (void)ipkon; (void)lakon; }
void bddc_elem(ITG e, ITG nope, ITG ndof, double *sm) { (void)e; (void)nope; (void)ndof; (void)sm; }
void ctelemcapture_(ITG *iel, ITG *nope, ITG *ndof, double *sm) { (void)iel; (void)nope; (void)ndof; (void)sm; }
void bddc_capture_end(ITG *kon, ITG *ipkon, ITG ne, ITG *nactdof, ITG *mi, ITG nk, double *co, ITG neq)
{
  (void)kon; (void)ipkon; (void)ne; (void)nactdof; (void)mi; (void)nk; (void)co; (void)neq;
}
int bddc_solve(double *b, ITG neq, double tol, int adjoint, double *ad, double *au, ITG *jq, ITG *irow)
{
  (void)b; (void)neq; (void)tol; (void)adjoint; (void)ad; (void)au; (void)jq; (void)irow;
  return 0;
}
#endif
