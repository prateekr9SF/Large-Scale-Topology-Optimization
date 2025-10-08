/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

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

int main(int argc,char *argv[])
{

  FILE *f1;

  char *sideload=NULL;    /**< load label for distributed loads. */
  char *set=NULL;         /**< node set or element set */
  char *matname=NULL;     /**< material name */
  char *orname=NULL;      /**< name of orientation of element */
  char *amname=NULL;      /**< amplitude name. */
  char *filab=NULL;       /**< field label (U, S,..) */
  char *lakon=NULL;       /**< element type string */
  char *labmpc=NULL;      /**< MPC label. */
  char *prlab=NULL;       /**< output request */
  char *prset=NULL;       /**< node or element set assiciated with output request.*/
  char  jobnamec[660]=""; /**< jobname */
  char  jobnamef[132]=""; /**< jobname for fluid analysis. */
  char  output[4]="asc";  /**< output type */
  char  *typeboun=NULL;   /**< boundary condition type */
  char  *inpc=NULL;       /**< processed inp file content */
  char  *tieset=NULL;     /**< tie constraints */
  char  *cbody=NULL;      /**< element set to which body load applies */
  char  fneig[132]="";    /**< file name for eigenvalue calculations */
  char  *sideloadtemp=NULL; /**< temperature value associated with a sideload. */
  char  kind1[2]="T";     /**< decription. */
  char  kind2[2]="T";     /**< decription. */
  char  *heading=NULL;    /**< decription. */
  char  *objectset=NULL;  /**< decription. */

  ITG *kon=NULL;          /**< element connectivity */
  ITG *nodeboun=NULL;     /**< boundary nodes (SPC) */
  ITG *ndirboun=NULL;     /**< directions of SPC at nodeboun */
  ITG *ipompc=NULL;       /**< index of node with MPCs */ 
  ITG *nodempc=NULL;      /**< node numbers with MPCs and DOFs */
  ITG *nodeforc=NULL;     /**< node numbers where forces are applied */
  ITG *ndirforc=NULL;     /**< direction of applied forces */
  ITG  *nelemload=NULL;   /**< elements associated with facial distributed loads */
  ITG im;                 /**< general purpose loop variable */
  ITG *inodesd=NULL;      /**< node numbers associated with active D.O.Fs */
  ITG nload1;             /**< number of facial distributed loads applied to elements */
  ITG *idefforc=NULL;     /**< indicates whether a load is effectively applied to  specific node and direction */
  ITG  *nactdof=NULL;     /**< active D.O.Fs for nodes in system */
  ITG *icol=NULL;         /**<stores the column numbers for the sparse matrix storage of D.O.Fs */
  ITG *ics=NULL;          /**< array for identifying whether a certain node belongs to a constraint set */
  ITG  *jq=NULL;          /**< stores the row numbers for the spatse matrix structure of global stiffness matrix */
  ITG *mast1=NULL;        /**< master D.O.Fs for the MPC system */
  ITG *irow=NULL;         /**< row indices of non-zero entries in a sparse matrix represenation */
  ITG *rig=NULL;          /**< rigid body modes */
  ITG *idefbody=NULL;     /**< body force definition */
  ITG  *ikmpc=NULL;       /**< global D.O.Fs of the dependent terms in MPCs */
  ITG *ilmpc=NULL;        /**< MPICs indices in relation to their dependent terms */
  ITG *ikboun=NULL;       /**< sorted D.O.Fs for SPCs */
  ITG *ilboun=NULL;       /**< B.C node indices for SPCs */
  ITG  *nreorder=NULL;    /**< Re-order node/element array for refinement */
  ITG *ipointer=NULL;     /**< pointer for sparse matrix storage */
  ITG *idefload=NULL;     /**< Load defined on same element/load label previously during analysis */
  ITG  *istartset=NULL;   /**< stores starting position of element or node sets in certain operations */
  ITG *iendset=NULL;      /**< stores starting position of element or node sets in certain operations */
  ITG *ialset=NULL;       /**< information about the lements or nodes that belong to a set */
  ITG *ielmat=NULL;       /**< material properties associated with each element */
  ITG *ielorien=NULL;     /**< oreientation data for each element */
  ITG *nrhcon=NULL;       /**< element heat conductivity values */
  ITG *nodebounold=NULL;  /**< node numbers of SPCs from previous step/increment */
  ITG *ndirbounold=NULL;  /**< node direction of SPCs from previous step/increment */
  ITG *nelcon=NULL;        /**< elasticity constants for each element */
  ITG *nalcon=NULL;       /**< load and boundary condition amplitude definition */
  ITG *iamforc=NULL;      /**< concentrated force definition */
  ITG *iamload=NULL;     /**< amplitude definition for distributed loads */
  ITG *iamt1=NULL;      /**< amplitude definition for temperature loads */
  ITG *namta=NULL;     /**< amplitude values for temperature definitions */
  ITG *ipkon=NULL;    /** < element pointer to kon */
  ITG *iamboun=NULL; /**< ampplitude definitins associated with SPCs */
  ITG *nplicon=NULL; /**< number of plasticity constants for material */
  ITG *nplkcon=NULL; /**< number of plasticity constants associated with elements */
  ITG *inotr=NULL;   /**< nodal transformation data */
  ITG *iponor=NULL; /**< node-element connectivity */
  ITG *knor=NULL;  /**< face-normal vectors */
  ITG *ikforc=NULL; /**< ordered D.O.Fs for concentrated loads */
  ITG *ilforc=NULL; /**< load indices for concentrated loads */
  ITG *iponoel=NULL; /**< node-element connectivity for each node */
  ITG *inoel=NULL; /**< node numbers associated with each element */
  ITG *nshcon=NULL; /**< number of shell conductivity values per element */
  ITG *ncocon=NULL; /**< number of convection-related constants for the material */
  ITG *ibody=NULL; /**< body-forces-relatred information for the elements */
  ITG *ielprop=NULL; /**< element properties */
  ITG *islavsurf=NULL; /**< surface-slave information for contact surfaces */
  ITG *ipoinpc=NULL; /**< MPC location in input deck */
  ITG mt; /**< number of D.O.Fs in the system */
  ITG nxstate; /**< number of state variabels */
  ITG nload0; /**< number of load aplications at the start of analysis */
  ITG iload; /**< loop counter for handling loads in system */
  ITG *iuel=NULL; /**< element numbers involved in user-defined subroutines */

  ITG nk; /**< number of nodes in the system */
  ITG ne; /**< number of elements in the system */
  ITG nboun; /**< number of boundary conditions (SPCs) in the system */
  ITG nmpc; /**< number of MPCs in the system */
  ITG nforc; /**< number of concentrated forces applied */
  ITG nload; /**< number of load cases to model */
  ITG nprint=0; /**< number of print requests to zero */
  ITG nset; /**< number of elements or node sets in specific set */
  ITG nalset; /**< number of elements or nodes in a specific set */
  ITG nentries=17; /**< number of entries ? */
  ITG nmethod; /**< numerical method */
  ITG neq[3]={0,0,0}; /**< number if equations in different stages of analysis */
  ITG i; /**< general-purpose integer variable */
  ITG mpcfree=1; /**< index for free MPIC slots */
  ITG mei[4]; /**< mechanical element indices used in assembly process */
  ITG j; /**< general-purpose integer variable */
  ITG nzl; /**< number of non-zero entries in sparse matrix  */
  ITG nam; /** number of amplitude definitions in the input */
  ITG nbounold=0; /**< number of boundary conditions from previous step */
  ITG nforcold=0; /**< old number of concentrated forces */ 
  ITG nloadold=0; /**< old number of load cases */
  ITG nbody,nbody_=0; /**< body forces applied in the model */
  ITG nbodyold=0; /**< body forces from previous step */
  ITG network=0; /**< network analysis */
  ITG nheading_=0; /**< heading of a particular input section */
  ITG k; /**< general-purpose integer variable */
  ITG nzs[3]; /**< number of non-zero entries in different stages of a sparse matrix */ 
  ITG nmpc_=0; /**< number of multiple point constraints */
  ITG nload_=0; /**< number of load cases to zero */
  ITG nforc_=0; /**< number of concentrated forces */
  ITG istep; /**< current analysis step number */
  ITG istat; /**< status variable used to track the success or failure for a specific operation */
  ITG nboun_=0; /**< number to boundary conditions to zero */
  ITG nintpoint=0; /**< number of integration points to zero */
  ITG iperturb[2]={0,0}; /**< array indicating whether a perturbation analysis is being conducted */
  ITG nmat; /**< number of materials added in the model */
  ITG ntmat_=0; /**< number of material types to zero */
  ITG norien; /**< number of orientations defined in the model */
  ITG ithermal[2]={0,0}; /**< array indicating whether thermal analysis is being conducted */
  ITG nmpcold; /**< old number of MPCs from previous step */
  ITG iprestr; /**< prestress */
  ITG kode; /**< unknown */
  ITG isolver=0; /**< specified the solver being used */
  ITG nslavs=0; /**< number of slave nodes in contact problems  */
  ITG nkon_=0; /**< total number of entries in k-vector (?) */
  ITG ne0; /**< original nnumber of elements */
  ITG nkon0; /**< original number of entries in k-vector(?) */
  ITG mortar=0; /**< mortar contact */
  ITG jout[2]={1,1}; /**< array for output control */
  ITG nlabel; /**< number of labels in the model */
  ITG nkon=0; /**< number ofd entries in the connectivity vector to zero */
  ITG idrct; /**< direction of loading or displacement constraints */
  ITG jmax[2]; /**< array for storing maximum values of certain parameters */
  ITG iexpl; /**< explicit analysis flag */
  ITG nevtot=0; /**< total number of events or load steps to zero */
  ITG ifacecount=0; /**< counter for contact faces to zero */
  ITG iplas=0; /**< plasticity flag to zero */
  ITG npmat_=0; /**< number of material points to zero */
  ITG mi[3]={0,3,1}; /**< array of ints storing model information such as D.O.Fs */
  ITG ntrans; /**< number of transformations applied in the analysis */
  ITG mpcend=-1; /**< index for the end of MPC list to -1 */
  ITG namtot_=0; /**< total number of amplitide defs to zero */
  ITG iumat=0; /**< flag for material subroutines */
  ITG icascade=0; /**< cascading load or boundary condition definitions to zero */
  ITG maxlenmpc; /**< max. length of the MPC definitions */
  ITG mpcinfo[4]; /**< MPIC information */
  ITG ne1d=0; /**< number of 1D elements in the model to zero */
  ITG ne2d=0; /**< number of 2D elements in the model to zero */
  ITG infree[4]={0,0,0,0}; /**< array indicating free indices for node or element arrays */ 
  ITG callfrommain; /**< flag for subroutine being called from the main program */
  ITG nflow=0; /**< number of flow conditions to zero */
  ITG jin=0; /**< general-purpose integer variable */
  ITG irstrt[2]={0,0}; /**< flag for restart analysis */
  ITG nener=0; /**< number of energy variables to zero */
  ITG jrstrt=0; /**< flag for restart analysis */
  ITG nenerold; /**< number of energy variables from the previous step */
  ITG nline; /**< number of lines in the output or input */
  ITG *ipoinp=NULL; /**< position of input data in the input deck */
  ITG *inp=NULL; /**< input data array */
  ITG ntie; /**< number of tie constraints *(?) */
  ITG ntie_=0; /**< number of tie constraints to zero */
  ITG mcs=0; /**< number ofd constraint sets in the model */
  ITG nprop_=0; /**< number ofd property defintions to zero */
  ITG nprop=0; /**< unknown (double definiton ?) */
  ITG itpamp=0; /**< index for time-amplitude definitions */
  ITG iviewfile; /**< flag for a view file is being used */
  ITG nkold; /**< number of nodes in the previous step of analysis */
  ITG nevdamp_=0; /**< number of viscous damping elements to zero */ 
  ITG npt_=0; /**< number of plasticity related terms to zero */
  ITG cyclicsymmetry; /**< flag for using cyclic symmetry */
  ITG nmethodl; /**< flag for number of methods being used in the analysis */
  ITG iaxial=1; /**< flag for axial symmetry */
  ITG inext=0; /**< next iteration number initialized to zero */
  ITG icontact=0; /**< number of contact definitions to zero */
  ITG nobject=0; /**< number of objects defined in the model to zero */
  ITG nobject_=0; /**< number of objects to zero (?) */ 
  ITG iit=-1; /**< current iteration number, initialized to -1 */
  ITG nzsprevstep[3]; /**< number of non-zero entries in the previous step's matrix */
  ITG memmpcref_; /**< memoery allocated for MPC references */
  ITG mpcfreeref=-1; /**< initializes the free index for MPC referecs to -1 */
  ITG maxlenmpcref; /**< max. length of the MPC references */
  ITG *nodempcref=NULL; /**< node numbers for the MPC references */
  ITG *ikmpcref=NULL; /**< sorted D.O.Fs for MPC references */
  ITG isens=0; /**< flag for sensitivity analysis to zero */
  ITG namtot=0; /**< total number of amplitude definitions to zero */
  ITG nstam=0; /**< number of sate varibles in  the model to zero */
  ITG ndamp=0; /**< number of damping elements in the model */
  ITG nef=0; /**< number of F.E element entities in the model to zero */

  ITG *meminset=NULL; /**< pointer to memeory locations related to element or node sets */
  ITG *rmeminset=NULL; /**< real memory locations associated with element or node sets */

  ITG nzs_;    /**< number of non-zero entries in a sparse matix */
  ITG nk_=0;  /**< number of nodes to zero */
  ITG ne_=0;  /** number of elements to zero */
  ITG nset_=0; /**< number of element or node sets to zero */
  ITG nalset_=0; /**< number of elements or nodes in a set to zero */
  ITG nmat_=0; /**< number of material definitions to zero */
  ITG norien_=0; /**< orientation to zero */
  ITG nam_=0;  /**< number of material definitions to zero */
  ITG ntrans_=0; /**< number ofd transformations to zero */
  ITG ncs_=0; /**< number of constraint sets to zero */
  ITG nstate_=0; /**< number of state variables to zero */
  ITG ncmat_=0; /**< number of constitutive material models to zero */
  ITG memmpc_=0; /**< memoery allocated for MPCs to zero */
  ITG nprint_=0; /**< print flags set to zero*/
  ITG nuel_=0; /**< number of user-defined elements to zero */


  double *co=NULL;   /**< coordinates of all nodes */
  double  *xboun=NULL; /**< values of SPCS applied at nodes */
  double  *coefmpc=NULL; /**< coeffs for MPCs */
  double  *xforc=NULL; /**< magnitudes of concentrated forces applied to nodes */
  double *clearini=NULL; /**< initial clearance values between slave nodes and master surfaces in contact analysis */
  double  *xload=NULL; /**< distriubuted load magnitudes */
  double  *xbounold=NULL; /**< values of SPCs at the start of a step */
  double  *xforcold=NULL; /**< magnitudes of concentrated forces applied to nodes at an old step */
	double *vold=NULL; /**< Displacement vector */
  double  *sti=NULL; /**< Unknown (?) */
  double  *xloadold=NULL; /**< magnitudes of distributed loads at an old step */
  double  *xnor=NULL; /**< Unknown (?) */
  double *reorder=NULL; /**< Unknown (?) */
  double *dcs=NULL; /**< Unknown (?) */
  double  *thickn=NULL; /**< shell element thickness vales */
  double  *thicke=NULL; /**< element thickness values */
  double  *offset=NULL; /**< offsets in elementrs like beams or shells where nodes are not necessarily positioned at the centeriod */
  double *elcon=NULL; /**< material properties related elasticity for the elements */
  double  *rhcon=NULL; /**< material properties related to thermal conductivity */
  double  *alcon=NULL; /**< material properties related to specific material type */
  double  *alzero=NULL; /**< material properties related to thermal expansion for isotropic materials */
  double  *t0=NULL; /**< initial temperature values for nodes or elemenrts */
  double  *t1=NULL; /**< final temperature values at the end of a step */
	double *prestr=NULL; /**< prestress values in elements */
  double  *orab=NULL; /**< orientation information for material anisotropy */
  double  *amta=NULL; /**< amplitude definitions that scale loads and other parameters */
  double *veold=NULL; /**< velocity values from the previous increment */
  double  *accold=NULL; /**< acceleration values from previous increment */
  double  *t1old=NULL; /**< temperature values at the start of a step */
  double  *eme=NULL; /**< strain values at the integration points in each element */
  double  *plicon=NULL; /**< plasticity-related constants for materials undergoing pastic deformation */
  double  *pslavsurf=NULL; /**< slace surface properties in contact analysis */
  double  *plkcon=NULL; /**< plasticity-related constants */
	double *xstate=NULL; /**< internal state variables of the elements */
  double *trab=NULL; /**< coordinates for transformations between global and local coordinate systems */
  double *ener=NULL; /**< energy values at integration points */
  double *shcon=NULL; /**< shell conductivty values */
  double *cocon=NULL; /**< convective heat transfer properties for elements */
  double *cs=NULL; /**< transformation matrices for coordinate systems */
  double *tietol=NULL; /**, tie constraint tolerance */
  double *fmpc=NULL; /**< values of MPC at each increment */
  double *prop=NULL; /**< element property values */
  double *t0g=NULL; /**< global temperature values at start */
  double *t1g=NULL; /**< global temperature values at the end of step */
  double *xbody=NULL; /**< body force values from current step */
  double *xbodyold=NULL; /**< body force values from old step */
  double *coefmpcref=NULL; /**< coefficients of reference MPCs */
  double *dacon=NULL; /**< damping constants for dynamic analysis */
  double *vel=NULL;  /**< velocity values from current step */
  double *velo=NULL; /**< velocity values from a previous step */
  double *veloo=NULL; /**< velocity values from old-old step */
  double *design=NULL; /**< design variables for topology optimization */
  double *rhoPhys=NULL; /**< phyiscal element densities */
  double *stx = NULL;   /**< 1D array for element stress */

  double *gradCompl=NULL;   /**< compliance gradient */
  double *elCompl=NULL; /**<  element complaince */
  double *elCG=NULL; /**< element CG */
  double *eleVol=NULL; /**< element voluime */
  double *designFiltered=NULL; /**< filtered densitites */
  double *gradComplFiltered=NULL; /**< filtered compliance  sensitivities */ 
  double *eleVolFiltered=NULL; /**< filtered element volume sensitivities */
  int *passiveIDs = NULL;
  int numPassive = 0;


  double ctrl[56]={4.5,8.5,9.5,16.5,10.5,4.5,0.,5.5,0.,0.,0.25,0.5,0.75,0.85,0.,0.,1.5,0.,0.005,0.01,0.,0.,0.02,1.e-5,1.e-3,1.e-8,1.e30,1.5,0.25,1.01,1.,1.,5.e-7,5.e-7,5.e-7,5.e-7,5.e-7,5.e-7,5.e-7,-1.,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.e20,1.5,0.5,20.5,1.5,1.5,0.001,0.1,100.5,60.5};

  double fei[3];
  double *xmodal=NULL;
  double timepar[5];

  double  alpha;
  double ttime=0.;
  double qaold[2]={0.,0.};
  double physcon[13]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double pSupplied=0.0;   /**< user-defined penalty paramter */
  double pstiff=0.0;      /**< some penalty paramter */
  double  rmin=0.000000000001; /**< minimum radius */
  double  volfrac=1.00; /**< volume fraction */
  double  qfilter = 3; /**< q-filter value */
  double Pnorm = 0;

  double sigma0 = 0.0;
  double eps_relax = 1e-03;
  double rhomin = 1e-04;
  double pexp = 0.0;

  ITG itertop= 1; /**<iteration counter in topology optimization */
  ITG fnnzassumed = 500; /**< assume 500 non zeros in each row of filtermatrix */ 
  //filternnz=total number of nonzeros in filtermatrix,filternnzElem=no of nonzeros in each row of filtermatrix
  ITG filternnz=0;  /**< actual nnz values in filter matrix */

  ITG *filternnzElems=NULL;

  double *FilterMatrixs=NULL; /**<filter matrix */

  ITG *rowFilters=NULL; /**<row index */
  ITG *colFilters=NULL; /**<column matrix */

  #ifdef CALCULIX_MPI
  MPI_Init(&argc, &argv) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
  #endif

  if(argc==1)
  { 
    /* Inadequate input arguments */
    printf("Usage: Flags: -i jobname -p PENALTY -r RADIUS -f FILTERNNZ --pexp PNORM EXPONENT --sigmin SIGMA MINIMUM --sigrelax STRESS RELAXATION\n");
    FORTRAN(stop,());
  }

else
{
  /* Loop over all input arguments */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"-i")==0)
    {
      /* Input file found, allocate  forjob name */
      strcpy(jobnamec,argv[i+1]);
      strcpy1(jobnamef,argv[i+1],132);
      jin++;
      break;
    }

    if(strcmp1(argv[i],"-v")==0)
    {
      /* Print version */
	    printf("\nThis is CALTOP v1.0: A CalculiX-based large-scale topology optimization framework \n");
      printf(" Authors: Prateek Ranjan, Massachusetts Institute of Technology \n");
      printf(" Authors: Ghanendra Das, Georgia Institute of Technology \n");
      printf(" Authors: Wanzheng Zheng, University of Illinois at Urbana-Champaign \n");
      printf(" Kai. A James, Georgia Institute of Technology \n");
	    FORTRAN(stop,());
    }
  }

  if(jin==0)
  {
    strcpy(jobnamec,argv[1]);strcpy1(jobnamef,argv[1],132);
  }

  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"-o")==0)
    {
      strcpy(output,argv[i+1]);
      break;
    }
  }

  /* Read penalty parameter */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"-p")==0) 
    {
      pSupplied=atof(argv[i+1]);
      break;
    }
  }

  /* Read filter radius */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"-r")==0)
    {
      rmin=atof(argv[i+1]);
      break;
    }
  }

  /* Read P-norm exponent */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"--pexp")==0) 
    {
      pexp=atof(argv[i+1]);
      break;
    }
  }  

  /* Read Sigmal */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"--sigmin")==0) 
    {
      sigma0=atof(argv[i+1]);
      break;
    }
  }

  /* Stress relaxation factor */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"--sigrelax")==0) 
    {
      eps_relax=atof(argv[i+1]);
      break;
    }
  }

  /* Read expected non-zeros in filter matrix */
  for(i=1; i<argc; i++)
  {
    if(strcmp1(argv[i],"-f")==0) 
    {
      fnnzassumed=atoi(argv[i+1]);
      break;
    }
  }

}

//printf("\n */*/*/* rmin= %f\n",rmin);

/* Assign default value for penalty parameter */
if (pSupplied==0.0) 
  {
    pstiff=1.0;
  }
else 
  {
    pstiff=pSupplied;
  }

//printf("\n*/*/*/*pstiff= %f\n",pstiff);


//setenv("CCX_JOBNAME_GETJOBNAME",jobnamec,1);
putenv("CCX_JOBNAME_GETJOBNAME=jobnamec");

#ifdef BAM
ITG lop=0,lrestart=0,kstep=1,kinc=1;
double time[2],dtime;
FORTRAN(uexternaldb,(&lop,&lrestart,time,&dtime,&kstep,&kinc));
#endif

//FORTRAN(openfile,(jobnamef,output));

printf("\n");

printf("  #####    #####   ##        #######  #####   ######  \n");
printf(" ##       ##   ##  ##           ##   ##   ##  ##   ## \n");
printf(" ##       #######  ##           ##   ##   ##  ######  \n");
printf(" ##       ##   ##  ##           ##   ##   ##  ##      \n");
printf("  #####   ##   ##  #######      ##    #####   ##      \n");
    
printf("\n");

printf("* Contributors:\n");
printf("* Prateek Ranjan, Dept. of Aeronautics & Astronautics,\n");
printf("* Massachusetts Institute of Technology \n");
printf("* Wanzheng Zheng, Dept. of Aerospace Engineering,\n");
printf("* University of Illinois at Urbana Champaign \n");
printf("* Ghanendra Kumar Das, Dept. of Aerospace Engineering,\n");
printf("* Georgia Institute of Technology\n");
printf("* Kai A. James, Dept. of Aerospace Engineering,\n");
printf("* Georgia Institute of Technology\n");
printf("************************************************************\n\n");
printf("Original FEA source code: CalculiX by Guido Dhondt\n");
printf("************************************************************\n\n");



istep=0;
istat=0;
iprestr=0;
kode=0;

/* default solver */

#if defined(SGI)
 isolver=4;
#elif defined(PARDISO)
 isolver=7;
#elif defined(SPOOLES)
 isolver=0;
#elif defined(TAUCS)
 isolver=5;
#else
 isolver=3;
#endif

NNEW(ipoinp,ITG,2*nentries);

/* conservative estimate of the fields to be allocated */

readinput(jobnamec,&inpc,&nline,&nset_,ipoinp,&inp,&ipoinpc,ithermal,&nuel_);

NNEW(set,char,81*nset_);
NNEW(meminset,ITG,nset_);
NNEW(rmeminset,ITG,nset_);
NNEW(iuel,ITG,4*nuel_);


FORTRAN(allocation,(&nload_,&nforc_,&nboun_,&nk_,&ne_,&nmpc_,&nset_,&nalset_,
   &nmat_,&ntmat_,&npmat_,&norien_,&nam_,&nprint_,mi,&ntrans_,
   set,meminset,rmeminset,&ncs_,&namtot_,&ncmat_,&memmpc_,&ne1d,
   &ne2d,&nflow,jobnamec,irstrt,ithermal,&nener,&nstate_,&istep,
   inpc,ipoinp,inp,&ntie_,&nbody_,&nprop_,ipoinpc,&nevdamp_,&npt_,&nslavs,
   &nkon_,&mcs,&mortar,&ifacecount,&nintpoint,infree,&nheading_,&nobject_,
   iuel,&iprestr,&nstam,&ndamp,&nef));


SFREE(set);
SFREE(meminset);
SFREE(rmeminset);
mt=mi[1]+1;
NNEW(heading,char,66*nheading_);

nzs_=20000000;

nload=0;
nbody=0;
nforc=0;
nboun=0;
nk=0;
nmpc=0;
nam=0;

/* caveat: if changing next line:
   - change noelfiles appropriately
   - change nlabel in geomview.f, expand.c, storeresidual.f
     and createmddof.f
   - change the dimension of label in geomview.f
   - change the documentation (tex-file)  */

nlabel=48;

/* While the status variable istat is non-negative */
while(istat>=0)
{
  //printf("iStep at begining of outer loop: %d ", istat);
  fflush(stdout);

  /* in order to reduce the number of variables to be transferred to
     the subroutines, the max. field sizes are (for most fields) copied
     into the real sizes */

  nzs[1]=nzs_;
  //  nprint=nprint_;

  if((istep==0)||(irstrt[0]<0)) 
  {
    ne=ne_;
    nset=nset_;
    nalset=nalset_;
    nmat=nmat_;
    norien=norien_;
    ntrans=ntrans_;
    ntie=ntie_;
  //    nobject=nobject_;

    /* allocating space before the first step */

    /* coordinates and topology */

    NNEW(co,double,3*nk_);
    NNEW(kon,ITG,nkon_);
    NNEW(ipkon,ITG,ne_);
    NNEW(lakon,char,8*ne_);
    NNEW(design,double,ne_);



    /* property cards */

    if(nprop_>0)
    {
	    NNEW(ielprop,ITG,ne_);

	    for(i=0;i<ne_;i++) ielprop[i]=-1;

	    NNEW(prop,double,nprop_);
    }

    /* fields for 1-D and 2-D elements */

    if((ne1d!=0)||(ne2d!=0))
    {
	    NNEW(iponor,ITG,2*nkon_);
	    for(i=0;i<2*nkon_;i++) iponor[i]=-1;
	    NNEW(xnor,double,36*ne1d+24*ne2d);
	    NNEW(knor,ITG,24*(ne1d+ne2d)*(mi[2]+1));
	    NNEW(thickn,double,2*nk_);
	    NNEW(thicke,double,mi[2]*nkon_);
	    NNEW(offset,double,2*ne_);
	    NNEW(iponoel,ITG,nk_);
	    NNEW(inoel,ITG,9*ne1d+24*ne2d);
	    NNEW(rig,ITG,nk_);

	    if(infree[2]==0)infree[2]=1;
    }

    /* SPC's */

    NNEW(nodeboun,ITG,nboun_);
    NNEW(ndirboun,ITG,nboun_);
    NNEW(typeboun,char,nboun_+1);

    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamboun,ITG,nboun_);

    NNEW(xboun,double,nboun_);
    NNEW(ikboun,ITG,nboun_);
    NNEW(ilboun,ITG,nboun_);

    /* MPC's */

    NNEW(ipompc,ITG,nmpc_);
    NNEW(nodempc,ITG,3*memmpc_);

    for(i=0;i<3*memmpc_;i+=3)
    {
      nodempc[i+2]=i/3+2;
    }

    nodempc[3*memmpc_-1]=0;
    NNEW(coefmpc,double,memmpc_);
    NNEW(labmpc,char,20*nmpc_+1);
    NNEW(ikmpc,ITG,nmpc_);
    NNEW(ilmpc,ITG,nmpc_);
    NNEW(fmpc,double,nmpc_);

    /* nodal loads */

    NNEW(nodeforc,ITG,2*nforc_);
    NNEW(ndirforc,ITG,nforc_);

    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamforc,ITG,nforc_);
    
    NNEW(idefforc,ITG,nforc_);
    NNEW(xforc,double,nforc_);
    NNEW(ikforc,ITG,nforc_);
    NNEW(ilforc,ITG,nforc_);

    /* distributed facial loads */

    NNEW(nelemload,ITG,2*nload_);

    if((istep == 0)||((irstrt[0]<0)&&(nam_>0)))NNEW(iamload,ITG,2*nload_);
    
    NNEW(idefload,ITG,nload_);
    NNEW(sideload,char,20*nload_);
    NNEW(xload,double,2*nload_);

    /* distributed volumetric loads */

    NNEW(cbody,char,81*nbody_);
    NNEW(idefbody,ITG,nbody_);
    NNEW(ibody,ITG,3*nbody_);
    NNEW(xbody,double,7*nbody_);
    NNEW(xbodyold,double,7*nbody_);

    /* printing output */

    NNEW(prlab,char,6*nprint_);
    NNEW(prset,char,81*nprint_);

    /* set definitions */

    NNEW(set,char,81*nset);
    NNEW(istartset,ITG,nset);
    NNEW(iendset,ITG,nset);
    NNEW(ialset,ITG,nalset);

    /* (hyper)elastic constants */

    NNEW(elcon,double,(ncmat_+1)*ntmat_*nmat);
    NNEW(nelcon,ITG,2*nmat);

    /* density */

    NNEW(rhcon,double,2*ntmat_*nmat);
    NNEW(nrhcon,ITG,nmat);

    /* damping */

    if(ndamp>0)
    {
      NNEW(dacon,double,nmat);
    }

    /* specific heat */

    NNEW(shcon,double,4*ntmat_*nmat);
    NNEW(nshcon,ITG,nmat);

    /* thermal expansion coefficients */

    NNEW(alcon,double,7*ntmat_*nmat);
    NNEW(nalcon,ITG,2*nmat);
    NNEW(alzero,double,nmat);

    /* conductivity */

    NNEW(cocon,double,7*ntmat_*nmat);
    NNEW(ncocon,ITG,2*nmat);

    /* isotropic and kinematic hardening coefficients*/

    if(npmat_>0)
    {
	    NNEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	    NNEW(nplicon,ITG,(ntmat_+1)*nmat);
	    NNEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
	    NNEW(nplkcon,ITG,(ntmat_+1)*nmat);
    }

    /* linear dynamic properties */

    NNEW(xmodal,double,11+nevdamp_);
    xmodal[10]=nevdamp_+0.5;

    /* internal state variables (nslavs is needed for restart
       calculations) */

    if(mortar!=1)
    {
	    NNEW(xstate,double,nstate_*mi[0]*(ne+nslavs));
	    nxstate=nstate_*mi[0]*(ne+nslavs);
    }
    else if(mortar==1)
    {
	    NNEW(xstate,double,nstate_*mi[0]*(ne+nintpoint));
	    nxstate=nstate_*mi[0]*(ne+nintpoint);
    }

    /* material orientation */

    if((istep == 0)||((irstrt[0]<0)&&(norien>0))) 
    {
	    NNEW(orname,char,80*norien);
	    NNEW(orab,double,7*norien);
	    NNEW(ielorien,ITG,mi[2]*ne_);
    }

    /* transformations */

    if((istep == 0)||((irstrt[0]<0)&&(ntrans>0))) 
    {
	    NNEW(trab,double,7*ntrans);
	    NNEW(inotr,ITG,2*nk_);
    }

    /* amplitude definitions */

    if((istep == 0)||((irstrt[0]<0)&&(nam_>0))) 
    {
	    NNEW(amname,char,80*nam_);
	    NNEW(amta,double,2*namtot_);
	    NNEW(namta,ITG,3*nam_);
    }

    if((istep == 0)||((irstrt[0]<0)&&(ithermal[0]>0))) 
    {
	    NNEW(t0,double,nk_);
	    NNEW(t1,double,nk_);

	    if((ne1d!=0)||(ne2d!=0))
      {
	      NNEW(t0g,double,2*nk_);
	      NNEW(t1g,double,2*nk_);
	    }
    }

    /* the number in next line is NOT 1.2357111317 -> points
       to user input; instead it is a generic nonzero
       initialization */

    if(istep==0)
    {
	    DMEMSET(t0,0,nk_,1.2357111319);
	    DMEMSET(t1,0,nk_,1.2357111319);
    }

    if((istep == 0)||((irstrt[0]<0)&&(ithermal[0]>0)&&(nam_>0)))
    {
      NNEW(iamt1,ITG,nk_);
    }

    if((istep==0)||((irstrt[0]<0)&&(iprestr>0)))
    {
      NNEW(prestr,double,6*mi[0]*ne_);
    }

    NNEW(vold,double,mt*nk_);
    NNEW(veold,double,mt*nk_);

    /* CFD-results */

    NNEW(vel,double,8*nef);
    NNEW(velo,double,8*nef);
    NNEW(veloo,double,8*nef);

    NNEW(ielmat,ITG,mi[2]*ne_);

    NNEW(matname,char,80*nmat);

    NNEW(filab,char,87*nlabel);

    /* tied constraints */

    if(ntie_>0)
    {
      NNEW(tieset,char,243*ntie_);
      NNEW(tietol,double,3*ntie_);
      NNEW(cs,double,17*ntie_);
    }

    /* objectives for sensitivity analysis */

    if(nobject_>0)
    {
      NNEW(objectset,char,324*nobject_);
      for(i=0;i<324*nobject_;i++){objectset[i]=' ';}
    }

    /* temporary fields for cyclic symmetry calculations */

    if((ncs_>0)||(npt_>0))
    {
      if(2*npt_>24*ncs_)
      {
	      NNEW(ics,ITG,2*npt_);
      }
      else
      {
	      NNEW(ics,ITG,24*ncs_);
      }
      if(npt_>30*ncs_)
      {
	      NNEW(dcs,double,npt_);
      }
      else
      {
	      NNEW(dcs,double,30*ncs_);
      }
    }

    /* slave faces */

    NNEW(islavsurf,ITG,2*ifacecount+2);

  }
  else 
  {

    /* allocating and reallocating space for subsequent steps */

    if((nmethod!=4)&&(nmethod!=5)&&(nmethod!=8)&&(nmethod!=9)&&
       ((abs(nmethod)!=1)||(iperturb[0]<2)))
       {
          NNEW(veold,double,mt*nk_);
       }
    else
    {
      RENEW(veold,double,mt*nk_);
      DMEMSET(veold,mt*nk,mt*nk_,0.);
    }

    RENEW(vold,double,mt*nk_);
    DMEMSET(vold,mt*nk,mt*nk_,0.);

    RENEW(nodeboun,ITG,nboun_);
    RENEW(ndirboun,ITG,nboun_);
    RENEW(typeboun,char,nboun_+1);
    RENEW(xboun,double,nboun_);
    RENEW(ikboun,ITG,nboun_);
    RENEW(ilboun,ITG,nboun_);

    RENEW(nodeforc,ITG,2*nforc_);
    RENEW(ndirforc,ITG,nforc_);
    NNEW(idefforc,ITG,nforc_);
    RENEW(xforc,double,nforc_);
    RENEW(ikforc,ITG,nforc_);
    RENEW(ilforc,ITG,nforc_);

    RENEW(nelemload,ITG,2*nload_);
    NNEW(idefload,ITG,nload_);
    RENEW(sideload,char,20*nload_);
    RENEW(xload,double,2*nload_);

    RENEW(cbody,char,81*nbody_);
    NNEW(idefbody,ITG,nbody_);
    RENEW(ibody,ITG,3*nbody_);
    RENEW(xbody,double,7*nbody_);
    RENEW(xbodyold,double,7*nbody_);

    for(i=7*nbodyold;i<7*nbody_;i++) xbodyold[i]=0;

    if(nam > 0) 
    {
      RENEW(iamforc,ITG,nforc_);
      RENEW(iamload,ITG,2*nload_);
      RENEW(iamboun,ITG,nboun_);
      RENEW(amname,char,80*nam_);
      RENEW(amta,double,2*namtot_);
      RENEW(namta,ITG,3*nam_);
    }

    RENEW(ipompc,ITG,nmpc_);

    RENEW(labmpc,char,20*nmpc_+1);
    RENEW(ikmpc,ITG,nmpc_);
    RENEW(ilmpc,ITG,nmpc_);
    RENEW(fmpc,double,nmpc_);

    if(ntrans > 0)
    {
      RENEW(inotr,ITG,2*nk_);DMEMSET(inotr,2*nk,2*nk_,0);
    }

    RENEW(co,double,3*nk_);DMEMSET(co,3*nk,3*nk_,0.);

    if(ithermal[0] != 0)
    {
	    RENEW(t0,double,nk_);DMEMSET(t0,nk,nk_,0.);
	    RENEW(t1,double,nk_);DMEMSET(t1,nk,nk_,0.);

	    if((ne1d!=0)||(ne2d!=0))
      {
	      RENEW(t0g,double,2*nk_);DMEMSET(t0g,2*nk,2*nk_,0.);
	      RENEW(t1g,double,2*nk_);DMEMSET(t1g,2*nk,2*nk_,0.);
	    }

	    if(nam > 0) 
      {
        RENEW(iamt1,ITG,nk_);
      }
    }
  } //end else

  /* allocation of fields in the restart file */

  if(irstrt[0]<0)
  {
    NNEW(nodebounold,ITG,nboun_);
    NNEW(ndirbounold,ITG,nboun_);
    NNEW(xbounold,double,nboun_);
    NNEW(xforcold,double,nforc_);
    NNEW(xloadold,double,2*nload_);

    if(ithermal[0]!=0) NNEW(t1old,double,nk_);

    NNEW(sti,double,6*mi[0]*ne);
    NNEW(eme,double,6*mi[0]*ne);

    if(nener==1)NNEW(ener,double,mi[0]*ne*2);
    if(mcs>ntie_) RENEW(cs,double,17*mcs);

    if(mortar==1)
    {
	    NNEW(pslavsurf,double,3*nintpoint);
	    NNEW(clearini,double,3*9*ifacecount);
    }
  }

  nenerold=nener;
  nkold=nk;

  /* opening the eigenvalue file and checking for cyclic symmetry */
  strcpy(fneig,jobnamec);
  strcat(fneig,".eig");
  cyclicsymmetry=0;

  if((f1=fopen(fneig,"rb"))!=NULL)
  {
    if(fread(&cyclicsymmetry,sizeof(ITG),1,f1)!=1)
    {
	    printf("*ERROR reading the information whether cyclic symmetry is involved in the eigenvalue file");
	    exit(0);
    }
    fclose(f1);
  }

  nmpcold=nmpc;

  /* reading the input file */
  if(istep==0)mortar=-1;


  FORTRAN(calinput,(co,&nk,kon,ipkon,lakon,&nkon,&ne,
            nodeboun,ndirboun,xboun,&nboun,
	    ipompc,nodempc,coefmpc,&nmpc,&nmpc_,nodeforc,ndirforc,xforc,&nforc,
	    &nforc_,nelemload,sideload,xload,&nload,&nload_,
	    &nprint,prlab,prset,&mpcfree,&nboun_,mei,set,istartset,iendset,
	    ialset,&nset,&nalset,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
	    alzero,t0,t1,matname,ielmat,orname,orab,ielorien,amname,
            amta,namta,&nam,&nmethod,iamforc,iamload,iamt1,
	    ithermal,iperturb,&istat,&istep,&nmat,&ntmat_,&norien,prestr,
	    &iprestr,&isolver,fei,veold,timepar,
	    xmodal,filab,jout,&nlabel,&idrct,
	    jmax,&iexpl,&alpha,iamboun,plicon,nplicon,
	    plkcon,nplkcon,&iplas,&npmat_,mi,&nk_,trab,inotr,&ntrans,
	    ikboun,ilboun,ikmpc,ilmpc,ics,dcs,&ncs_,&namtot_,cs,&nstate_,
	    &ncmat_,&iumat,&mcs,labmpc,iponor,xnor,knor,thickn,thicke,
	    ikforc,ilforc,offset,iponoel,inoel,rig,infree,nshcon,shcon,
            cocon,ncocon,physcon,&nflow,
            ctrl,&maxlenmpc,&ne1d,&ne2d,&nener,vold,nodebounold,
            ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
            sti,ener,xstate,jobnamec,irstrt,&ttime,
            qaold,output,typeboun,inpc,ipoinp,inp,tieset,tietol,
            &ntie,fmpc,cbody,ibody,xbody,&nbody,&nbody_,xbodyold,&nam_,
	    ielprop,&nprop,&nprop_,prop,&itpamp,&iviewfile,ipoinpc,
	    &nslavs,t0g,t1g,&network,&cyclicsymmetry,idefforc,idefload,
	    idefbody,&mortar,&ifacecount,islavsurf,pslavsurf,clearini,
	    heading,&iaxial,&nobject,objectset,&nprint_,iuel,&nuel_,
	    nodempcref,coefmpcref,ikmpcref,&memmpcref_,&mpcfreeref,
	    &maxlenmpcref,&memmpc_,&isens,&namtot,&nstam,dacon,vel,&nef,
	    velo,veloo));

    struct stat buffer;

    if (stat("skinElementList.nam", &buffer) != 0) 
    {
      printf("File 'skinElementList.nam' not found.\n");
    } 
    else
    /* Read the element indies from skinElementList.nam */ 
    {
      passiveIDs = passiveElements("skinElementList.nam", &numPassive);
      printf("File 'skinElementList.nam' found. \n");
      printf("Read %d passive elements.\n", numPassive);

      printf("First five skin element IDs => \n");
      for (int i = 0; i < 5; i++) 
      {
        printf("Passive Element ID: %d\n", passiveIDs[i]);
      }

      printf(".\n");
      printf(".\n");
      printf(".\n");

      printf("Last five skin element IDs => \n");
      for (int i = numPassive - 5; i < numPassive; i++) 
      {
        printf("Passive Element ID: %d\n", passiveIDs[i]);
      }

    }

    printf("\n");


  /* Read element desitiies from .dat file, if absent, initialize the design to one */    
  rho(design,ne);

  /* Set thhe element density of passive elements to 1 */
  filterOutPassiveElems_density(design, ne, passiveIDs, numPassive);

  /* Define stress array here, pass to linstatic and then plot in write_vtu.c */
  NNEW(stx,double,6*mi[0]*ne);

  /* Define sensitivity vector for Pnorm and pass to linstatic()*/
  double *dPnorm_drho = NULL;
  double *dPnorm_drhoFiltered = NULL;
  NNEW(dPnorm_drho, double, ne);


  #ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT
  for(i=0;i!=nmat;++i)
  {
    calculix_registerExternalBehaviour(matname+80*i);
  }
  #endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */

  if((istep==1)&&(mortar==-1))
  {
    mortar=0;
  }
  else
  {
    icontact=1;
  }

  nload0=nload;
  SFREE(idefforc);
  SFREE(idefload);
  SFREE(idefbody);

  /*  <--- Use this for debugging
  printf("Sleeping noe...\n");
  sleep(5);
  printf("sleep over!");
  export OMP_NUM_THREADS=4  -< copy to terminal
  */
  if(nheading_>=0)
  {
     /* do nothing!
      writeheading(jobnamec,heading,&nheading_);
      SFREE(heading);
      nheading_=-1;

      */
  }



  if((abs(nmethod)!=1)||(iperturb[0]<2))icascade=0;

  //	FORTRAN(writeboun,(nodeboun,ndirboun,xboun,typeboun,&nboun));

  if(istat<0) break;

  /* first step of F.E analysis */
  if(istep == 1) 
  {
    SFREE(iuel);

    /* tied contact constraints: generate appropriate MPC's */

    tiedcontact(&ntie, tieset, &nset, set,istartset, iendset, ialset,
       lakon, ipkon, kon,tietol,&nmpc, &mpcfree, &memmpc_,
       &ipompc, &labmpc, &ikmpc, &ilmpc,&fmpc, &nodempc, &coefmpc,
       ithermal, co, vold,&nef,&nmpc_,mi,&nk,&istep,ikboun,&nboun,
       kind1,kind2);

    /* reallocating space in the first step */

    /* allocating and initializing fields pointing to the previous step */

    RENEW(vold,double,mt*nk);
    NNEW(sti,double,6*mi[0]*ne);

    /* strains */

    NNEW(eme,double,6*mi[0]*ne);

    /* residual stresses/strains */

    if(iprestr==1) 
    {
      printf("CAUTION: Structure is pre-stressed! \n");
	    RENEW(prestr,double,6*mi[0]*ne);

	    for(i=0;i<ne;i++)
      {
	      for(j=0;j<mi[0];j++)
        {
		      for(k=0;k<6;k++)
          {
		        sti[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		      }
	      }
	    }
    }
    else if(iprestr==2)
    {
	    RENEW(prestr,double,6*mi[0]*ne);

	    for(i=0;i<ne;i++)
      {
	      for(j=0;j<mi[0];j++)
        {
		      for(k=0;k<6;k++)
          {
		        eme[6*mi[0]*i+6*j+k]=prestr[6*mi[0]*i+6*j+k];
		      }
	      }
	    }
    }

    /* Structure is not pre-stressed, free variable */
    else 
    {
	    SFREE(prestr);
    }

    NNEW(nodebounold,ITG,nboun);
    NNEW(ndirbounold,ITG,nboun);
    NNEW(xbounold,double,nboun);
    NNEW(xforcold,double,nforc);
    NNEW(xloadold,double,2*nload);

    /* initial temperatures: store in the "old" boundary conditions */

    if(ithermal[0]>1)
    {
	    for(i=0;i<nboun;i++)
      {
	      if(strcmp1(&typeboun[i],"F")==0) continue;
	      if(ndirboun[i]==0)
        {
		      xbounold[i]=vold[mt*(nodeboun[i]-1)];
	      }
	    }
    }

    /* initial temperatures: store in the "old" temperature field */
    if(ithermal[0]!=0)
    {
      NNEW(t1old,double,nk);
      for(i=0;i<nk;i++) t1old[i]=t0[i];
    }

    /* element definition */

    RENEW(kon,ITG,nkon);
    RENEW(ipkon,ITG,ne);
    RENEW(lakon,char,8*ne);

    /* property cards */

    if(nprop>0)
    {
	    RENEW(ielprop,ITG,ne);
	    RENEW(prop,double,nprop);
    }
    else
    {
	    SFREE(ielprop);SFREE(prop);
    }

    /* fields for 1-D and 2-D elements */
    if((ne1d!=0)||(ne2d!=0))
    {
	    RENEW(iponor,ITG,2*nkon);
	    RENEW(xnor,double,infree[0]);
	    RENEW(knor,ITG,infree[1]);
	    SFREE(thickn);
	    RENEW(thicke,double,mi[2]*nkon);
	    RENEW(offset,double,2*ne);
	    RENEW(inoel,ITG,3*(infree[2]-1));
	    RENEW(iponoel,ITG,infree[3]);
	    RENEW(rig,ITG,infree[3]);
    }

    /* set definitions */

    RENEW(set,char,81*nset);
    RENEW(istartset,ITG,nset);
    RENEW(iendset,ITG,nset);
    RENEW(ialset,ITG,nalset);

    /* material properties */

    RENEW(elcon,double,(ncmat_+1)*ntmat_*nmat);
    RENEW(nelcon,ITG,2*nmat);

    RENEW(rhcon,double,2*ntmat_*nmat);
    RENEW(nrhcon,ITG,nmat);

    if(ndamp>0)
    {
      RENEW(dacon,double,nmat);
    }

    RENEW(shcon,double,4*ntmat_*nmat);
    RENEW(nshcon,ITG,nmat);

    RENEW(cocon,double,7*ntmat_*nmat);
    RENEW(ncocon,ITG,2*nmat);

    RENEW(alcon,double,7*ntmat_*nmat);
    RENEW(nalcon,ITG,2*nmat);
    RENEW(alzero,double,nmat);

    RENEW(matname,char,80*nmat);
    RENEW(ielmat,ITG,mi[2]*ne);

    /* allocating space for the state variables */

    if(mortar!=1)
    {
  
	    RENEW(xstate,double,nstate_*mi[0]*(ne+nslavs));

	    for(i=nxstate;i<nstate_*mi[0]*(ne+nslavs);i++)
      {
        xstate[i]=0.;
      }
    }
    else if(mortar==1)
    {
	    RENEW(xstate,double,nstate_*mi[0]*(ne+nintpoint));

	    for(i=nxstate;i<nstate_*mi[0]*(ne+nintpoint);i++)
      {
        xstate[i]=0.;
      }
    }

    /* next statements for plastic materials and nonlinear springs */
    if(npmat_>0)
    {
	    RENEW(plicon,double,(2*npmat_+1)*ntmat_*nmat);
	    RENEW(nplicon,ITG,(ntmat_+1)*nmat);
	    RENEW(plkcon,double,(2*npmat_+1)*ntmat_*nmat);
	    RENEW(nplkcon,ITG,(ntmat_+1)*nmat);
    }

    /* material orientation */
    if(norien > 0) 
    {
      RENEW(orname,char,80*norien);
      RENEW(ielorien,ITG,mi[2]*ne);
      RENEW(orab,double,7*norien);
    }
    else 
    {
	    SFREE(orname);
	    SFREE(ielorien);
	    SFREE(orab);
    }

    /* amplitude definitions */
    if(nam > 0) 
    {
      RENEW(amname,char,80*nam);
      RENEW(namta,ITG,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }
    else 
    {
      SFREE(amname);
      SFREE(amta);
      SFREE(namta);
      SFREE(iamforc);
      SFREE(iamload);
      SFREE(iamboun);
    }

    if(ntrans > 0)
    {
      RENEW(trab,double,7*ntrans);
    }
    else
    {
      SFREE(trab);
      SFREE(inotr);
    }

    if(ithermal[0] == 0)
    {
      SFREE(t0);
      SFREE(t1);
      SFREE(t0g);
      SFREE(t1g);
    }

    if((ithermal[0] == 0)||(nam<=0))
    {
      SFREE(iamt1);
    }

    if(ncs_>0)
    {
      RENEW(ics,ITG,ncs_);
      //      SFREE(dcs);
    }
    else if(npt_>0)
    {
      SFREE(ics);
    }

    SFREE(dcs);

    if(mcs>0)
    {
	    RENEW(cs,double,17*mcs);
    }
    else
    {
	    SFREE(cs);
    }

  } // end if(istep == 1)
  else
  {
    printf("CAUTION: Analysis is in second step. Not valid for static analysis!");

    /* reallocating space in all but the first step (>1) */
    RENEW(vold,double,mt*nk);

    /* if the SPC boundary conditions were changed in the present step,
       they have to be rematched with those in the last step. Removed SPC
       boundary conditions do not appear any more (this is different from
       forces and loads, where removed forces or loads are reset to zero;
       a removed SPC constraint does not have a numerical value any more) */

    NNEW(reorder,double,nboun);
    NNEW(nreorder,ITG,nboun);

    if(nbounold<nboun)
    {
      RENEW(xbounold,double,nboun);
      RENEW(nodebounold,ITG,nboun);
      RENEW(ndirbounold,ITG,nboun);
    }

    FORTRAN(spcmatch,(xboun,nodeboun,ndirboun,&nboun,xbounold,nodebounold,
		      ndirbounold,&nbounold,ikboun,ilboun,vold,reorder,nreorder,
                      mi));
    RENEW(xbounold,double,nboun);
    RENEW(nodebounold,ITG,nboun);
    RENEW(ndirbounold,ITG,nboun);
    SFREE(reorder); SFREE(nreorder);

    /* for additional forces or loads in the present step, the
       corresponding slots in the force and load fields of the
       previous steps are initialized */

    RENEW(xforcold,double,nforc);

    for(i=nforcold;i<nforc;i++) xforcold[i]=0;

    RENEW(xloadold,double,2*nload);

    for(i=2*nloadold;i<2*nload;i++) xloadold[i]=0;

    if(ithermal[0]!=0)
    {
      RENEW(t1old,double,nk);
    }

    if(nam > 0) 
    {
      RENEW(amname,char,80*nam);
      RENEW(namta,ITG,3*nam);
      RENEW(amta,double,2*namta[3*nam-2]);
    }

  } /* End of istep > 1 */

  /* reallocating fields for all steps (>=1) */
  RENEW(co,double,3*nk);

  RENEW(nodeboun,ITG,nboun);
  RENEW(ndirboun,ITG,nboun);
  RENEW(typeboun,char,nboun+1);
  RENEW(xboun,double,nboun);
  RENEW(ikboun,ITG,nboun);
  RENEW(ilboun,ITG,nboun);

  RENEW(nodeforc,ITG,2*nforc);
  RENEW(ndirforc,ITG,nforc);
  RENEW(xforc,double,nforc);
  RENEW(ikforc,ITG,nforc);
  RENEW(ilforc,ITG,nforc);

  /* temperature loading */

  if(ithermal[0] != 0)
  {
      RENEW(t0,double,nk);
      RENEW(t1,double,nk);

      if((ne1d!=0)||(ne2d!=0))
      {
	      RENEW(t0g,double,2*nk);
	      RENEW(t1g,double,2*nk);
      }

      if(nam > 0)
      {
        RENEW(iamt1,ITG,nk);
      }
  }

  RENEW(nelemload,ITG,2*nload);
  RENEW(sideload,char,20*nload);
  RENEW(xload,double,2*nload);

  RENEW(cbody,char,81*nbody);
  RENEW(ibody,ITG,3*nbody);
  RENEW(xbody,double,7*nbody);
  RENEW(xbodyold,double,7*nbody);

  RENEW(ipompc,ITG,nmpc);
  RENEW(labmpc,char,20*nmpc+1);
  RENEW(ikmpc,ITG,nmpc);
  RENEW(ilmpc,ITG,nmpc);
  RENEW(fmpc,double,nmpc);

  /* energy */

  if((nener==1)&&(nenerold==0))
  {
    NNEW(ener,double,mi[0]*ne*2);

    if((istep>1)&&(iperturb[0]>1))
    {
      printf(" *ERROR in CalculiX: in nonlinear calculations\n");
      printf("        energy output must be selected in the first step\n\n");
      FORTRAN(stop,());
    }
  }

  /* following segment executed if solving:

      Linear dynamic analysis
      Steady state dynamics
      Magnetostatics
      Static analysis with material and geometric non-linearities
  */

  /* initial velocities and accelerations */
  if((nmethod==4)||(nmethod==5)||(nmethod==8)||(nmethod==9)||
     ((abs(nmethod)==1)&&(iperturb[0]>=2)))
     {
        RENEW(veold,double,mt*nk);
     }

  /* not being used so free memory */
  //else 
  //{
  //  SFREE(veold);
 // }

  /* if solving linear dynamic system with geometric non-linearities */

  if((nmethod == 4)&&(iperturb[0]>1)) 
  {
    NNEW(accold,double,mt*nk);
  }

  /* if non-zero amplitude definition */
  if(nam > 0) 
  {
    RENEW(iamforc,ITG,nforc);
    RENEW(iamload,ITG,2*nload);
    RENEW(iamboun,ITG,nboun);
  }

  /* generate force convection elements */

  //  if(network==1){
  if(network>0)  /* Not used in CalTop */
  {
    ne0=ne;nkon0=nkon;nload1=nload;
    RENEW(ipkon,ITG,ne+nload);
    RENEW(lakon,char,8*(ne+nload));
    RENEW(kon,ITG,nkon+9*nload);
    NNEW(inodesd,ITG,nk);
    RENEW(nelemload,ITG,4*nload);
    RENEW(sideload,char,40*nload);

    FORTRAN(genadvecelem,(inodesd,ipkon,&ne,lakon,kon,&nload,
			    sideload,nelemload,&nkon,&network));

    SFREE(inodesd);
    RENEW(ipkon,ITG,ne);
    RENEW(lakon,char,8*ne);
    RENEW(kon,ITG,nkon);
    RENEW(sti,double,6*mi[0]*ne);
    RENEW(eme,double,6*mi[0]*ne);

    if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
    if(nprop>0) RENEW(ielprop,ITG,ne);
    if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);

    RENEW(nelemload,ITG,2*nload);
    RENEW(sideload,char,20*nload);
    RENEW(xload,double,2*nload);
    RENEW(xloadold,double,2*nload);

    if(nam>0)
    {
	    RENEW(iamload,ITG,2*nload);

	    for(i=2*nload1;i<2*nload;i++)iamload[i]=0;
    }

    if(nener==1)RENEW(ener,double,mi[0]*ne*2);
    if(norien>0)RENEW(ielorien,ITG,mi[2]*ne);
    RENEW(ielmat,ITG,mi[2]*ne);
    for(i=mi[2]*ne0;i<mi[2]*ne;i++)ielmat[i]=1;
  } // end network > 0

  if(ntrans > 0)
  {
    RENEW(inotr,ITG,2*nk);
  }

  /*   calling the user routine ufaceload (can be empty) */

  if(ithermal[1]>=2) /* Not used in CalTop */
  {

    NNEW(sideloadtemp,char,20*nload);
    for(i=0;i<nload;i++)
    {
	    strcpy1(&sideloadtemp[20*i],&sideload[20*i],20);
	    if((strcmp1(&sideload[20*i]," ")==0)&&
	     (strcmp1(&sideload[20*i+1]," ")!=0))
       {
	        strcpy1(&sideloadtemp[20*i],"F",1);
	      }
    }

    FORTRAN(ufaceload,(co,ipkon,kon,lakon,&nboun,nodeboun,
              nelemload,sideloadtemp,&nload,&ne,&nk));
      SFREE(sideloadtemp);
  } // end ithermal[1] >=2

  /* storing the undecascaded MPC's */
  if(mpcfreeref==-1)
  {
    printf("Storing uncascaded MPCs...");
    memmpcref_=memmpc_;mpcfreeref=mpcfree;maxlenmpcref=maxlenmpc;
    NNEW(nodempcref,ITG,3*memmpc_);memcpy(nodempcref,nodempc,sizeof(ITG)*3*memmpc_);
    NNEW(coefmpcref,double,memmpc_);memcpy(coefmpcref,coefmpc,sizeof(double)*memmpc_);
    NNEW(ikmpcref,ITG,nmpc);memcpy(ikmpcref,ikmpc,sizeof(ITG)*nmpc);
    printf("done! \n");
  }

  /* decascading MPC's only necessary if MPC's changed */
  /* Necesary for all analysis types */
  if(((istep == 1)||(ntrans>0)||(mpcend<0)||(nk!=nkold)||(nmpc!=nmpcold))&&(icascade==0)) 
  {
    //  if(icascade==0) {

    /* decascading the MPC's */
    printf("Decascading the MPC's...");

    callfrommain=1;

    cascade(ipompc,&coefmpc,&nodempc,&nmpc,
	    &mpcfree,nodeboun,ndirboun,&nboun,ikmpc,
	    ilmpc,ikboun,ilboun,&mpcend,
	    labmpc,&nk,&memmpc_,&icascade,&maxlenmpc,
            &callfrommain,iperturb,ithermal);

    printf("done \n");
  }

  /* determining the matrix structure: changes if SPC's have changed */

  if((icascade==0)&&(nmethod<8))

  NNEW(nactdof,ITG,mt*nk);
  NNEW(mast1,ITG,nzs[1]);
  NNEW(irow,ITG,1);

  if((mcs==0)||(cs[1]<0))
  {

    NNEW(icol,ITG,mt*nk);
    NNEW(jq,ITG,mt*nk+1);
    NNEW(ipointer,ITG,mt*nk);

    if((icascade==0)&&((nmethod<8)||(nmethod>10)))
    {
	    if(nmethod==11)
      {
        nmethodl=2;
      }
      else
      {
        nmethodl=nmethod;
      }

      printf("Analyzing matrix structure...done!\n");
	    mastruct(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,ipompc,
		   nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,neq,
		   ikmpc,ilmpc,ipointer,nzs,&nmethodl,ithermal,
                   ikboun,ilboun,iperturb,mi,&mortar,typeboun,labmpc,
		   &iit,&icascade,&network);
    }

    else
    {
      neq[0]=1;
      neq[1]=1;
      neq[2]=1;
    }
  }
  /* following not active for static aero-elasticity */
  else /* Not used in CalTop */
  {
    NNEW(icol,ITG,8*nk);
    NNEW(jq,ITG,8*nk+1);
    NNEW(ipointer,ITG,8*nk);

    mastructcs(&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,&nboun,
		 ipompc,nodempc,&nmpc,nactdof,icol,jq,&mast1,&irow,&isolver,
		 neq,ikmpc,ilmpc,ipointer,nzs,&nmethod,
		 ics,cs,labmpc,&mcs,mi,&mortar);
  }

  SFREE(ipointer);
  SFREE(mast1);

  if((icascade==0)&&(nmethod<8))
  {
    RENEW(irow,ITG,nzs[2]);
  }
  /* nmethod=1: static analysis   */
  /* nmethod=2: frequency analysis  */
  /* nmethod=3: buckling analysis */
  /* nmethod=4: (linear or nonlinear) dynamic analysis */
  /* nmethod=5: steady state dynamics analysis */
  /* nmethod=6: Coriolis frequency calculation */
  /* nmethod=7: flutter frequency calculation */
  /* nmethod=8:  magnetostatics */
  /* nmethod=9:  magnetodynamics */
  /* nmethod=10: electromagnetic eigenvalue problems */
  /* nmethod=11: superelement creation or Green function calculation */
  /* nmethod=12: sensitivity analysis  */


  /* if solving static analysis or green function calculation */
  if((nmethod<=1)||(nmethod==11)||((iperturb[0]>1)&&(nmethod<8)))
  {
	  if(iperturb[0]<2)
    {
	    mpcinfo[0]=memmpc_;
      mpcinfo[1]=mpcfree;
      mpcinfo[2]=icascade;
	    mpcinfo[3]=maxlenmpc;

	    if(icascade!=0)
      {
	      printf(" *ERROR in CalculiX: the matrix structure may");
	      printf("        change due to nonlinear equations;");
	      printf("        a purely linear calculation is not");
	      printf("        feasible; use NLGEOM on the *STEP card.");
	      FORTRAN(stop,());
	    }

      printf("\nTOPOLOGY OPTIMIZATION PARAMTERS----------------------------|\n");
      printf("SIMP \n");
      printf("  Penalty parameter                 %.2f\n", pstiff);
      printf("\n");
      printf("FILTER MATRIX \n");
      printf("  Density filter radius             %.2f\n", rmin);
      printf("  Non zeros in Filtermatrix         %d\n", fnnzassumed);
      printf("\n");
      printf("STRESS AGGREGATION \n");
      printf("  P-norm exponent                   %.2f\n", pexp);
      printf("  Stress relaxation                 %.2f\n", eps_relax);
      printf("  Stress threshold                  %.2f\n", sigma0);
    
      printf("|------------------------------------------------------------|\n");
     
      fflush(stdout);

      printf("\n");
      if(pSupplied!=0)
      {
        printf("Penalization parameter found --> will evaluate senitivities\n");
      
      }

      if(pSupplied!=0)
      {
        printf("Checking if filter matrix needs to be built \n");


        //NNEW(FilterMatrixs,double,fnnzassumed*ne_); //Sparse filter matrix stored as row,colum,value with fassumed nnzs per element assumed
    
        //NNEW(rowFilters,ITG,fnnzassumed*ne_);
    
        //NNEW(colFilters,ITG,fnnzassumed*ne_);
    
        NNEW(filternnzElems,ITG,ne_);
        NNEW(designFiltered,double,ne_);

        /* Create or assemble the density filter */
        //densityfilter(co,&nk,&kon,&ipkon,&lakon,&ne,&ttime,timepar,&mortar,
        //          &rmin,&filternnz,
        //          FilterMatrixs,rowFilters,colFilters,filternnzElems,itertop,&fnnzassumed);

        //printf("Read within densityfilter:%d \n", filternnz);

        
        /* apply the filter matrix on rho to get rhoPhys */ 
        //filterVector(&ipkon,design,designFiltered,FilterMatrixs,filternnzElems,rowFilters,colFilters,&ne,&ttime,timepar,&fnnzassumed, &qfilter, filternnz);

        /* Try the multi-threaded */
        printf("Filtering element densities...\n");
        //filterDensity_buffered_dat_mt(design, designFiltered, filternnzElems, &ne, &fnnzassumed, &qfilter, filternnz);
        filterDensity_buffered_bin_mt(design, designFiltered, filternnzElems, &ne, &fnnzassumed, &qfilter, filternnz);
        printf("Done!");

        // DEBUG: Print first five and last five values of designFiltered
        /*
        printf("First five designFiltered values:\n");
        for (int i = 0; i < 5 && i < ne_; ++i) 
        {
            printf("  designFiltered[%d] = %g\n", i, designFiltered[i]);
        }

        printf("Last five designFiltered values:\n");
        for (int i = (ne_ > 5 ? ne_ - 5 : 0); i < ne_; ++i) 
        {
            printf("  designFiltered[%d] = %g\n", i, designFiltered[i]);
        }

        */


        rhoPhys=designFiltered;
      }
      else
      {
        printf("No penalization parameter found, initializing all densities to one \n");
        /* design was initialized to 1.0 in rho.c */
        rhoPhys=design;

        // DEBUG only

        /*printf("First five rhoPhys values:\n");

        for (int i = 0; i < 5 && i < ne_; ++i) {
            printf("  rhoPhys[%d] = %g\n", i, rhoPhys[i]);
        }

        printf("Last five rhoPhys values:\n");
        for (int i = (ne_ > 5 ? ne_ - 5 : 0); i < ne_; ++i) 
        {
            printf(" rhoPhys[%d] = %g\n", i, rhoPhys[i]);
        }
        */

      }

      
      //printf("\n For stiffness, penalty considered= %f \n",pstiff);

      time_t startl, endl; 
	    startl = time(NULL);

      printf("\n\nSOLVING LINEAR SYSTEM----------------------------------------|\n\n");
      
  
	    linstatic(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun,
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload,
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,ikmpc,
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,&ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,timepar,thicke,jobnamec,tieset,&ntie,&istep,&nmat,
	     ielprop,prop,typeboun,&mortar,mpcinfo,tietol,ics,&icontact,
	     orname,rhoPhys,&pstiff, stx, &sigma0, &eps_relax, &rhomin, &pexp, &Pnorm, dPnorm_drho);

      endl = time(NULL);

	    printf("\n Time taken for linstatic.c is %.8f seconds \n", 
		  difftime(endl, startl)); 
      

      // NOTE: FILTER, WRITE AND FREE STRESS SENS to reduce memory signature downstrewam
      /*--------------------------------------STRESS SENSITIVITY FILTERING AND I/O -----------------------------------*/
      printf(" Filter element stress (P-norm) gradient ");
      /* Allocate memory for P-norm stress sensitivities */
      // NOTE: P-norm sensitivity initialized and passed to linstatic() earlier
      NNEW(dPnorm_drhoFiltered, double, ne_);
      filterSensitivity_bin_buffered_mts(dPnorm_drho, dPnorm_drhoFiltered, ne, filternnz);

      int rs = write_Stress_sens("stress_sens.csv", ne, dPnorm_drhoFiltered);
      if (rs != 0) 
      {
        printf(" Unable to write P-norm sensitivities to disk!\n");
      }
      SFREE(dPnorm_drho);
      SFREE(dPnorm_drhoFiltered);
      printf("done \n");
      printf("|------------------------------------------------------------|\n");

	    for(i=0;i<3;i++)
      {
        nzsprevstep[i]=nzs[i];
      }

	    memmpc_=mpcinfo[0];
      mpcfree=mpcinfo[1];
      icascade=mpcinfo[2];
      maxlenmpc=mpcinfo[3];

    } // end if(iperturb[0]<2)

    else /* non-linear analysis from here */
    {
	    mpcinfo[0]=memmpc_;
      mpcinfo[1]=mpcfree;
      mpcinfo[2]=icascade;
	    mpcinfo[3]=maxlenmpc;

	    nonlingeo(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun,
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc,&nelemload,&sideload,xload,&nload,
	     nactdof,&icol,jq,&irow,neq,&nzl,&nmethod,&ikmpc,
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,timepar,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,&iamload,iamt1,&alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,tieset,&itpamp,&iviewfile,jobnamec,tietol,&nslavs,thicke,
	     ics,&nintpoint,&mortar,
	     &ifacecount,typeboun,&islavsurf,&pslavsurf,&clearini,&nmat,
	     xmodal,&iaxial,&inext,&nprop,&network,orname,vel,&nef,
	     velo,veloo, rhoPhys,&pstiff);

	      memmpc_=mpcinfo[0];
        mpcfree=mpcinfo[1];
        icascade=mpcinfo[2];
        maxlenmpc=mpcinfo[3];

	    for(i=0;i<3;i++)
      {
        nzsprevstep[i]=nzs[i];
      }
    } // end else
  } //end if((nmethod<=1)||(nmethod==11)||((iperturb[0]>1)&&(nmethod<8)))
  
  
  else if(nmethod==2)
  {
    if((mcs==0)||(cs[1]<0))
    {
      #ifdef ARPACK
	    mpcinfo[0]=memmpc_;
      mpcinfo[1]=mpcfree;
      mpcinfo[2]=icascade;
	    mpcinfo[3]=maxlenmpc;

	    arpack(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,xboun,&nboun,
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload,
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc,
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     shcon,nshcon,cocon,ncocon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,vold,iperturb,sti,nzs,
	     &kode,mei,fei,filab,
	     &iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,&ncmat_,&nstate_,&ener,jobnamec,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,&nslavs,tietol,&nkon,mpcinfo,
	     &ntie,&istep,&mcs,ics,tieset,cs,&nintpoint,&mortar,&ifacecount,
	     &islavsurf,&pslavsurf,&clearini,&nmat,typeboun,ielprop,prop,
             orname,rhoPhys,&pstiff);

	    memmpc_=mpcinfo[0];
      mpcfree=mpcinfo[1];
      icascade=mpcinfo[2];
	    maxlenmpc=mpcinfo[3];

	    for(i=0;i<3;i++)
      {
        nzsprevstep[i]=nzs[i];
      }

      #else
	    printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	    FORTRAN(stop,());
      #endif

    }
    
    else
    {
      #ifdef ARPACK

	    mpcinfo[0]=memmpc_;
      mpcinfo[1]=mpcfree;
      mpcinfo[2]=icascade;
	    mpcinfo[3]=maxlenmpc;

	    arpackcs(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,
             xboun,&nboun,
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc, nelemload,sideload,xload,&nload,
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc,
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     &iexpl,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,matname,mi,ics,cs,&mpcend,&ncmat_,
             &nstate_,&mcs,&nkon,jobnamec,output,set,&nset,istartset,
             iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
             ibody,xbody,&nbody,&nevtot,thicke,&nslavs,tietol,mpcinfo,
	     &ntie,&istep,tieset,&nintpoint,&mortar,&ifacecount,&islavsurf,
	     &pslavsurf,&clearini,&nmat,typeboun,ielprop,prop,orname,rhoPhys,
       &pstiff);

	    memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
	    maxlenmpc=mpcinfo[3];

	    for(i=0;i<3;i++)
      {
        nzsprevstep[i]=nzs[i];
      }
      #else
	      printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
	      FORTRAN(stop,());
      #endif
    }
  } // end if nmethod ==2
  /*
  else if(nmethod==3)
  {

    #ifdef ARPACK
	  arpackbu(co,&nk,kon,ipkon,lakon,&ne,nodeboun,ndirboun,xboun,&nboun,
	     ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc,
	     nelemload,sideload,xload,&nload,
	     nactdof,icol,jq,irow,neq,&nzl,&nmethod,ikmpc,
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
             alcon,nalcon,alzero,ielmat,ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,
	     vold,iperturb,sti,nzs,&kode,mei,fei,filab,
	     eme,&iexpl,plicon,nplicon,plkcon,nplkcon,
	     xstate,&npmat_,matname,mi,&ncmat_,&nstate_,ener,output,
             set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,&isolver,trab,inotr,&ntrans,&ttime,fmpc,cbody,
	     ibody,xbody,&nbody,thicke,jobnamec,&nmat,ielprop,prop,
	     orname,typeboun,rhoPhys,&pstiff);
    #else
      printf("*ERROR in CalculiX: the ARPACK library is not linked\n\n");
      FORTRAN(stop,());
    #endif
  } */
  /*
  else if(nmethod==4)
  {
	  if((ne1d!=0)||(ne2d!=0))
    {
	    printf(" *WARNING: 1-D or 2-D elements may cause problems in modal dynamic calculations\n");
	    printf("           ensure that point loads defined in a *MODAL DYNAMIC step\n");
	    printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	    printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	    printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");}

      printf(" Composing the dynamic response from the eigenmodes\n\n");

      dyna(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,&sti,nzs,
	    timepar,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	    &ntie,tieset,&idrct,jmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,&nslavs,&nmat,typeboun,ielprop,prop,orname);
    } */
    /*
    else if(nmethod==5)
    {
	    if((ne1d!=0)||(ne2d!=0))
      {
	      printf(" *WARNING: 1-D or 2-D elements may cause problems in steady state calculations\n");
	      printf("           ensure that point loads defined in a *STEADY STATE DYNAMICS step\n");
	      printf("           and applied to nodes belonging to 1-D or 2-D elements have been\n");
	      printf("           applied to the same nodes in the preceding FREQUENCY step with\n");
	      printf("           magnitude zero; look at example shellf.inp for a guideline.\n\n");
      }

      printf(" Composing the steady state response from the eigenmodes\n\n");

      steadystate(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	      &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	      nelemload,sideload,xload,&nload,
	      &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	      &t1,ithermal,prestr,&iprestr,&vold,iperturb,sti,nzs,
	      timepar,xmodal,&veold,amname,amta,
	      namta,&nam,iamforc,iamload,&iamt1,
	      jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
	      xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,ics,cs,&mpcend,
	      ctrl,ikforc,ilforc,thicke,&nmat,typeboun,ielprop,prop,orname,
	      &ndamp,dacon);
    } */
    /*
    else if((nmethod==6)||(nmethod==7))
    {

      printf(" Composing the complex eigenmodes from the real eigenmodes\n\n");

      complexfreq(&co,&nk,&kon,&ipkon,&lakon,&ne,&nodeboun,&ndirboun,&xboun,&nboun,
	    &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,&nforc,
	    nelemload,sideload,xload,&nload,
	    &nactdof,neq,&nzl,icol,irow,&nmethod,&ikmpc,&ilmpc,&ikboun,&ilboun,
            elcon,nelcon,rhcon,nrhcon,cocon,ncocon,
            alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,&t0,
	    &t1,ithermal,prestr,&iprestr,&vold,iperturb,&sti,nzs,
	    timepar,xmodal,&veold,amname,amta,
	    namta,&nam,iamforc,iamload,&iamt1,
	    jout,&kode,filab,&eme,xforcold,xloadold,
            &t1old,&iamboun,&xbounold,&iexpl,plicon,
            nplicon,plkcon,nplkcon,xstate,&npmat_,matname,
            mi,&ncmat_,&nstate_,&ener,jobnamec,&ttime,set,&nset,
            istartset,iendset,&ialset,&nprint,prlab,
            prset,&nener,trab,&inotr,&ntrans,&fmpc,cbody,ibody,xbody,&nbody,
            xbodyold,&istep,&isolver,jq,output,&mcs,&nkon,&mpcend,ics,cs,
	    &ntie,tieset,&idrct,jmax,ctrl,&itpamp,tietol,&nalset,
	    ikforc,ilforc,thicke,jobnamef,mei,&nmat,ielprop,prop,orname,
            typeboun);
    } */
    /*
    else if((nmethod>7)&&(nmethod<12))
    {
	    mpcinfo[0]=memmpc_;
      mpcinfo[1]=mpcfree;
      mpcinfo[2]=icascade;
	    mpcinfo[3]=maxlenmpc;

	    electromagnetics(&co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,
             ndirboun,xboun,&nboun,
	     &ipompc,&nodempc,&coefmpc,&labmpc,&nmpc,nodeforc,ndirforc,xforc,
             &nforc,&nelemload,&sideload,xload,&nload,
	     nactdof,&icol,&jq,&irow,neq,&nzl,&nmethod,&ikmpc,
	     &ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr,
	     &vold,iperturb,sti,nzs,&kode,filab,&idrct,jmax,
	     jout,timepar,eme,xbounold,xforcold,xloadold,
	     veold,accold,amname,amta,namta,
	     &nam,iamforc,&iamload,iamt1,&alpha,
             &iexpl,iamboun,plicon,nplicon,plkcon,nplkcon,
	     &xstate,&npmat_,&istep,&ttime,matname,qaold,mi,
	     &isolver,&ncmat_,&nstate_,&iumat,cs,&mcs,&nkon,&ener,
	     mpcinfo,output,
             shcon,nshcon,cocon,ncocon,physcon,&nflow,ctrl,
             &set,&nset,&istartset,&iendset,&ialset,&nprint,prlab,
             prset,&nener,ikforc,ilforc,trab,inotr,&ntrans,&fmpc,
             cbody,ibody,xbody,&nbody,xbodyold,ielprop,prop,
	     &ntie,&tieset,&itpamp,&iviewfile,jobnamec,&tietol,&nslavs,thicke,
	     ics,&nalset,&nmpc_,&nmat,typeboun,&iaxial,&nload_,&nprop,
	     &network,orname,rhoPhys,&pstiff);

	    memmpc_=mpcinfo[0];
      mpcfree=mpcinfo[1];
      icascade=mpcinfo[2];
      maxlenmpc=mpcinfo[3];
    } */

    /* Write elastic fields to a vtu file */
    printf("\nWriting output fields...");
    tecplot_vtu(nk, ne, co, kon, ipkon, vold, stx, rhoPhys);
    printf("done!\n\n");

  
    /* adjoint sensitivity calculation */
    if(pSupplied!=0)
    {
      printf("SENSITIVITY ANALYSIS-----------------------------------------|\n\n");
      
      printf(" Allocating memory...");
      /* allocate memory for compliance gradient and initialize to zero */
      NNEW(gradCompl,double,ne_);

      /* allocate memory for element complaince and initialize to zero */
      NNEW(elCompl,double,ne_);

      /* allocate memory for element C.G and initialize to zero */
      NNEW(elCG,double,3*ne_);

      /* allocate memory for element volume and initialize to zero */
      NNEW(eleVol,double,ne_); 

      /* allocate memory for filtered compliance gradient and initialize to zero */
      NNEW(gradComplFiltered,double,ne_);  //allocate memory to gradcompliance, initialize to 0

      /* allocate memory for filtered volume gradient and initialize to zero */
      NNEW(eleVolFiltered,double,ne_);

      /* Allocate memory for CG sensitivities */
      double *dCGx = (double*)calloc(ne, sizeof(double));
      double *dCGy = (double*)calloc(ne, sizeof(double));
      double *dCGz = (double*)calloc(ne, sizeof(double));

      /* Allocate memory for filteredCG sensitivities */
      double *dCGxFiltered = (double*)calloc(ne, sizeof(double));
      double *dCGyFiltered = (double*)calloc(ne, sizeof(double));
      double *dCGzFiltered = (double*)calloc(ne, sizeof(double));



      printf("done! \n");

      time_t starts, ends; 
	    starts = time(NULL);


      /* Evaluate sensitivities */
      printf( "Evaluating compliance sensitivities...");
	    sensitivity(co,&nk,&kon,&ipkon,&lakon,&ne,nodeboun,ndirboun,
	     xboun,&nboun, ipompc,nodempc,coefmpc,labmpc,&nmpc,nodeforc,
             ndirforc,xforc,&nforc, nelemload,sideload,xload,&nload,
	     nactdof,icol,jq,&irow,neq,&nzl,&nmethod,ikmpc,
	     ilmpc,ikboun,ilboun,elcon,nelcon,rhcon,nrhcon,
	     alcon,nalcon,alzero,&ielmat,&ielorien,&norien,orab,&ntmat_,
             t0,t1,t1old,ithermal,prestr,&iprestr, vold,iperturb,sti,nzs,
	     &kode,filab,eme,&iexpl,plicon,
             nplicon,plkcon,nplkcon,&xstate,&npmat_,matname,
	     &isolver,mi,&ncmat_,&nstate_,cs,&mcs,&nkon,&ener,
             xbounold,xforcold,xloadold,amname,amta,namta,
             &nam,iamforc,iamload,iamt1,iamboun,&ttime,
             output,set,&nset,istartset,iendset,ialset,&nprint,prlab,
             prset,&nener,trab,inotr,&ntrans,fmpc,cbody,ibody,xbody,&nbody,
	     xbodyold,timepar,thicke,jobnamec,tieset,&ntie,&istep,&nmat,
	     ielprop,prop,typeboun,&mortar,mpcinfo,tietol,ics,&icontact,
	     &nobject,&objectset,&istat,orname,nzsprevstep,&nlabel,physcon,
             jobnamef,rhoPhys,&pstiff,gradCompl,elCompl,elCG,eleVol);

      //printf("done!\n");
      // Mass and C.G properties
      double M, cgx, cgy, cgz;
      
      /*---------------------------------C.G SENSITIVITY FILTERING AND I/O ----------------------------------------*/


      printf(" Evaluate and filter CG sensitivities...");
      compute_mass_cg_and_cg_sens(ne, eleVol, rhoPhys, elCG,
                            &M, &cgx, &cgy, &cgz,
                            dCGx, dCGy, dCGz);

      filterSensitivity_bin_buffered_mts3(dCGx, dCGy, dCGz, dCGxFiltered, dCGyFiltered, dCGzFiltered,ne, filternnz);
      printf("done \n");

      printf(" Writing CG sensitivities to disk...");
      /* ... after you fill dCGx, dCGy, dCGz ... */
      int rc = write_cg_sens("cg_sens.csv", ne, dCGxFiltered, dCGyFiltered, dCGzFiltered);
      if (rc != 0) 
      {
        printf("Unable to write CG sensitivities to disk!\n");
      }

      printf("done!\n");

      free(dCGx);
      free(dCGy);
      free(dCGz);
      free(dCGxFiltered);
      free(dCGyFiltered);
      free(dCGzFiltered);
      
      dCGx = NULL;
      dCGy = NULL; 
      dCGz = NULL; 
      dCGxFiltered = NULL;
      dCGyFiltered = NULL;
      dCGzFiltered = NULL;

      /*---------------------------------------------------------------------------------------------------------------*/      


      /*---------------------------------------------------------------------------------------------------------------*/


      /*-------------------------------------COMPLIANCE SENSITIVITY FILTERING AND I/O----------------------------------*/
      double compliance_sum=0;

      printf(" Filter compliance gradient ");
      filterSensitivity_bin_buffered_mts(gradCompl, gradComplFiltered, ne, filternnz);
      printf("done! \n");

      FILE *gradC;
      printf(" Writing compliance sensitivities...");
      write_compliance_sensitivities(ne,gradCompl,gradComplFiltered,elCompl,&compliance_sum);
      printf("done!\n");

      SFREE(gradCompl);
      SFREE(gradComplFiltered);
      
      /*---------------------------------------------------------------------------------------------------------------*/

      /*-------------------------------------VOLUME SENSITIVITY FILTERING AND I/O----------------------------------*/

      FILE *elV_file;
      printf(" Filter element volume gradient ");
      filterSensitivity_bin_buffered_mts(eleVol, eleVolFiltered, ne, filternnz);
      printf("done! \n");
      
      printf(" Writing volume sensitivities...");
      write_volume_sensitivities(ne, eleVol, rhoPhys, eleVolFiltered);
      printf("done!\n");

      SFREE(eleVol);
      SFREE(eleVolFiltered);

      ends = time(NULL);
	    //printf("Time taken for sensitivity calculation: %.2f seconds \n", 
		  //difftime(ends, starts)); 

      /* write compliance gradient to file */
     
      //FILE *elC_file;
      


      /* write compliance value */
      //elC_file=fopen("objectives.dat","w");

      /* set the filtered element densities of passive elements to 0 */
      //filterOutPassiveElems_density(rhoPhys, ne, passiveIDs, numPassive);

      /* set the filtered compliance sens of passive elements to 0 */
      //filterOutPassiveElems_sens(gradComplFiltered, ne, passiveIDs, numPassive);

      /* initialize for compliance */
      
      printf("\n\nOUTPUT FILEDS--------------------------------------------------------------|\n\n");
     
      /* set the filtered volume sens of passive elements to 0 */
      //filterOutPassiveElems_sens(eleVolFiltered, ne, passiveIDs, numPassive);

      
      printf(" Writing objectives...");
      write_objectives(ne, eleVol, rhoPhys, &compliance_sum, &M, &cgx, &cgy, &cgz, &Pnorm);
      printf("done!\n");

      /* initialize for total materal volume with rho = 1 */
      //double initialVol_sum=0;

      /* initialize for total material volume with optimized rho */
      //double designVol_sum=0;

      /* loop over all elements to compute summed values */
      //for (int iii=0;iii<ne;iii++)
      //{
        /* compute initial volume */
      //  initialVol_sum+=eleVol[iii];

        /* compute current design volume */
      //  designVol_sum+=(eleVol[iii]*rhoPhys[iii]);
                
      //}

      /* write summed compliance and volume fraction values to file */
      //fprintf(elC_file,"%.15f , %.15f , %.15f , %.15f \n",compliance_sum,designVol_sum-volfrac*initialVol_sum, initialVol_sum,designVol_sum);

      /* ensure any buffered data is written to file */
      //fflush(elC_file); 

      /* close all files */
      //fclose(elC_file);
   
      /* evaluate discreteness of the structure*/
      /* discreteness is a metric often used in topology optimization problems to assess hpw close the design
       is to a "0-1" or binary solution. */

      //double mnd = 0.0;
        


      /* loop over all element density values */
      //for (int iii=0;iii<ne;iii++)
      //{
      //  mnd+=(rhoPhys[iii]*(1-rhoPhys[iii]));               
     // }

     // Evaluate the discretness of the design


      /* Normalize and average */
      //mnd=(4*mnd*100/ne);
  
     /* print output */
      
      printf("\n Compliance:                 %.3f \n",compliance_sum);
      printf(" Mass:                       %.3f \n", M);
      printf(" Aggregated stress (P-norm): %.3f \n", Pnorm);
      //printf("Total domain volume:         %.6f \n",initialVol_sum);
      //printf("Current domain volume:       %.6f \n",designVol_sum);
      //printf("Volume constraint violation:: %.6f \n",designVol_sum-volfrac*initialVol_sum);
      //printf("Discreteness, mnd, percent:               %.6f \n",mnd);

    } // end adjoint calculation

    printf("\n Writing rhos.dat...");
    fflush(stdout);

    FILE *rho_file;
    rho_file=fopen("rhos.dat","w"); //open in write mode
    
    /* loop over all elements and write element densities to file */
    for (int iii=0;iii<ne;iii++)
    {
      fprintf(rho_file,"%.15f  ,  %.15f \n",design[iii],rhoPhys[iii]);            
    }

    /* ensure any buffered data is written to file */
    fflush(rho_file); 

    /* all operations done, close file */
    fclose(rho_file);

    printf("done!\n");


    SFREE(nactdof);
    SFREE(icol);
    SFREE(jq);
    SFREE(irow);

    /* deleting the perturbation loads and temperatures */

    if((iperturb[0] == 1)&&(nmethod==3)) 
    {
      nforc=0;
      nload=0;
      nbody=0;
      if(ithermal[0] == 1) 
      {
	      for(k=0;k<nk;++k)
        {
	        t1[k]=t0[k];
	      }
      }
    }
    else
    {
      nbounold=nboun;

      for (i=0;i<nboun;i++) 
      {
	      nodebounold[i]=nodeboun[i];
	      ndirbounold[i]=ndirboun[i];
      }
      nforcold=nforc;
      nloadold=nload;
      nbodyold=nbody;

      /* resetting the amplitude to none except for time=total time amplitudes */

      if(nam > 0) 
      {
	      for (i=0;i<nboun;i++) 
        {
	        if(iamboun[i]>0)
          {
		        if(namta[3*iamboun[i]-1]>0)
            {
		          iamboun[i]=0;
		          xboun[i]=xbounold[i];
            }
	        }
	      }

	      for (i=0;i<nforc;i++)
        {
	        if(iamforc[i]>0)
          {
		        if(namta[3*iamforc[i]-1]>0)
            {
		          iamforc[i]=0;
		          xforc[i]=xforcold[i];
            }
	        }
	      }

	      for (i=0;i<2*nload;i++)
        {
	        if(iamload[i]>0)
          {
		        if(namta[3*iamload[i]-1]>0)
            {
		          iamload[i]=0;
		          xload[i]=xloadold[i];
            }
	        }
	      }

	      for (i=1;i<3*nbody;i=i+3)
        {
	        if(ibody[i]>0)
          {
		        if(namta[3*ibody[i]-1]>0)
            {
		          ibody[i]=0;
		          xbody[7*(i-1)/3]=xbodyold[7*(i-1)/3];
            }
	        }
	      }
	      if(ithermal[0]==1) 
        {
	        if(iamt1[i]>0)
          {
		        if(namta[3*iamt1[i]-1]>0)
            {
		          iamt1[i]=0;
		          t1[i]=t1old[i];
            }
	        }
	      }
      }
    }

    /* removing the advective elements, if any */
    if(network>0)
    {
      ne=ne0;nkon=nkon0;
      RENEW(ipkon,ITG,ne);
      RENEW(lakon,char,8*ne);
      RENEW(kon,ITG,nkon);
      RENEW(sti,double,6*mi[0]*ne);
      RENEW(eme,double,6*mi[0]*ne);

      if(iprestr>0) RENEW(prestr,double,6*mi[0]*ne);
      if(nprop>0) RENEW(ielprop,ITG,ne);
      if((ne1d!=0)||(ne2d!=0)) RENEW(offset,double,2*ne);
      if(nener==1)RENEW(ener,double,mi[0]*ne*2);
      if(norien>0)RENEW(ielorien,ITG,mi[2]*ne);
      RENEW(ielmat,ITG,mi[2]*ne);

      /* reactivating the original load labels */

      for(i=nload-1;i>=nload0;i--)
      {
	      if(strcmp2(&sideload[20*i],"                    ",20)==0)
        {
	        iload=nelemload[2*i+1];
	        strcpy1(&sideload[20*(iload-1)],"F",1);
	      }
      }
    }

  nload=nload0;

  if((nmethod == 4)&&(iperturb[0]>1))
  {
    SFREE(accold);
  }
  if(irstrt[0]>0)
  {
    jrstrt++;

    if(jrstrt>=irstrt[0])
    {
      jrstrt=0;
      FORTRAN(restartwrite,(&istep,&nset,&nload,&nforc,&nboun,&nk,&ne,
        &nmpc,&nalset,&nmat,&ntmat_,&npmat_,&norien,&nam,&nprint,
        mi,&ntrans,&ncs_,&namtot,&ncmat_,&mpcend,&maxlenmpc,&ne1d,
        &ne2d,&nflow,&nlabel,&iplas,&nkon,ithermal,&nmethod,iperturb,
        &nstate_,&nener,set,istartset,iendset,ialset,co,kon,ipkon,
        lakon,nodeboun,ndirboun,iamboun,xboun,ikboun,ilboun,ipompc,
        nodempc,coefmpc,labmpc,ikmpc,ilmpc,nodeforc,ndirforc,iamforc,
        xforc,ikforc,ilforc,nelemload,iamload,sideload,xload,
        elcon,nelcon,rhcon,nrhcon,alcon,nalcon,
        alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,ielorien,
        trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
        ielmat,matname,prlab,prset,filab,vold,nodebounold,
        ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
        iponor,xnor,knor,thicke,offset,iponoel,inoel,rig,
        shcon,nshcon,cocon,ncocon,ics,
	      sti,ener,xstate,jobnamec,infree,prestr,&iprestr,cbody,
	      ibody,xbody,&nbody,xbodyold,&ttime,qaold,cs,&mcs,output,
	      physcon,ctrl,typeboun,fmpc,tieset,&ntie,tietol,&nslavs,t0g,t1g,
	      &nprop,ielprop,prop,&mortar,&nintpoint,&ifacecount,islavsurf,
	      pslavsurf,clearini,irstrt,vel,&nef,velo,veloo));
    }
  }
  if(istep == 1)
  {
    printf("Linear analysis complete!\n");
    break;
  }

 } // end while(istat>=0)

  printf("De-allocate memory...");

  FORTRAN(closefile,());

  /*
  strcpy(fneig,jobnamec);
  strcat(fneig,".frd");

  printf("I am trying to write the frd file.\n");
  if((f1=fopen(fneig,"ab"))==NULL)
  {
    printf("*ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }

  fprintf(f1," 9999\n");
  fclose(f1);
  printf("FRD file sucessfully closed!\n");
  */
/* deallocating the fields
   this section is addressed immediately after leaving calinput */

  

//free up space
/* if(pSupplied!=0){SFREE(filternnzElem);SFREE(FilterMatrix);SFREE(rowFilter);SFREE(colFilter);
   SFREE(gradCompl);SFREE(elCompl);SFREE(eleVol);SFREE(elCG);
   SFREE(gradComplFiltered);SFREE(eleVolFiltered);SFREE(designFiltered);
 }*/



  /* Free topolology-optimization related fields */
  SFREE(filternnzElems);
  SFREE(FilterMatrixs);
  SFREE(rowFilters);
  SFREE(colFilters);
 
  
  
  SFREE(elCG);
  //SFREE(gradComplFiltered);
  //SFREE(eleVolFiltered);
  SFREE(designFiltered);

  rhoPhys = NULL;
  free(passiveIDs);





  SFREE(ipoinpc);
  SFREE(inpc);
  SFREE(inp);
  SFREE(ipoinp);
  if(ncs_>0) SFREE(ics);
  if(mcs>0) SFREE(cs);
  SFREE(tieset);
  SFREE(tietol);
  SFREE(co);
  SFREE(kon);
  SFREE(ipkon);
  SFREE(lakon);
  SFREE(design);
  SFREE(nodeboun);
  SFREE(ndirboun);
  SFREE(typeboun);
  SFREE(xboun);
  SFREE(ikboun);
  SFREE(ilboun);
  SFREE(nodebounold);
  SFREE(ndirbounold);
  SFREE(xbounold);
  SFREE(ipompc);
  SFREE(labmpc);
  SFREE(ikmpc);
  SFREE(ilmpc);
  SFREE(fmpc);
  SFREE(nodempc);
  SFREE(coefmpc);
  SFREE(nodempcref);
  SFREE(coefmpcref);
  SFREE(ikmpcref);
  SFREE(nodeforc);
  SFREE(ndirforc);
  SFREE(xforc);
  SFREE(ikforc);
  SFREE(ilforc);
  SFREE(xforcold);
  SFREE(nelemload);
  SFREE(sideload);
  SFREE(xload);
  SFREE(xloadold);
  SFREE(cbody);
  SFREE(ibody);
  SFREE(xbody);
  SFREE(xbodyold);
  SFREE(stx);
  if(nam>0)
  {
    SFREE(iamboun);
    SFREE(iamforc);
    SFREE(iamload);
    SFREE(amname);
    SFREE(amta);
    SFREE(namta);
  }
  SFREE(set);
  SFREE(istartset);
  SFREE(iendset);
  SFREE(ialset);
  SFREE(elcon);
  SFREE(nelcon);
  SFREE(rhcon);
  SFREE(nrhcon);
  SFREE(shcon);
  SFREE(nshcon);
  SFREE(cocon);
  SFREE(ncocon);
  SFREE(alcon);
  SFREE(nalcon);
  SFREE(alzero);
  if(nprop>0)
  {
    SFREE(ielprop);SFREE(prop);
  }

  if(npmat_>0)
  {
    SFREE(plicon);
    SFREE(nplicon);
    SFREE(plkcon);
    SFREE(nplkcon);
  }

  if(ndamp>0)
  {
    SFREE(dacon);
  }

  if(norien>0)
  {
    SFREE(orname);
    SFREE(orab);
    SFREE(ielorien);
  }

  if(ntrans>0)
  {
    SFREE(trab);
    SFREE(inotr);
  }

  if(iprestr>0)
  {
    SFREE(prestr);
  }

  if(ithermal[0]!=0)
  {
    SFREE(t0);
    SFREE(t1);
    SFREE(t1old);

    if(nam>0) SFREE(iamt1);

    if((ne1d!=0)||(ne2d!=0))
    {
      SFREE(t0g);
      SFREE(t1g);
    }
  }

  SFREE(prlab);
  SFREE(prset);
  SFREE(filab);
  SFREE(xmodal);
  SFREE(ielmat);
  SFREE(matname);
  SFREE(sti);
  SFREE(eme);
  SFREE(ener);
  SFREE(xstate);
  SFREE(vold);
  SFREE(veold);
  SFREE(vel);
  SFREE(velo);
  SFREE(veloo);


  if((ne1d!=0)||(ne2d!=0))
  {
    SFREE(iponor);
    SFREE(xnor);
    SFREE(knor);
    SFREE(thicke);
    SFREE(offset);
    SFREE(iponoel);
    SFREE(inoel);
    SFREE(rig);
  }

  SFREE(islavsurf);

  if(mortar==1)
  {
    SFREE(pslavsurf);
    SFREE(clearini);
  }

  if(nobject_>0)
  {
    SFREE(objectset);
  }

  #ifdef CALCULIX_MPI
  MPI_Finalize();
  #endif

  #ifdef CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT
  calculix_freeExternalBehaviours();
  #endif /* CALCULIX_EXTERNAL_BEHAVIOURS_SUPPORT */
  printf("done! \n");
  return 0;
}
