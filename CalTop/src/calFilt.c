
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

int main(int argc,char *argv[])
{

  FILE *f1;

  char *sideload=NULL;    /**< load label for distributed loads. */
  char *set=NULL;         /**< node set or element set */
  char *matname=NULL;     /**< material name */
  char *orname=NULL;      /**< name of orientation of element */
  char *amname=NULL;      /**< amplitude name. */
  char *filab=NULL;       /**< field label (U, S,..) */
  char *lakon=NULL;       /**< element label */
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

  ITG *kon=NULL;          /**< Element topology vector */
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
  ITG *ipkon=NULL;    /** < element connectivity */
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
	double *vold=NULL; /**< Unknown (?) */
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

  if (argc ==1)
  {
    /* Inadequate input arguments */

    printf("Usage: Flags: -i jobname -r LENGTH SCALE -f FILTERNNZ \n");
    FORTRAN(stop,());
  }

  else
  {
    /* Loop over all input arguments */
    for (i = 1; i < argc; i++)
    {
        if(strcmp1(argv[i],"-i")==0)
        {
            /* Input file found, allocate  forjob name */
            strcpy(jobnamec,argv[i+1]);
            strcpy1(jobnamef,argv[i+1],132);
            jin++;
            break;
        }

        if(jin==0)
        {
            strcpy(jobnamec,argv[1]);strcpy1(jobnamef,argv[1],132);
        }

        /* Read length scale */
        for(i=1; i<argc; i++)
        {
            if(strcmp1(argv[i],"-r")==0)
            {
                rmin=atof(argv[i+1]);
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
  } // DOne reading all input flags

  /* Assign a default value for penalty parameter */
  if (pSupplied == 0.0)
  {
    pstiff = 1.0;
  }
  else
  {
    pstiff = pSupplied;
  }

  putenv("CCX_JOBNAME_GETJOBNAME=jobnamec");

    #ifdef BAM
    ITG lop=0,lrestart=0,kstep=1,kinc=1;
    double time[2],dtime;
    FORTRAN(uexternaldb,(&lop,&lrestart,time,&dtime,&kstep,&kinc));
    #endif

    printf("\n");

    printf("  #####    #####   ##        #######  #####   ######  \n");
    printf(" ##       ##   ##  ##           ##   ##   ##  ##   ## \n");
    printf(" ##       #######  ##           ##   ##   ##  ######  \n");
    printf(" ##       ##   ##  ##           ##   ##   ##  ##      \n");
    printf("  #####   ##   ##  #######      ##    #####   ##      \n");
    
    printf("\n");

    printf("* Contributors:\n");
    printf("* Prateek Ranjan, Dept. of Aerospace Engineering,\n");
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

    /* Conservative estimate of the fielts to be allocated */
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

    nlabel=48;

    /* if solving static analysis problem */
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

            printf("\n#------------------------------DENSITY FILTER PARAMTERS----------------------------#\n");
            printf("Length scale             %.2f\n", rmin);
            printf("Non zeros in Filtermatrix         %d\n", fnnzassumed);
            printf("#-------------------------------------------------------------------------------------#\n");

            fflush(stdout);

            printf("Building filter matrix \n");
            
            /* Sparse filter matrix stored as row,colum,value with fassumed nnzs per element assumed */
            NNEW(FilterMatrixs,double,fnnzassumed*ne_);

            NNEW(rowFilters,ITG,fnnzassumed*ne_);

            NNEW(colFilters,ITG,fnnzassumed*ne_);

            NNEW(filternnzElems,ITG,ne_);

            NNEW(designFiltered,double,ne_);

            /* Build the filter matrix */

            densityfilter(co,&nk,&kon,&ipkon,&lakon,&ne,&ttime,timepar,&mortar,
                  &rmin,&filternnz,
                  FilterMatrixs,rowFilters,colFilters,filternnzElems,itertop,&fnnzassumed);

            printf("Filter matrix files written to disk! \n");


        } 
    }

    // Begin de-allocating memory 
    SFREE(filternnzElems);
    SFREE(FilterMatrixs);
    SFREE(rowFilters);
    SFREE(colFilters);
    SFREE(gradCompl);
    SFREE(elCompl);
    SFREE(eleVol);
    SFREE(elCG);
    SFREE(gradComplFiltered);
    SFREE(eleVolFiltered);
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
    return 0
}  // main end