#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include <time.h>

int main(int argc, char *argv[])
{
    FILE *f1;

    char *sideload = NULL;  /*!< \brief load label.*/
    char *set = NULL;       /*!< \brief magnitude of load at end of a step.*/
    char *matname = NULL;   /*!< \brief name of material.*/
    char *orname = NULL;    /*!< \brief name of orientation.*/
    char *amname = NULL;    /*!< \brief name of the amplitude.*/
    char *filab = NULL;     /*!< \brief label corresponding to outfield.*/
    char *lakon = NULL;     /*!< \brief element label.*/
    char *labmpc = NULL;    /*!< \brief label of MPI i.*/
    char *prlab = NULL;     /*!< \brief label for output request.*/
    char  *prset = NULL;    /*!< \brief node or element set corresponding to output request i.*/
    char jobnamec[660] = "";/*!< \brief jobname for C environment. */
    char jobnamef[132] = "";/*!< \brief jobname for FORTRAN environment. */
    char output[4] = "asc"; /*!< \brief Output filename character array.*/
    char *typeboun = NULL;  /*!< \brief SPC type (A=acceleration, B = prescribed b.c, R= Rigidbody)/*/
    char *inpc = NULL;      /*!< \brief Undefined.*/
    char *tieset = NULL;    /*!< \brief name of the tie constraint.*/
    char *cbody = NULL;     /*!< \brief element number or element set to which load i applies.*/
    char fneig[132] = "";   /*!< \brief Undefined.*/
    char *sideloadtemp = NULL;/*!< \brief load label temp.*/
    char  kind1[2] = "T";   /*!< \brief Undefined.*/
    char kind2[2] = "T";    /*!< \brief Undefined.*/
    char *heading = NULL;   /*!< \brief Undefined.*/
    char *objectset = NULL; /*!< \brief Undefined.*/

      ITG nk, ne, nboun, nmpc, nforc, nload, nprint = 0, nset, nalset, nentries = 17,
                                         nmethod, neq[3] = {0, 0, 0}, i, mpcfree = 1, mei[4], j, nzl, nam, nbounold = 0,
                                         nforcold = 0, nloadold = 0, nbody, nbody_ = 0, nbodyold = 0, network = 0, nheading_ = 0,
                                         k, nzs[3], nmpc_ = 0, nload_ = 0, nforc_ = 0, istep, istat, nboun_ = 0, nintpoint = 0,
                                         iperturb[2] = {0, 0}, nmat, ntmat_ = 0, norien, ithermal[2] = {0, 0}, nmpcold,
                                         iprestr, kode, isolver = 0, nslavs = 0, nkon_ = 0, ne0, nkon0, mortar = 0,
                                         jout[2] = {1, 1}, nlabel, nkon = 0, idrct, jmax[2], iexpl, nevtot = 0, ifacecount = 0,
                                         iplas = 0, npmat_ = 0, mi[3] = {0, 3, 1}, ntrans, mpcend = -1, namtot_ = 0, iumat = 0,
                                         icascade = 0, maxlenmpc, mpcinfo[4], ne1d = 0, ne2d = 0, infree[4] = {0, 0, 0, 0},
                                         callfrommain, nflow = 0, jin = 0, irstrt[2] = {0, 0}, nener = 0, jrstrt = 0, nenerold,
                                         nline, *ipoinp = NULL, *inp = NULL, ntie, ntie_ = 0, mcs = 0, nprop_ = 0,
                                         nprop = 0, itpamp = 0, iviewfile, nkold, nevdamp_ = 0, npt_ = 0, cyclicsymmetry,
                                         nmethodl, iaxial = 1, inext = 0, icontact = 0, nobject = 0, nobject_ = 0, iit = -1,
                                         nzsprevstep[3], memmpcref_, mpcfreeref = -1, maxlenmpcref, *nodempcref = NULL,
                                         *ikmpcref = NULL, isens = 0, namtot = 0, nstam = 0, ndamp = 0, nef = 0;   

    ITG *kon = NULL;/*!< \brief Field containing connectivity lists of the elements in sucessive order.*/
    ITG *nodeboun = NULL;/*!< \brief SPC node.*/
    ITG *ndirboun = NULL; /*!< \brief SPC direction.*/
    ITG *ipompc = NULL;/*!< \brief starting location in nodempc and coefmpc of MPC i*/
    ITG *nodempc = NULL;/*!< \brief node of first term of MPC i.*/
    ITG *nodeforc = NULL;/*!< \brief node in which force is applied.*/
    ITG *ndirforc = NULL; /*!< \brief direction of force.*/
    ITG *nelemload = NULL;/*!< \brief element to which distributed load is applied.*/
    ITG im;/*!< \brief Undefined.*/
    ITG *inodesd = NULL;/*!< \brief Undefined.*/              
    ITG nload1;/*!< \brief # of facial distributed loads.*/
    ITG *idefforc = NULL;/*!< \brief 0: no force was defined for this node. 1: at least one force was applied.*/
    ITG *nactdof = NULL;/*!< \brief actual degree of freedom (in the system of equations) of DOF i of node j*/
    ITG *icol = NULL;/*!< \brief # of subdiagonal nonzero's in column i(only for symmteric matrices)*/
    ITG *ics = NULL;/*!< \brief one-dimensional field; contains all independent nodes, one part after the other, and stored within each part.*/  
    ITG *jq = NULL; /*!< \brief location in field irow of the first subdiagonal nonzero in column i (only for symmetric matrices)*/
    ITG *mast1 = NULL;/*!< \brief Undefined.*/
    ITG *irow = NULL;/*!< \brief row of element i in field au.*/
    ITG *rig = NULL;/*!< \brief integer field indicating whether node i is a rigid node (non-zero value) or not (zero value)*/
    ITG *idefbody = NULL;/*!< \brief 0: no body load was defined on the same set with the same code and the same load case number before within the actual step.*/
    ITG *ikmpc = NULL;/*!< \brief ordered array of the dependent DOFs corresping to the MPCs*/
    ITG *ilmpc = NULL;/*!< \brief original MPC number for ikmpc(i)*/
    ITG *ikboun = NULL;/*!< \brief ordered array of the DOFs corresponding to the SPC's.*/          
    ITG *ilboun = NULL;/*!< \brief original SPC number for ikboun(i).*/
    ITG *nreorder = NULL;/*!< \brief Undefined.*/
    ITG *ipointer = NULL;/*!< \brief Undefined.*/
    ITG *idefload = NULL;/*!< \brief 0: no load was defned on the same element with the same label and on the same sector 1: at least one load was defined on the same element with the same label.*/
    ITG *istartset = NULL;/*!< \brief pointer into ialset containing the first set member.*/        
    ITG *iendset = NULL; /*!< \brief pointer into ialset containing the last set member. */
    ITG *ialset = NULL;/*!< \brief member of a set or surface: this is a node for a node or nodal surface, element for an element set.*/
    ITG *ielmat = NULL;/*!< \brief material number of layer j.*/
    ITG *ielorien = NULL;/*!< \brief orientation number of layer j. */
    ITG *nrhcon = NULL;/*!< \brief temperature data points for the density of material i.*/
    ITG *nodebounold = NULL;/*!< \brief SPC node old.*/
    ITG *ndirbounold = NULL;/*!< \brief SPC node direction old.*/  
    ITG *nelcon = NULL;/*!< \brief hyperelastic constants for material i (negative kode for nonlinear elastic constants.)*/
    ITG *nalcon = NULL;/*!< \brief # of expansion constants for material i.*/
    ITG *iamforc = NULL;/*!< \brief amplitude number.*/   
    ITG *iamload = NULL;/*!< \brief amplitude number for xload(1,i).*/    
    ITG *iamt1 = NULL;/*!< \brief amplitude number.*/      
    ITG *namta = NULL;/*!< \brief location of first (time, amplitude) pair in amta. */ 
    ITG *ipkon = NULL;/*!< \brief points to the location in field kon preceding the topology of element i.*/ 
    ITG *iamboun = NULL;/*!< \brief amplitude number.*/ 
    ITG *nplicon = NULL;/*!< \brief temperature data points for the isotropic hardening cruve of material i.*/ 
    ITG *nplkcon = NULL; /*!< \brief temperature data points for the kinematic hardening curve of material i.*/ 
    ITG *inotr = NULL;/*!< \brief transformation number applicable in node j.*/ 
    ITG *iponor = NULL;/*!< \brief two pointers for entry i of kon. 1st ptr-> location of xnor precceding normals of entry i. 2nd ptr-> location of kow of newly generated dependent nodes of entry i.*/                                  
    ITG *knor = NULL;/*!< \brief field containing extra nodes needed to expand the shell and beam elements to volume elements.*/ 
    ITG *ikforc = NULL;/*!< \brief ordered array of the DOFs corresponding to the point loads*/ 
    ITG *ilforc = NULL;/*!< \brief original SPC number for ikforc(i)*/ 
    ITG *iponoel = NULL;/*!< \brief ptr to node i into field inoel, which stores the 1D and 2D elements belonging to the node.*/ 
    ITG *inoel = NULL;/*!< \brief field containing an element number, a local node number within this element and a pointer to another array (or zero if there is no other.)*/ 
    ITG *nshcon = NULL;/*!< \brief # temperature data points for the spaecific heat of material i.*/ 
    ITG *ncocon = NULL;/*!< \brief # of conductivity constants for material i.*/ 
    ITG *ibody = NULL;/*!< \brief code identifying the kind of body load: 1: centrifugal 2: gravity, 3: load case.*/ 
    ITG *ielprop = NULL;/*!< \brief ptr to the position in field prop after which the properties for element i start.*/ 
    ITG *islavsurf = NULL;/*!< \brief Undefined.*/
    ITG *ipoinpc = NULL;/*!< \brief index of the first col in field inp containing info on a block of lines in the inputdeck corresponding to a fundamental key.*/     
    ITG mt;/*!< \brief Undefined.*/
    ITG nxstate;/*!< \brief Undefined.*/
    ITG nload0;/*!< \brief # of facial distributed loads. */
    ITG iload;/*!< \brief Undefined.*/
    ITG *iuel = NULL;/*!< \brief type number of the user element, # of integration pts, max fof, # of nodes in any element.*/              

    ITG *meminset = NULL; /*!< \brief Undefined.*/
    ITG *rmeminset = NULL; /*!< \brief Undefined.*/

    ITG nzs_, nk_ = 0, ne_ = 0, nset_ = 0, nalset_ = 0, nmat_ = 0, norien_ = 0, nam_ = 0,
            ntrans_ = 0, ncs_ = 0, nstate_ = 0, ncmat_ = 0, memmpc_ = 0, nprint_ = 0, nuel_ = 0;

    double *co = NULL, *xboun = NULL, *coefmpc = NULL, *xforc = NULL, *clearini = NULL,
         *xload = NULL, *xbounold = NULL, *xforcold = NULL,
         *vold = NULL, *sti = NULL, *xloadold = NULL, *xnor = NULL,
         *reorder = NULL, *dcs = NULL, *thickn = NULL, *thicke = NULL, *offset = NULL,
         *elcon = NULL, *rhcon = NULL, *alcon = NULL, *alzero = NULL, *t0 = NULL, *t1 = NULL,
         *prestr = NULL, *orab = NULL, *amta = NULL, *veold = NULL, *accold = NULL,
         *t1old = NULL, *eme = NULL, *plicon = NULL, *pslavsurf = NULL, *plkcon = NULL,
         *xstate = NULL, *trab = NULL, *ener = NULL, *shcon = NULL, *cocon = NULL,
         *cs = NULL, *tietol = NULL, *fmpc = NULL, *prop = NULL, *t0g = NULL, *t1g = NULL,
         *xbody = NULL, *xbodyold = NULL, *coefmpcref = NULL, *dacon = NULL, *vel = NULL,
         *velo = NULL, *veloo = NULL, *design = NULL, *rhoPhys = NULL;

    double *gradCompl = NULL, *elCompl = NULL, *elCG = NULL, *eleVol = NULL; // gradCompl=compliance gradient of each element,elcompl=compliance of elements
    double *designFiltered = NULL, *gradComplFiltered = NULL, *eleVolFiltered = NULL, *gradVol = NULL;

    double ctrl[56] = {4.5, 8.5, 9.5, 16.5, 10.5, 4.5, 0., 5.5, 0., 0., 0.25, 0.5, 0.75, 0.85, 0., 0., 1.5, 0., 0.005, 0.01, 0., 0., 0.02, 1.e-5, 1.e-3, 1.e-8, 1.e30, 1.5, 0.25, 1.01, 1., 1., 5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, 5.e-7, -1., 1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.e20, 1.5, 0.5, 20.5, 1.5, 1.5, 0.001, 0.1, 100.5, 60.5};

    double fei[3], *xmodal = NULL, timepar[5],
                 alpha, ttime = 0., qaold[2] = {0., 0.}, physcon[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    double pSupplied = 0.0, pstiff = 0.0, rmin = 0.000000000001, volfrac = 1.00, qfilter = 3;

    ITG itertop = 1;       // counter for iteration number in topology optimization
    ITG fnnzassumed = 300; // We assume 500 non zeros in each row of filtermatrix
    // filternnz=total number of nonzeros in filtermatrix,filternnzElem=no of nonzeros in each row of filtermatrix
    ITG filternnz = 0; // For actual nnz values in filter matrix

    ITG *filternnzElems = NULL;


    double *FilterMatrixs = NULL; // Pointer to filter matrix

    ITG *rowFilters = NULL; // Pointer to row index
    ITG *colFilters = NULL; // Pointer to column matrix

    ITG caltop_mode = 3; // 2: write densityfsi.dat and exit


    /* Parse input arguments */
    if (argc == 1)
    {
        printf("Usage: Flags: -i jobname -p PENALTY -q HIGHERORDERFILTER -r RADIUS -v VOLUMEFRACTION -s ITERATION -f FILTERNNZ -m CALTOP_mode\n");
        FORTRAN(stop, ());
    }

    else
    {
        for (i = 1; i < argc; i++)
        {   
            if (strcmp1(argv[i], "-i") == 0)
            {
                strcpy(jobnamec, argv[i + 1]);
                strcpy1(jobnamef, argv[i + 1], 132);
                jin++;
                break;
            }
            if (strcmp1(argv[i], "-v") == 0)
            {
                printf("\nThis is Version 2.15 modified for Topology Optimization with SIMP method: Ghanendra Kumar Das,CDILab, AE, UIUC\n\n");
                printf("Last edit:15 April, 2022 to add passive element sets \n\n");
                FORTRAN(stop, ());
            }
        }
        if (jin == 0)
        {
            strcpy(jobnamec, argv[1]);
            strcpy1(jobnamef, argv[1], 132);
        }

        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-o") == 0)
            {
                strcpy(output, argv[i + 1]);
                break;
            }
        }

        // Read penalty p
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-p") == 0)
            {
                pSupplied = atof(argv[i + 1]);
                break;
            }
        }

        // Read filter radius r
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-r") == 0)
            {
                rmin = atof(argv[i + 1]);
                break;
            }
        }

        // Read filter order q
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-q") == 0)
            {
                qfilter = atof(argv[i + 1]);
                break;
            }
        }

        // Read volume fraction v
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-v") == 0)
            {
                volfrac = atof(argv[i + 1]);
                break;
            }
        }

        // Read expected non zeros in filtermatrix
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-f") == 0)
            {
                fnnzassumed = atoi(argv[i + 1]);
                break;
            }
        }

        // Read iteration number
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-s") == 0)
            {
                itertop = atoi(argv[i + 1]);
                break;
            }
        }

        // Read mode
        for (i = 1; i < argc; i++)
        {
            if (strcmp1(argv[i], "-m") == 0)
            {
                caltop_mode = atoi(argv[i + 1]);
                break;
            }
        }
    }

    if (pSupplied == 0.0)
    {
        pstiff = 1.0;
    }
    else
    {
        pstiff = pSupplied;
    }

    // setenv("CCX_JOBNAME_GETJOBNAME",jobnamec,1);
    putenv("CCX_JOBNAME_GETJOBNAME=jobnamec");

    FORTRAN(openfile, (jobnamef, output));

    printf("\n************************************************************\n\n");
    printf("CalculiX Version 2.15 Topology Optimization with SIMP Method, Ghanendra Kumar Das, CDILab, Aerospace Dept, UIUC \n");
    printf("Original FEA source code:CalculiX by G Dhond.\n");
    printf("************************************************************\n\n");
    printf("\n\n");
    fflush(stdout);

    istep = 0;
    istat = 0;
    iprestr = 0;
    kode = 0;

    #if defined(SGI)
    isolver = 4;
    #elif defined(PARDISO)
    isolver = 7;
    #elif defined(SPOOLES)
    isolver = 0;
    #elif defined(TAUCS)
    isolver = 5;
    #else
    isolver = 3;
    #endif

    NNEW(ipoinp, ITG, 2 * nentries);

}