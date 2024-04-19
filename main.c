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

    ITG nk, ne, nboun, nmpc, nforc, nload, nprint = 0, nset, nalset, 
    nmethod, neq[3] = {0, 0, 0}, i, mpcfree = 1, mei[4], j, nzl, nam, nbounold = 0,
                                         nforcold = 0, nloadold = 0, nbody, nbody_ = 0, nbodyold = 0, network = 0, nheading_ = 0,
                                         k, nzs[3], nmpc_ = 0, nload_ = 0, nforc_ = 0, istep, istat, nboun_ = 0, nintpoint = 0,
                                         iperturb[2] = {0, 0}, nmat, ntmat_ = 0, norien, ithermal[2] = {0, 0}, nmpcold,
                                         iprestr, kode, isolver = 0, nslavs = 0, nkon_ = 0, ne0, nkon0, mortar = 0,
                                         jout[2] = {1, 1}, nlabel, nkon = 0, idrct, jmax[2], iexpl, nevtot = 0, ifacecount = 0,
                                         iplas = 0, npmat_ = 0, mi[3] = {0, 3, 1}, ntrans, mpcend = -1, namtot_ = 0, iumat = 0,
                                         icascade = 0, maxlenmpc, mpcinfo[4], ne1d = 0, ne2d = 0, infree[4] = {0, 0, 0, 0},
                                         callfrommain, nflow = 0, jin = 0, irstrt[2] = {0, 0}, nener = 0, jrstrt = 0, nenerold,
                                         nline, *inp = NULL, ntie, ntie_ = 0, mcs = 0, nprop_ = 0,
                                         nprop = 0, itpamp = 0, iviewfile, nkold, nevdamp_ = 0, npt_ = 0, cyclicsymmetry,
                                         nmethodl, iaxial = 1, inext = 0, icontact = 0, nobject = 0, nobject_ = 0, iit = -1,
                                         nzsprevstep[3], memmpcref_, mpcfreeref = -1, maxlenmpcref, *nodempcref = NULL,
                                         *ikmpcref = NULL, isens = 0, namtot = 0, nstam = 0, ndamp = 0, nef = 0;   

    ITG *ipoinp = NULL;
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

    ITG nentries = 17;      

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
    printf("Calling readinput.c \n");
    /* conservative estimate of the fields to be allocated */
    readinput(jobnamec, &inpc, &nline, &nset_, ipoinp, &inp, &ipoinpc, ithermal, &nuel_);

    printf("Done! \n");

    SFREE(ipoinp);

    NNEW(set, char, 81 * nset_);
    NNEW(meminset, ITG, nset_);
    NNEW(rmeminset, ITG, nset_);
    NNEW(iuel, ITG, 4 * nuel_);

    FORTRAN(allocation, (&nload_, &nforc_, &nboun_, &nk_, &ne_, &nmpc_, &nset_, &nalset_,
                       &nmat_, &ntmat_, &npmat_, &norien_, &nam_, &nprint_, mi, &ntrans_,
                       set, meminset, rmeminset, &ncs_, &namtot_, &ncmat_, &memmpc_, &ne1d,
                       &ne2d, &nflow, jobnamec, irstrt, ithermal, &nener, &nstate_, &istep,
                       inpc, ipoinp, inp, &ntie_, &nbody_, &nprop_, ipoinpc, &nevdamp_, &npt_, &nslavs,
                       &nkon_, &mcs, &mortar, &ifacecount, &nintpoint, infree, &nheading_, &nobject_,
                       iuel, &iprestr, &nstam, &ndamp, &nef));

    SFREE(set);
    SFREE(meminset);
    SFREE(rmeminset);
    mt = mi[1] + 1;
    NNEW(heading, char, 66 * nheading_);

    nzs_ = 20000000;

    nload = 0;
    nbody = 0;
    nforc = 0;
    nboun = 0;
    nk = 0;
    nmpc = 0;
    nam = 0;

    nlabel = 48;

    while (istat >= 0)
    {
        fflush(stdout);

        /* in order to reduce the number of variables to be transferred to
       the subroutines, the max. field sizes are (for most fields) copied
       into the real sizes */

       nzs[1] = nzs_;

       if ((istep == 0) || (irstrt[0] < 0))
       {
            ne = ne_;
            nset = nset_;
            nmat = nmat_;
            norien = norien_;
            ntrans = ntrans_;
            ntie = ntie_;

            /* allocating space before the first step */
            /* coordinates and topology */

            NNEW(co, double, 3 * nk_);
            NNEW(kon, ITG, nkon_);
            NNEW(ipkon, ITG, ne_);
            NNEW(lakon, char, 8 * ne_);
            NNEW(design, double, ne_);

            /* property cards */
            if (nprop_ > 0)
            {
                NNEW(ielprop, ITG, ne_);
                for (i = 0; i < ne_; i++)
                {
                    ielprop[i] = -1;
                }
                NNEW(prop, double, nprop_);
            }

            /* fields for 1-D and 2-D elements */
            if ((ne1d != 0) || (ne2d != 0))
                {
                    NNEW(iponor, ITG, 2 * nkon_);
                    for (i = 0; i < 2 * nkon_; i++)
                        iponor[i] = -1;
                    NNEW(xnor, double, 36 * ne1d + 24 * ne2d);
                    NNEW(knor, ITG, 24 * (ne1d + ne2d) * (mi[2] + 1));
                    NNEW(thickn, double, 2 * nk_);
                    NNEW(thicke, double, mi[2] * nkon_);
                    NNEW(offset, double, 2 * ne_);
                    NNEW(iponoel, ITG, nk_);
                    NNEW(inoel, ITG, 9 * ne1d + 24 * ne2d);
                    NNEW(rig, ITG, nk_);
                    if (infree[2] == 0)
                        infree[2] = 1;
                 }

            /* SPC's */
            NNEW(nodeboun, ITG, nboun_);
            NNEW(ndirboun, ITG, nboun_);
            NNEW(typeboun, char, nboun_ + 1);
            if ((istep == 0) || ((irstrt[0] < 0) && (nam_ > 0)))
                NNEW(iamboun, ITG, nboun_);
            NNEW(xboun, double, nboun_);
            NNEW(ikboun, ITG, nboun_);
            NNEW(ilboun, ITG, nboun_);

            /* MPC's */
            NNEW(ipompc, ITG, nmpc_);
            NNEW(nodempc, ITG, 3 * memmpc_);

            for (i = 0; i < 3 * memmpc_; i += 3)
            {
                nodempc[i + 2] = i / 3 + 2;
            }
            nodempc[3 * memmpc_ - 1] = 0;
            NNEW(coefmpc, double, memmpc_);
            NNEW(labmpc, char, 20 * nmpc_ + 1);
            NNEW(ikmpc, ITG, nmpc_);
            NNEW(ilmpc, ITG, nmpc_);
            NNEW(fmpc, double, nmpc_);

            /* nodal loads */
            NNEW(nodeforc, ITG, 2 * nforc_);
            NNEW(ndirforc, ITG, nforc_);
            if ((istep == 0) || ((irstrt[0] < 0) && (nam_ > 0)))
                NNEW(iamforc, ITG, nforc_);

            NNEW(idefforc, ITG, nforc_);
            NNEW(xforc, double, nforc_);
            NNEW(ikforc, ITG, nforc_);
            NNEW(ilforc, ITG, nforc_);

            /* distributed facial loads */
            NNEW(nelemload, ITG, 2 * nload_);
            if ((istep == 0) || ((irstrt[0] < 0) && (nam_ > 0)))
                NNEW(iamload, ITG, 2 * nload_);

            NNEW(idefload, ITG, nload_);
            NNEW(sideload, char, 20 * nload_);
            NNEW(xload, double, 2 * nload_);

            /* distributed volumetric loads */
            NNEW(cbody, char, 81 * nbody_);
            NNEW(idefbody, ITG, nbody_);
            NNEW(ibody, ITG, 3 * nbody_);
            NNEW(xbody, double, 7 * nbody_);
            NNEW(xbodyold, double, 7 * nbody_);

            /* printing output */
            NNEW(prlab, char, 6 * nprint_);
            NNEW(prset, char, 81 * nprint_);

            /* set definitions */
            NNEW(set, char, 81 * nset);
            NNEW(istartset, ITG, nset);
            NNEW(iendset, ITG, nset);
            NNEW(ialset, ITG, nalset);

            /* (hyper)elastic constants */
            NNEW(elcon, double, (ncmat_ + 1) * ntmat_ * nmat);
            NNEW(nelcon, ITG, 2 * nmat);

            /* density */
            NNEW(rhcon, double, 2 * ntmat_ * nmat);
            NNEW(nrhcon, ITG, nmat);

            /* damping */
            if (ndamp > 0)
            {
                NNEW(dacon, double, nmat);
            }

            /* specific heat */
            NNEW(shcon, double, 4 * ntmat_ * nmat);
            NNEW(nshcon, ITG, nmat);

            /* thermal expansion coefficients */
            NNEW(alcon, double, 7 * ntmat_ * nmat);
            NNEW(nalcon, ITG, 2 * nmat);
            NNEW(alzero, double, nmat);

            /* conductivity */
            NNEW(cocon, double, 7 * ntmat_ * nmat);
            NNEW(ncocon, ITG, 2 * nmat);

            /* isotropic and kinematic hardening coefficients*/

            if (npmat_ > 0)
            {
                NNEW(plicon, double, (2 * npmat_ + 1) * ntmat_ * nmat);
                NNEW(nplicon, ITG, (ntmat_ + 1) * nmat);
                NNEW(plkcon, double, (2 * npmat_ + 1) * ntmat_ * nmat);
                NNEW(nplkcon, ITG, (ntmat_ + 1) * nmat);
            }

            /* linear dynamic properties */
            NNEW(xmodal, double, 11 + nevdamp_);
            xmodal[10] = nevdamp_ + 0.5;


            /* internal state variables (nslavs is needed for restart
            calculations) */

            if (mortar != 1)
            {
                NNEW(xstate, double, nstate_ *mi[0] * (ne + nslavs));
                nxstate = nstate_ * mi[0] * (ne + nslavs);
            }
            else if (mortar == 1)
            {
                NNEW(xstate, double, nstate_ *mi[0] * (ne + nintpoint));
                nxstate = nstate_ * mi[0] * (ne + nintpoint);
            }

            /* material orientation */
            if ((istep == 0) || ((irstrt[0] < 0) && (norien > 0)))
            {
                NNEW(orname, char, 80 * norien);
                NNEW(orab, double, 7 * norien);
                NNEW(ielorien, ITG, mi[2] * ne_);
            }

            /* transformations */
            if ((istep == 0) || ((irstrt[0] < 0) && (ntrans > 0)))
            {
                NNEW(trab, double, 7 * ntrans);
                NNEW(inotr, ITG, 2 * nk_);
            }

            /* amplitude definitions */
            if ((istep == 0) || ((irstrt[0] < 0) && (nam_ > 0)))
            {
                NNEW(amname, char, 80 * nam_);
                NNEW(amta, double, 2 * namtot_);
                NNEW(namta, ITG, 3 * nam_);
            }

            if ((istep == 0) || ((irstrt[0] < 0) && (ithermal[0] > 0)))
            {
                NNEW(t0, double, nk_);
                NNEW(t1, double, nk_);
                if ((ne1d != 0) || (ne2d != 0))
                {
                    NNEW(t0g, double, 2 * nk_);
                    NNEW(t1g, double, 2 * nk_);
                }
            }

            /* the number in next line is NOT 1.2357111317 -> points
            to user input; instead it is a generic nonzero
            initialization */

            if (istep == 0)
            {
                DMEMSET(t0, 0, nk_, 1.2357111319);
                DMEMSET(t1, 0, nk_, 1.2357111319);
            }

            if ((istep == 0) || ((irstrt[0] < 0) && (ithermal[0] > 0) && (nam_ > 0)))
                NNEW(iamt1, ITG, nk_);

            if ((istep == 0) || ((irstrt[0] < 0) && (iprestr > 0)))
                NNEW(prestr, double, 6 * mi[0] * ne_);

            NNEW(vold, double, mt *nk_);
            NNEW(veold, double, mt *nk_);

            /* CFD-results */
            NNEW(vel, double, 8 * nef);
            NNEW(velo, double, 8 * nef);
            NNEW(veloo, double, 8 * nef);

            NNEW(ielmat, ITG, mi[2] * ne_);

            NNEW(matname, char, 80 * nmat);

            NNEW(filab, char, 87 * nlabel);

            /* tied constraints */
            if (ntie_ > 0)
            {
                NNEW(tieset, char, 243 * ntie_);
                NNEW(tietol, double, 3 * ntie_);
                NNEW(cs, double, 17 * ntie_);
            }

            /* objectives for sensitivity analysis */
            if ((ncs_ > 0) || (npt_ > 0))
            {
                if (2 * npt_ > 24 * ncs_)
                {
                    NNEW(ics, ITG, 2 * npt_);
                }
                else
                {
                    NNEW(ics, ITG, 24 * ncs_);
                }
                if (npt_ > 30 * ncs_)
                {
                NNEW(dcs, double, npt_);
                }
                else
                {
                NNEW(dcs, double, 30 * ncs_);
                }
            }


            /* temporary fields for cyclic symmetry calculations */

            if ((ncs_ > 0) || (npt_ > 0))
            {
                if (2 * npt_ > 24 * ncs_)
                {       
                    NNEW(ics, ITG, 2 * npt_);
                }
                else
                {
                NNEW(ics, ITG, 24 * ncs_);
                }
                if (npt_ > 30 * ncs_)
                {
                    NNEW(dcs, double, npt_);
                }
                else
                {
                    NNEW(dcs, double, 30 * ncs_);
                }
            }

            /* slave faces */
            NNEW(islavsurf, ITG, 2 * ifacecount + 2);
       }
       else
       {
            /* allocating and reallocating space for subsequent steps */
            if ((nmethod != 4) && (nmethod != 5) && (nmethod != 8) && (nmethod != 9) &&
            ((abs(nmethod) != 1) || (iperturb[0] < 2)))
            {
                NNEW(veold, double, mt *nk_);
            }
            else
            {
                RENEW(veold, double, mt *nk_);
                DMEMSET(veold, mt * nk, mt * nk_, 0.);
            }

            RENEW(vold, double, mt *nk_);
            DMEMSET(vold, mt * nk, mt * nk_, 0.);

            RENEW(nodeboun, ITG, nboun_);
            RENEW(ndirboun, ITG, nboun_);
            RENEW(typeboun, char, nboun_ + 1);
            RENEW(xboun, double, nboun_);
            RENEW(ikboun, ITG, nboun_);
            RENEW(ilboun, ITG, nboun_);

            RENEW(nodeforc, ITG, 2 * nforc_);
            RENEW(ndirforc, ITG, nforc_);
            NNEW(idefforc, ITG, nforc_);
            RENEW(xforc, double, nforc_);
            RENEW(ikforc, ITG, nforc_);
            RENEW(ilforc, ITG, nforc_);

            RENEW(nelemload, ITG, 2 * nload_);
            NNEW(idefload, ITG, nload_);
            RENEW(sideload, char, 20 * nload_);
            RENEW(xload, double, 2 * nload_);

            RENEW(cbody, char, 81 * nbody_);
            NNEW(idefbody, ITG, nbody_);
            RENEW(ibody, ITG, 3 * nbody_);
            RENEW(xbody, double, 7 * nbody_);
            RENEW(xbodyold, double, 7 * nbody_);

            for (i = 7 * nbodyold; i < 7 * nbody_; i++)
                xbodyold[i] = 0;

            if (nam > 0)
            {
                RENEW(iamforc, ITG, nforc_);
                RENEW(iamload, ITG, 2 * nload_);
                RENEW(iamboun, ITG, nboun_);
                RENEW(amname, char, 80 * nam_);
                RENEW(amta, double, 2 * namtot_);
                RENEW(namta, ITG, 3 * nam_);
            }

            RENEW(ipompc, ITG, nmpc_);
            RENEW(labmpc, char, 20 * nmpc_ + 1);
            RENEW(ikmpc, ITG, nmpc_);
            RENEW(ilmpc, ITG, nmpc_);
            RENEW(fmpc, double, nmpc_);

            if (ntrans > 0)
            {
                RENEW(inotr, ITG, 2 * nk_);
                DMEMSET(inotr, 2 * nk, 2 * nk_, 0);
            }

            RENEW(co, double, 3 * nk_);
            DMEMSET(co, 3 * nk, 3 * nk_, 0.);

            if (ithermal[0] != 0)
            {
                RENEW(t0, double, nk_);
                DMEMSET(t0, nk, nk_, 0.);
                RENEW(t1, double, nk_);
                DMEMSET(t1, nk, nk_, 0.);

                if ((ne1d != 0) || (ne2d != 0))
                {
                    RENEW(t0g, double, 2 * nk_);
                    DMEMSET(t0g, 2 * nk, 2 * nk_, 0.);
                    RENEW(t1g, double, 2 * nk_);
                    DMEMSET(t1g, 2 * nk, 2 * nk_, 0.);
                }

                if (nam > 0)
                {
                    RENEW(iamt1, ITG, nk_);
                }
            }
       }

       /* allocation of fields in the restart file */

       if (irstrt[0] < 0)
       {
            NNEW(nodebounold, ITG, nboun_);
            NNEW(ndirbounold, ITG, nboun_);
            NNEW(xbounold, double, nboun_);
            NNEW(xforcold, double, nforc_);
            NNEW(xloadold, double, 2 * nload_);

            if (ithermal[0] != 0)
                NNEW(t1old, double, nk_);

            NNEW(sti, double, 6 * mi[0] * ne);
            NNEW(eme, double, 6 * mi[0] * ne);

            if (nener == 1)
                NNEW(ener, double, mi[0] * ne * 2);
            if (mcs > ntie_)
                RENEW(cs, double, 17 * mcs);
            if (mortar == 1)
            {
                NNEW(pslavsurf, double, 3 * nintpoint);
                NNEW(clearini, double, 3 * 9 * ifacecount);
            }
       }

        nenerold = nener;
        nkold = nk;

        /* opening the eigenvalue file and checking for cyclic symmetry */
        strcpy(fneig, jobnamec);
        strcat(fneig, ".eig");
        cyclicsymmetry = 0;
        if ((f1 = fopen(fneig, "rb")) != NULL)
        {
            if (fread(&cyclicsymmetry, sizeof(ITG), 1, f1) != 1)
            {
                printf("*ERROR reading the information whether cyclic symmetry is involved in the eigenvalue file");
                exit(0);
            }
            fclose(f1);
        }

        nmpcold = nmpc;

        /* reading the input file */

        if (istep == 0)
            mortar = -1;

        FORTRAN(betacalinput, (co, &nk, kon, ipkon, lakon, &nkon, &ne,
                       nodeboun, ndirboun, xboun, &nboun,
                       ipompc, nodempc, coefmpc, &nmpc, &nmpc_, nodeforc, ndirforc, xforc, &nforc,
                       &nforc_, nelemload, sideload, xload, &nload, &nload_,
                       &nprint, prlab, prset, &mpcfree, &nboun_, mei, set, istartset, iendset,
                       ialset, &nset, &nalset, elcon, nelcon, rhcon, nrhcon, alcon, nalcon,
                       alzero, t0, t1, matname, ielmat, orname, orab, ielorien, amname,
                       amta, namta, &nam, &nmethod, iamforc, iamload, iamt1,
                       ithermal, iperturb, &istat, &istep, &nmat, &ntmat_, &norien, prestr,
                       &iprestr, &isolver, fei, veold, timepar,
                       xmodal, filab, jout, &nlabel, &idrct,
                       jmax, &iexpl, &alpha, iamboun, plicon, nplicon,
                       plkcon, nplkcon, &iplas, &npmat_, mi, &nk_, trab, inotr, &ntrans,
                       ikboun, ilboun, ikmpc, ilmpc, ics, dcs, &ncs_, &namtot_, cs, &nstate_,
                       &ncmat_, &iumat, &mcs, labmpc, iponor, xnor, knor, thickn, thicke,
                       ikforc, ilforc, offset, iponoel, inoel, rig, infree, nshcon, shcon,
                       cocon, ncocon, physcon, &nflow,
                       ctrl, &maxlenmpc, &ne1d, &ne2d, &nener, vold, nodebounold,
                       ndirbounold, xbounold, xforcold, xloadold, t1old, eme,
                       sti, ener, xstate, jobnamec, irstrt, &ttime,
                       qaold, output, typeboun, inpc, ipoinp, inp, tieset, tietol,
                       &ntie, fmpc, cbody, ibody, xbody, &nbody, &nbody_, xbodyold, &nam_,
                       ielprop, &nprop, &nprop_, prop, &itpamp, &iviewfile, ipoinpc,
                       &nslavs, t0g, t1g, &network, &cyclicsymmetry, idefforc, idefload,
                       idefbody, &mortar, &ifacecount, islavsurf, pslavsurf, clearini,
                       heading, &iaxial, &nobject, objectset, &nprint_, iuel, &nuel_,
                       nodempcref, coefmpcref, ikmpcref, &memmpcref_, &mpcfreeref,
                       &maxlenmpcref, &memmpc_, &isens, &namtot, &nstam, dacon, vel, &nef,
                       velo, veloo));










            








       


    }
}