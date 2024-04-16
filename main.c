#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"
#include <time.h>

int main(int argc, char *argv[])
{
    FILE *f1;

    char *sideload = NULL, *set = NULL, *matname = NULL, *orname = NULL, *amname = NULL,
       *filab = NULL, *lakon = NULL, *labmpc = NULL, *prlab = NULL, *prset = NULL,
       jobnamec[660] = "", jobnamef[132] = "", output[4] = "asc", *typeboun = NULL,
       *inpc = NULL, *tieset = NULL, *cbody = NULL, fneig[132] = "", *sideloadtemp = NULL,
       kind1[2] = "T", kind2[2] = "T", *heading = NULL, *objectset = NULL;

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

    ITG *kon = NULL, *nodeboun = NULL, *ndirboun = NULL, *ipompc = NULL,
      *nodempc = NULL, *nodeforc = NULL, *ndirforc = NULL,
      *nelemload = NULL, im, *inodesd = NULL, nload1, *idefforc = NULL,
      *nactdof = NULL, *icol = NULL, *ics = NULL,
      *jq = NULL, *mast1 = NULL, *irow = NULL, *rig = NULL, *idefbody = NULL,
      *ikmpc = NULL, *ilmpc = NULL, *ikboun = NULL, *ilboun = NULL,
      *nreorder = NULL, *ipointer = NULL, *idefload = NULL,
      *istartset = NULL, *iendset = NULL, *ialset = NULL, *ielmat = NULL,
      *ielorien = NULL, *nrhcon = NULL, *nodebounold = NULL, *ndirbounold = NULL,
      *nelcon = NULL, *nalcon = NULL, *iamforc = NULL, *iamload = NULL,
      *iamt1 = NULL, *namta = NULL, *ipkon = NULL, *iamboun = NULL,
      *nplicon = NULL, *nplkcon = NULL, *inotr = NULL, *iponor = NULL, *knor = NULL,
      *ikforc = NULL, *ilforc = NULL, *iponoel = NULL, *inoel = NULL, *nshcon = NULL,
      *ncocon = NULL, *ibody = NULL, *ielprop = NULL, *islavsurf = NULL,
      *ipoinpc = NULL, mt, nxstate, nload0, iload, *iuel = NULL;

    ITG *meminset = NULL, *rmeminset = NULL;

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

}