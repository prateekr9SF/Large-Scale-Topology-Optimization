/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "mkl_spblas.h"
#include "CalculiX.h"
#include "pardiso.h"

#ifdef PROFILING_ON
#include <TAU.h>
#endif

ITG *icolpardiso=NULL,*pointers=NULL,iparm[64];
long long pt[64];
double *aupardiso=NULL;
/* double dparm[64];  not used */
ITG nthread_mkl=0;
int g_useNK=0;   /* --NK: Newton-Krylov in nonlingeo (see end of file) */

/* factor kept by pardiso_factor until pardiso_cleanup */
static ITG factor_alive=0,alive_neq=0,alive_symmetryflag=0;

/* Newton-Krylov state (separate PARDISO handle, see end of file) */
static ITG nk_alive=0;
static void nk_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs,ITG L);
static void nk_release(void);
#define NK_TOL_NEWTON    3   /* FGMRES rel. residual 1e-3, Newton solves      */
#define NK_TOL_SOLVE     8   /* FGMRES rel. residual 1e-8, adjoint solves     */
#define NK_KRYLOV_MAXIT 30   /* FGMRES iterations before refactoring          */
#define NK_KRYLOV_M     15   /* restart length (memory ~ (2m+1)*neq doubles)  */
#define NK_ORDERING      3   /* iparm[1]: 2 METIS, 3 parallel METIS           */
/* char envMKL[32];   moved to pardiso.h */
//char envMKL[32];

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_factor
 * ----------------------------------------------------------------------
 * Purpose:
 *    Prepares and factors the global linear system matrix using the
 *    Intel MKL PARDISO sparse direct solver. This function converts
 *    the CalculiX internal sparse matrix format into the CSR format
 *    required by PARDISO, applies an optional shift (K - sigma*M),
 *    configures threading, initializes PARDISO control parameters,
 *    and performs the symbolic and numerical factorization (phase = 12).
 *
 * Arguments:
 *    ad       - [in]  Diagonal entries of stiffness matrix K.
 *    au       - [in]  Off-diagonal entries of stiffness matrix K
 *                     (format depends on symmetryflag and inputformat).
 *    adb      - [in]  Diagonal entries of auxiliary matrix (e.g., mass) M.
 *    aub      - [in]  Off-diagonal entries of M.
 *    sigma    - [in]  Scalar shift; if non-zero, factorizes K - sigma*M.
 *
 *    icol     - [in]  Number of off-diagonal entries in each column of K.
 *    irow     - [in]  Row indices for off-diagonal terms in K.
 *
 *    neq      - [in]  Number of equations (matrix dimension).
 *    nzs      - [in]  Number of off-diagonal nonzeros.
 *
 *    symmetryflag - [in]  0 = symmetric, 1 = unsymmetric matrix.
 *    inputformat   - [in]  Storage scheme used by CalculiX
 *                          (1 = split triangular, 3 = column-based unsymmetric).
 *
 *    jq       - [in]  Indices of upper triangular entries (inputformat = 1).
 *    nzs3     - [in]  Auxiliary offset for upper triangular storage.
 *
 * Output:
 *    No explicit return. Allocates and fills global arrays:
 *       pointers     - CSR row pointer array (1-based)
 *       icolpardiso  - CSR column indices
 *       aupardiso    - CSR matrix values
 *       pt[], iparm[] configured for later solve call
 *
 * Notes:
 *    - Sets MKL_NUM_THREADS based on CCX_NPROC_EQUATION_SOLVER,
 *      MKL_NUM_THREADS, or OMP_NUM_THREADS (first call only).
 *    - Performs PARDISO symbolic + numeric factorization (phase = 12).
 *    - Must be followed by pardiso_solve() and pardiso_cleanup().
 *
 * ----------------------------------------------------------------------
 */

static void pardiso_threads(void)
	{
  		char *env;
  		ITG nthread,nthread_v;

		/* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
   			must be done once  */
  		if (nthread_mkl == 0) 
		{
    		nthread=1;
    		env=getenv("MKL_NUM_THREADS");
    		if(env) 
			{
      			nthread=atoi(env);
			}
    		else 
			{
      			env=getenv("OMP_NUM_THREADS");
      			if(env) 
				{
					nthread=atoi(env);
				}
    		}
    		env=getenv("CCX_NPROC_EQUATION_SOLVER");
    		if(env) 
			{
      			nthread_v=atoi(env);
      			if (nthread_v <= nthread) 
				{
					nthread=nthread_v;
				}
    		}
    		if (nthread < 1) 
			{
				nthread=1;
			}
    		sprintf(envMKL,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);  
    		putenv(envMKL);
    		nthread_mkl=nthread;
  		}
    
//  		printf(" number of threads =% d\n\n",nthread_mkl);
	}

/* CalculiX storage -> PARDISO CSR (pointers, icolpardiso, aupardiso);
   returns the PARDISO matrix type */
static ITG pardiso_csr(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
	{
  		ITG i,j,k,l,mtype,*irowpardiso=NULL,kflag,kstart,n,k2,ndim;

  		if(*symmetryflag==0)
		{

      		/* symmetric matrix; the subdiagonal entries are stored
         	column by column in au, the diagonal entries in ad;
         	pardiso needs the entries row per row */      

      		mtype=-2;

      
      		ndim=*neq+*nzs;
      
      		NNEW(pointers,ITG,*neq+1);
      		NNEW(icolpardiso,ITG,ndim);
      		NNEW(aupardiso,double,ndim);
      
      		k=ndim;
      		l=*nzs;
      
      		if(*sigma==0.)
			{
	  			pointers[*neq]=ndim+1;

	  			for(i=*neq-1;i>=0;--i)
				{
	      			for(j=0;j<icol[i];++j)
					{
		  				icolpardiso[--k]=irow[--l];
		  				aupardiso[k]=au[l];
	      			}
	      			pointers[i]=k--;
	      			icolpardiso[k]=i+1;
	      			aupardiso[k]=ad[i];
	  			}
      		}
      		else
			{
	  			pointers[*neq]=ndim+1;
	  			for(i=*neq-1;i>=0;--i)
				{
	      			for(j=0;j<icol[i];++j)
					{
		  				icolpardiso[--k]=irow[--l];
		  				aupardiso[k]=au[l]-*sigma*aub[l];
	      			}
	      			pointers[i]=k--;
	      			icolpardiso[k]=i+1;
	      			aupardiso[k]=ad[i]-*sigma*adb[i];
	  			}
      		}
  		}
		else
		{
      		mtype=11;

      		if(*inputformat==3)
			{

          		/* off-diagonal terms  are stored column per
            	 column from top to bottom in au;
             	diagonal terms are stored in ad  */

	  			ndim=*neq+*nzs;
	  			NNEW(pointers,ITG,*neq+1);
	  			NNEW(irowpardiso,ITG,ndim);	  
	  			NNEW(icolpardiso,ITG,ndim);
	  			NNEW(aupardiso,double,ndim);
	  
	  			k=0;
	  			k2=0;
	  			for(i=0;i<*neq;i++)
				{
	      			for(j=0;j<icol[i];j++)
					{
		  				if(au[k]>1.e-12||au[k]<-1.e-12)
						{
		      				icolpardiso[k2]=i+1;
		      				irowpardiso[k2]=irow[k];
		      				aupardiso[k2]=au[k];
		      				k2++;		  
		  				}
		  				k++;	      
	      			}	  
	  			}  
	  			/* diagonal terms */  
	  			for(i=0;i<*neq;i++)
				{
	      			icolpardiso[k2]=i+1;
	      			irowpardiso[k2]=i+1;
	      			aupardiso[k2]=ad[i];
	      			k2++;	  
	  			}
	  			ndim=k2;
	  
	  			/* pardiso needs the entries row per row; so sorting is
	     		needed */ 
	  
	  			kflag=2;
	  			FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			    	&ndim,&kflag));
	  
	  			/* sorting each row */
	  
	  			k=0;
	  			pointers[0]=1;
	  			for(i=0;i<*neq;i++)
				{
	      			j=i+1;
	      			kstart=k;
	      			do
					{
		  				if(irowpardiso[k]!=j )
						{
		      				n=k-kstart;		  
		      				FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				       		&n,&kflag));
		     	 			pointers[i+1]=k+1;
		      				break;  
		  				}
						else
						{
		      				if(k+1==ndim)
							{
			  					n=k-kstart+1;	  
			  					FORTRAN(isortid,(&icolpardiso[kstart],
                            	      &aupardiso[kstart],&n,&kflag));
			  					break;	       
		      				}
							else
							{
			  					k++;	       
		      				}  
		  				}
	      			}
					while(1);
	  			}
	  			pointers[*neq]=ndim+1;
	  			SFREE(irowpardiso);

      		}
			else if(*inputformat==1)
			{
	  
          		/* lower triangular matrix is stored column by column in
            	 au, followed by the upper triangular matrix row by row;
             		the diagonal terms are stored in ad */

          		/* reordering lower triangular matrix */

	  			ndim=*nzs;
	  			NNEW(pointers,ITG,*neq+1);
	  			NNEW(irowpardiso,ITG,ndim);
	  			NNEW(icolpardiso,ITG,ndim);
	  			NNEW(aupardiso,double,ndim);
				
	  			k=0;
	  			for(i=0;i<*neq;i++)
				{
	      			for(j=0;j<icol[i];j++)
					{
		  				icolpardiso[k]=i+1;
		  				irowpardiso[k]=irow[k];
		  				aupardiso[k]=au[k];
		  				k++;
	      			}
	  			}
	  
	  			/* pardiso needs the entries row per row; so sorting is
	     			needed */
	  
	  			kflag=2;
	  			FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
	  				&ndim,&kflag));
	  
	  			/* sorting each row */
	  
	  			k=0;
	  			pointers[0]=1;
	  			if(ndim>0)
				{
	      			for(i=0;i<*neq;i++)
					{
		  				j=i+1;
		  				kstart=k;
		  				do
						{
	  	      				/* end of row reached */

		      				if(irowpardiso[k]!=j)
							{
			  					n=k-kstart;
			  					FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					   				&n,&kflag));
			  					pointers[i+1]=k+1;
			  					break;
		      				}
							else
							{
		          				/* end of last row reached */
			  					if(k+1==ndim)
								{
			      					n=k-kstart+1;
			      					FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					       				&n,&kflag));
			      					break;
			  					}
								else
								{
			      					/* end of row not yet reached */
			      					k++;
			  					}
		      				
							}
		  				}
						while(1);
	      			}
	  			}
	  			pointers[*neq]=ndim+1;
	  			SFREE(irowpardiso);

	  			/* composing the matrix: lower triangle + diagonal + upper triangle */

	  			ndim=*neq+2**nzs;
	  			RENEW(icolpardiso,ITG,ndim);
	  			RENEW(aupardiso,double,ndim);
	  			k=ndim;

	  			for(i=*neq-1;i>=0;i--)
				{
	      			l=k+1;
	      			for(j=jq[i+1]-1;j>=jq[i];j--)
					{
		  				icolpardiso[--k]=irow[j-1];
		  				aupardiso[k]=au[j+*nzs3-1];
	      			}

	      			icolpardiso[--k]=i+1;
	      			aupardiso[k]=ad[i];

	      			for(j=pointers[i+1]-1;j>=pointers[i];j--)
					{
		  				icolpardiso[--k]=icolpardiso[j-1];
		  				aupardiso[k]=aupardiso[j-1];
	      			}

	      			pointers[i+1]=l;
	  			}
	  			pointers[0]=1;
      		}
  		}

  		return mtype;
	}

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
	{

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_TIMER(t_pardiso_factor,"PARDISO: Factorization","",TAU_USER);
	//		TAU_PROFILE_START(t_pardiso_factor);
    //	#endif

  		ITG i,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
      		msglvl=0,error=0;
  		double *b=NULL,*x=NULL;

  		iparm[0]=0;
  		pardiso_threads();

  		for(i=0;i<64;i++)
		{
			pt[i]=0;
		}

  		mtype=pardiso_csr(ad,au,adb,aub,sigma,icol,irow,neq,nzs,
  			  symmetryflag,inputformat,jq,nzs3);

  		FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
			   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  		factor_alive=1;
  		alive_neq=*neq;
  		alive_symmetryflag=*symmetryflag;

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_STOP(t_pardiso_factor);
    //	#endif
				   
  		return;
	}
/*
 * ----------------------------------------------------------------------
 * Function: pardiso_solve
 * ----------------------------------------------------------------------
 * Purpose:
 *    Solves the linear system A * x = b using the factorization
 *    produced earlier by pardiso_factor(). The solution overwrites b.
 *    This function invokes PARDISO with phase = 33 (solve only).
 *
 * Arguments:
 *    b            - [in/out]  Right-hand side vector(s). On return,
 *                              contains the computed solution(s).
 *    neq          - [in]      Number of equations.
 *    symmetryflag - [in]      0 = symmetric system, 1 = unsymmetric.
 *    nrhs         - [in]      Number of right-hand sides.
 *
 * Output:
 *    - Solution is written into b.
 *    - No matrix memory is modified; uses the existing PARDISO pt[] data.
 *
 * Notes:
 *    - Requires that pardiso_factor() has already been called.
 *    - Allocates a temporary solution vector x and copies results into b.
 *    - Thread count is inherited from the first call to pardiso_factor().
 *
 * ----------------------------------------------------------------------
 */

void pardiso_solve(double *b, ITG *neq,ITG *symmetryflag,ITG *nrhs)
{

//	#ifdef PROFILING_ON
//        TAU_PROFILE_TIMER(t_pardiso_solve,"PARDISO: Solve","",TAU_USER);
//		TAU_PROFILE_START(t_pardiso_solve);
//    #endif

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
    	msglvl=0,i,error=0;
  double *x=NULL;

  /* matrix stored by pardiso_set_matrix / pardiso_factor_solve */
  if(!factor_alive&&nk_alive)
  {
    nk_solve(b,neq,symmetryflag,nrhs,NK_TOL_SOLVE);
    return;
  }

  if(*symmetryflag==0)
  {
    mtype=-2;
  }
  else
  {
    mtype=11;
  }

  iparm[0]=0;

	/* pardiso_factor has been called befor, MKL_NUM_THREADS=nthread_mkl is set*/
  NNEW(x,double,*nrhs**neq);

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,nrhs,iparm,&msglvl,
                   b,x,&error));

  for(i=0;i<*nrhs**neq;i++){b[i]=x[i];}
  SFREE(x);

//	#ifdef PROFILING_ON
//    	TAU_PROFILE_STOP(t_pardiso_solve);
//	#endif

  return;
}

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_cleanup
 * ----------------------------------------------------------------------
 * Purpose:
 *    Releases all PARDISO internal memory associated with the matrix
 *    factorization and frees the CSR arrays allocated during
 *    pardiso_factor(). Uses PARDISO phase = -1 (memory release).
 *
 * Arguments:
 *    neq          - [in]  Number of equations.
 *    symmetryflag - [in]  0 = symmetric matrix, 1 = unsymmetric.
 *
 * Output:
 *    No return value. Frees global arrays:
 *        icolpardiso
 *        aupardiso
 *        pointers
 *    and clears all PARDISO internal structures in pt[].
 *
 * Notes:
 *    - Must be called after all solves are complete.
 *    - Safe to call even in one-shot solve (used by pardiso_main()).
 *
 * ----------------------------------------------------------------------
 */


void pardiso_cleanup(ITG *neq,ITG *symmetryflag)
{

//	#ifdef PROFILING_ON
//        TAU_PROFILE_TIMER(t_pardiso_cleanup,
//                          "PARDISO: Cleanup",
//                          "",
//                          TAU_USER);

//        TAU_PROFILE_START(t_pardiso_cleanup);
//    #endif

  ITG maxfct=1,mnum=1,phase=-1,*perm=NULL,nrhs=1,mtype,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0)
  {
    mtype=-2;
  }
  else
  {
    mtype=11;
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  SFREE(icolpardiso);
  SFREE(aupardiso);
  SFREE(pointers);

  factor_alive=0;

//    #ifdef PROFILING_ON
//        TAU_PROFILE_STOP(t_pardiso_cleanup);
//    #endif

  return;
}

/* 1 if a matrix is held for pardiso_solve (pardiso_factor or NK) */
ITG pardiso_factor_alive(void)
{
  return factor_alive||nk_alive;
}

/* release the held factorization(s), if any */
void pardiso_release(void)
{
  if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
  nk_release();
}

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_main
 * ----------------------------------------------------------------------
 * Purpose:
 *    Convenience wrapper that performs the complete PARDISO workflow
 *    for a single solve:
 *       1) Assemble and factor the system matrix (pardiso_factor),
 *       2) Solve A * x = b for the given right-hand side(s)
 *          (pardiso_solve),
 *       3) Release all PARDISO and matrix-related memory
 *          (pardiso_cleanup).
 *
 * Arguments:
 *    ad       - [in]  Diagonal entries of stiffness matrix K.
 *    au       - [in]  Off-diagonal entries of stiffness matrix K.
 *    adb      - [in]  Diagonal entries of auxiliary matrix (e.g., mass) M.
 *    aub      - [in]  Off-diagonal entries of M.
 *    sigma    - [in]  Scalar shift; if non-zero, K - sigma*M is factored.
 *
 *    b        - [in/out]  Right-hand side(s). On return, contains the
 *                          solution vector(s).
 *
 *    icol     - [in]  Number of off-diagonal entries in each column of K.
 *    irow     - [in]  Row indices for off-diagonal terms.
 *
 *    neq      - [in]  Number of equations (matrix dimension).
 *    nzs      - [in]  Number of off-diagonal nonzeros.
 *
 *    symmetryflag - [in]  0 = symmetric system, 1 = unsymmetric system.
 *    inputformat   - [in]  CalculiX matrix storage format
 *                          (1 = split triangular, 3 = column-based).
 *
 *    jq       - [in]  Column pointer array for upper triangle
 *                     (used when inputformat = 1).
 *    nzs3     - [in]  Offset into au for upper triangular entries.
 *
 *    nrhs     - [in]  Number of right-hand sides.
 *
 * Output:
 *    - Overwrites b with the solution(s) of the linear system.
 *    - Internally allocates and frees all data structures for PARDISO.
 *
 * Notes:
 *    - Returns immediately if *neq == 0 (empty system).
 *    - Suitable for single load case / one-shot solves.
 *      For repeated solves with the same matrix and different RHS,
 *      call pardiso_factor(), then pardiso_solve() repeatedly,
 *      and finally pardiso_cleanup().
 *
 * ----------------------------------------------------------------------
 */
void pardiso_main(double *ad, double *au, double *adb, double *aub, 
         double *sigma,double *b, ITG *icol, ITG *irow, 
	ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	ITG *jq, ITG *nzs3,ITG *nrhs)
	{

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_TIMER(t_pardiso_total,"PARDISO: Total","",TAU_USER);
    //    	TAU_PROFILE_START(t_pardiso_total);
    //	#endif

  		if(*neq==0)
		{
		//	#ifdef PROFILING_ON
      //  		TAU_PROFILE_STOP(t_pardiso_total);
    	//	#endif
			return;
		}
  		pardiso_factor(ad,au,adb,aub,sigma,icol,irow, 
			neq,nzs,symmetryflag,inputformat,jq,nzs3);

  		pardiso_solve(b,neq,symmetryflag,nrhs);

  		pardiso_cleanup(neq,symmetryflag);

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_STOP(t_pardiso_total);
    //	#endif

  		return;
	}

/*
 * ----------------------------------------------------------------------
 * Newton-Krylov (--NK true, nonlingeo only)
 * ----------------------------------------------------------------------
 * Separate PARDISO handle; the stock functions above are not affected.
 *  - symbolic analysis (phase 11) is redone only when the CSR pattern
 *    changes (FNV-1a hash); numerical factorization is phase 22
 *  - a new matrix is solved by right-preconditioned FGMRES(m) with the
 *    last factorization as preconditioner (phase 33); refactor only if
 *    FGMRES does not converge in NK_KRYLOV_MAXIT iterations
 *  - pardiso_factor_solve: Newton step (tol 1e-NK_TOL_NEWTON)
 *    pardiso_set_matrix:   store the converged K_T without factoring;
 *    pardiso_solve then solves with it to 1e-NK_TOL_SOLVE
 * ----------------------------------------------------------------------
 */

static long long nk_pt[64];
static ITG nk_iparm[64],*nk_pointers=NULL,*nk_icol=NULL;
static double *nk_au=NULL;
static ITG nk_valid=0,nk_stale=0,nk_neq=-1,nk_mtype=0,nk_nnz=-1;
static unsigned long long nk_hash=0ULL;
static long n_symbolic=0,n_factor=0,n_kry_ok=0,n_kry_fail=0,n_kry_iter=0,
            n_backsolve=0;
static double t_symbolic=0.,t_factor=0.,t_krylov=0.,t_backsolve=0.;

static void nk_call(ITG phase,ITG mtype,ITG *neq,ITG nrhs,double *b,
                    double *x,const char *where)
{
  ITG maxfct=1,mnum=1,*perm=NULL,msglvl=0,error=0;
  FORTRAN(pardiso,(nk_pt,&maxfct,&mnum,&mtype,&phase,neq,nk_au,
                   nk_pointers,nk_icol,perm,&nrhs,nk_iparm,&msglvl,
                   b,x,&error));
  if(error!=0){
    printf(" *ERROR in %s: pardiso phase %d returned error=%d\n",
           where,(int)phase,(int)error);
    FORTRAN(stop,());
  }
}

static void nk_set_iparm(ITG mtype)
{
  ITG i;
  for(i=0;i<64;i++) nk_iparm[i]=0;
  nk_iparm[0]=1;
  nk_iparm[1]=NK_ORDERING;
  if(mtype==-2){              /* MKL defaults, real symmetric indefinite */
    nk_iparm[9]=8;nk_iparm[10]=0;nk_iparm[12]=0;nk_iparm[20]=1;
  }else{                      /* MKL defaults, real unsymmetric */
    nk_iparm[9]=13;nk_iparm[10]=1;nk_iparm[12]=1;
  }
  nk_iparm[17]=-1;            /* report nnz(L) */
}

static void nk_release(void)
{
  ITG neq=nk_neq;
  if(!nk_alive) return;
  printf(" PARDISO (NK): symbolic=%ld factorizations=%ld back-solves=%ld"
         " FGMRES ok=%ld (avg %.1f it) failed->refactor=%ld\n"
         "               time: symbolic %.1f s, factor %.1f s, FGMRES %.1f s,"
         " back-solve %.1f s\n",
         n_symbolic,n_factor,n_backsolve,n_kry_ok,
         n_kry_ok>0 ? (double)n_kry_iter/(double)n_kry_ok : 0.,n_kry_fail,
         t_symbolic,t_factor,t_krylov,t_backsolve);
  nk_call(-1,nk_mtype,&neq,1,NULL,NULL,"nk_release");
  SFREE(nk_pointers);SFREE(nk_icol);SFREE(nk_au);
  nk_pointers=NULL;nk_icol=NULL;nk_au=NULL;
  nk_alive=0;nk_valid=0;nk_stale=0;nk_neq=-1;nk_mtype=0;nk_nnz=-1;
  nk_hash=0ULL;
  n_symbolic=n_factor=n_kry_ok=n_kry_fail=n_kry_iter=n_backsolve=0;
  t_symbolic=t_factor=t_krylov=t_backsolve=0.;
}

static unsigned long long nk_pattern_hash(ITG neq,ITG nnz)
{
  unsigned long long h=1469598103934665603ULL;
  ITG i;
  for(i=0;i<=neq;i++){h^=(unsigned long long)nk_pointers[i];h*=1099511628211ULL;}
  for(i=0;i<nnz;i++){h^=(unsigned long long)nk_icol[i];h*=1099511628211ULL;}
  return h;
}

static void nk_numeric(ITG *neq)
{
  double t0=omp_get_wtime(),dt;
  nk_call(22,nk_mtype,neq,1,NULL,NULL,"pardiso factor (NK)");
  nk_valid=1;nk_stale=0;n_factor++;
  dt=omp_get_wtime()-t0;t_factor+=dt;
  printf(" PARDISO: factorization %.2f s\n",dt);
}

/* store a new matrix in the NK handle; symbolic analysis only if the
   pattern changed; returns 1 if it was redone (no usable factor) */
static ITG nk_store(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
{
  ITG i,mtype,nnz,*sp=pointers,*si=icolpardiso;
  double *sa=aupardiso,t0;
  unsigned long long h;

  if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
  pardiso_threads();

  /* pardiso_csr fills the stock arrays; move them to the NK handle */
  mtype=pardiso_csr(ad,au,adb,aub,sigma,icol,irow,neq,nzs,
		    symmetryflag,inputformat,jq,nzs3);
  if(nk_pointers!=NULL){SFREE(nk_pointers);SFREE(nk_icol);SFREE(nk_au);}
  nk_pointers=pointers;nk_icol=icolpardiso;nk_au=aupardiso;
  pointers=sp;icolpardiso=si;aupardiso=sa;

  nnz=nk_pointers[*neq]-1;
  h=nk_pattern_hash(*neq,nnz);
  if(nk_alive&&nk_neq==*neq&&nk_mtype==mtype&&nk_nnz==nnz&&nk_hash==h)
    return 0;

  if(nk_alive){
    ITG neq0=nk_neq;
    nk_call(-1,nk_mtype,&neq0,1,NULL,NULL,"nk_store");
  }
  for(i=0;i<64;i++) nk_pt[i]=0;
  nk_set_iparm(mtype);
  nk_mtype=mtype;
  t0=omp_get_wtime();
  nk_call(11,mtype,neq,1,NULL,NULL,"pardiso symbolic (NK)");
  t_symbolic+=omp_get_wtime()-t0;n_symbolic++;
  nk_alive=1;nk_valid=0;nk_stale=0;
  nk_neq=*neq;nk_nnz=nnz;nk_hash=h;
  return 1;
}

static double vdot(ITG n,const double *x,const double *y)
{
  double s=0.;
  ITG i;
#pragma omp parallel for reduction(+:s)
  for(i=0;i<n;i++) s+=x[i]*y[i];
  return s;
}

/* y = A x with the stored CSR (upper triangle for mtype -2) */
static void nk_matvec(ITG n,const double *x,double *y)
{
  sparse_matrix_t A;
  struct matrix_descr d;
  sparse_status_t st;
  st=mkl_sparse_d_create_csr(&A,SPARSE_INDEX_BASE_ONE,n,n,nk_pointers,
                             nk_pointers+1,nk_icol,nk_au);
  if(st!=SPARSE_STATUS_SUCCESS){
    printf(" *ERROR in nk_matvec: mkl_sparse_d_create_csr=%d\n",(int)st);
    FORTRAN(stop,());
  }
  if(nk_mtype==-2){
    d.type=SPARSE_MATRIX_TYPE_SYMMETRIC;
    d.mode=SPARSE_FILL_MODE_UPPER;
    d.diag=SPARSE_DIAG_NON_UNIT;
  }else{
    d.type=SPARSE_MATRIX_TYPE_GENERAL;
  }
  st=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.0,A,d,x,0.0,y);
  mkl_sparse_destroy(A);
  if(st!=SPARSE_STATUS_SUCCESS){
    printf(" *ERROR in nk_matvec: mkl_sparse_d_mv=%d\n",(int)st);
    FORTRAN(stop,());
  }
}

/*
 * Right-preconditioned flexible GMRES(m) for A x = b, x0 = 0, with the
 * existing factor as M (flexible: MKL may apply iterative refinement in
 * phase 33). Returns the iteration count on ||b-Ax|| <= 10^-L ||b||
 * (true residual), or -1.
 */
static ITG nk_fgmres(ITG *neq,const double *b,double *x,ITG L)
{
  ITG n=*neq,m=NK_KRYLOV_M,maxit=NK_KRYLOV_MAXIT,it=0,i,j,k,jl;
  double tol=pow(10.,-(double)L),bnorm,beta,rnorm,tmp,*V=NULL,*Z=NULL,
         *H=NULL,*cs=NULL,*sn=NULL,*g=NULL,*y=NULL,*w=NULL;

  bnorm=sqrt(vdot(n,b,b));
  for(i=0;i<n;i++) x[i]=0.;
  if(bnorm==0.) return 1;

  NNEW(V,double,(size_t)(m+1)*n);
  NNEW(Z,double,(size_t)m*n);
  NNEW(H,double,(m+1)*m);
  NNEW(cs,double,m);NNEW(sn,double,m);NNEW(g,double,m+1);NNEW(y,double,m);
  NNEW(w,double,n);

  memcpy(V,b,sizeof(double)*n);
  beta=bnorm;

  while(1){
    for(i=0;i<n;i++) V[i]/=beta;
    for(i=0;i<=m;i++) g[i]=0.;
    g[0]=beta;
    jl=-1;

    for(j=0;j<m&&it<maxit;j++){
      double *vj=&V[(size_t)j*n],*zj=&Z[(size_t)j*n],*vj1=&V[(size_t)(j+1)*n];
      nk_call(33,nk_mtype,neq,1,vj,zj,"pardiso preconditioner (NK)");
      nk_matvec(n,zj,w);
      for(k=0;k<=j;k++){                      /* modified Gram-Schmidt */
        double *vk=&V[(size_t)k*n],hk=vdot(n,w,vk);
        H[k+j*(m+1)]=hk;
#pragma omp parallel for
        for(i=0;i<n;i++) w[i]-=hk*vk[i];
      }
      tmp=sqrt(vdot(n,w,w));
      H[j+1+j*(m+1)]=tmp;
      if(tmp>0.){
#pragma omp parallel for
        for(i=0;i<n;i++) vj1[i]=w[i]/tmp;
      }
      for(k=0;k<j;k++){                       /* previous Givens */
        double h0=H[k+j*(m+1)],h1=H[k+1+j*(m+1)];
        H[k+j*(m+1)]  = cs[k]*h0+sn[k]*h1;
        H[k+1+j*(m+1)]=-sn[k]*h0+cs[k]*h1;
      }
      {                                       /* new Givens */
        double h0=H[j+j*(m+1)],h1=H[j+1+j*(m+1)],r=sqrt(h0*h0+h1*h1);
        if(r==0.){cs[j]=1.;sn[j]=0.;}else{cs[j]=h0/r;sn[j]=h1/r;}
        H[j+j*(m+1)]=r;H[j+1+j*(m+1)]=0.;
        g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
      }
      it++;jl=j;
      if(fabs(g[j+1])<=tol*bnorm||tmp==0.) break;
    }

    for(k=jl;k>=0;k--){                       /* x += Z y,  H y = g */
      double sk=g[k];
      for(i=k+1;i<=jl;i++) sk-=H[k+i*(m+1)]*y[i];
      y[k]=(H[k+k*(m+1)]!=0.) ? sk/H[k+k*(m+1)] : 0.;
    }
    for(k=0;k<=jl;k++){
      double *zk=&Z[(size_t)k*n],yk=y[k];
#pragma omp parallel for
      for(i=0;i<n;i++) x[i]+=yk*zk[i];
    }

    nk_matvec(n,x,w);                         /* true residual */
#pragma omp parallel for
    for(i=0;i<n;i++) V[i]=b[i]-w[i];
    rnorm=sqrt(vdot(n,V,V));
    if(rnorm<=tol*bnorm) break;
    if(it>=maxit||rnorm!=rnorm){it=-1;break;}
    beta=rnorm;                               /* restart */
  }

  SFREE(V);SFREE(Z);SFREE(H);SFREE(cs);SFREE(sn);SFREE(g);SFREE(y);SFREE(w);
  return it;
}

/* solve with the stored matrix: back-solve if the factor is current,
   FGMRES on the stale factor otherwise (refactor if FGMRES fails) */
static void nk_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs,ITG L)
{
  ITG i,r,it,mtype=(*symmetryflag==0) ? -2 : 11;
  double *x=NULL,t0,dt;

  if(!nk_alive||!nk_valid||mtype!=nk_mtype||*neq!=nk_neq){
    printf(" *ERROR in pardiso_solve: no matching NK matrix stored\n");
    FORTRAN(stop,());
  }
  NNEW(x,double,*nrhs**neq);
  for(r=0;r<*nrhs;r++){
    double *br=&b[(size_t)r**neq],*xr=&x[(size_t)r**neq];
    if(nk_stale){
      t0=omp_get_wtime();
      it=nk_fgmres(neq,br,xr,L);
      dt=omp_get_wtime()-t0;t_krylov+=dt;
      if(it>0){
        n_kry_ok++;n_kry_iter+=it;
        printf(" PARDISO: old factor reused, FGMRES %d it to 1e-%d, %.2f s\n",
               (int)it,(int)L,dt);
        continue;
      }
      n_kry_fail++;
      printf(" PARDISO: FGMRES not converged in %d it (%.2f s) -> refactor\n",
             (int)NK_KRYLOV_MAXIT,dt);
      nk_numeric(neq);
    }
    t0=omp_get_wtime();
    nk_call(33,nk_mtype,neq,1,br,xr,"pardiso_solve (NK)");
    t_backsolve+=omp_get_wtime()-t0;n_backsolve++;
  }
  for(i=0;i<*nrhs**neq;i++) b[i]=x[i];
  SFREE(x);
}

/* Newton step: store the current K_T and solve it (FGMRES on the last
   factor; factor only if there is none or the pattern changed) */
void pardiso_factor_solve(double *ad, double *au, double *adb, double *aub,
                double *sigma,double *b,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3,ITG *nrhs)
{
  if(nk_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
	      inputformat,jq,nzs3)||!nk_valid){
    nk_numeric(neq);
  }else{
    nk_stale=1;
  }
  nk_solve(b,neq,symmetryflag,nrhs,NK_TOL_NEWTON);
}

/* store a matrix without factoring it; pardiso_solve then solves with
   it by FGMRES on the last factor (factors if there is none) */
void pardiso_set_matrix(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
{
  if(nk_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
	      inputformat,jq,nzs3)||!nk_valid){
    nk_numeric(neq);
  }else{
    nk_stale=1;
  }
}

#endif
