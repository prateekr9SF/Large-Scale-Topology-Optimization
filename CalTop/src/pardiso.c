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
/* char envMKL[32];   moved to pardiso.h */

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

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
	{

		#ifdef PROFILING_ON
        	TAU_PROFILE_TIMER(t_pardiso_factor,"PARDISO: Factorization","",TAU_USER);
			TAU_PROFILE_START(t_pardiso_factor);
    	#endif

  		char *env;
		/*  char env1[32]; */
  		ITG i,j,k,l,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
      		msglvl=0,error=0,*irowpardiso=NULL,kflag,kstart,n,ifortran,
      		lfortran,index,id,k2;
  		ITG ndim,nthread,nthread_v;

  		double *b=NULL,*x=NULL;


  		iparm[0]=0;
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

  		for(i=0;i<64;i++)
		{
			pt[i]=0;
		}

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


  		FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
			   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

		#ifdef PROFILING_ON
        	TAU_PROFILE_STOP(t_pardiso_factor);
    	#endif
				   
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

	#ifdef PROFILING_ON
        TAU_PROFILE_TIMER(t_pardiso_solve,"PARDISO: Solve","",TAU_USER);
		TAU_PROFILE_START(t_pardiso_solve);
    #endif

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
    	msglvl=0,i,error=0;
  double *x=NULL;

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

	#ifdef PROFILING_ON
    	TAU_PROFILE_STOP(t_pardiso_solve);
	#endif

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

	#ifdef PROFILING_ON
        TAU_PROFILE_TIMER(t_pardiso_cleanup,
                          "PARDISO: Cleanup",
                          "",
                          TAU_USER);

        TAU_PROFILE_START(t_pardiso_cleanup);
    #endif

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

    #ifdef PROFILING_ON
        TAU_PROFILE_STOP(t_pardiso_cleanup);
    #endif

  return;
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

		#ifdef PROFILING_ON
        	TAU_PROFILE_TIMER(t_pardiso_total,"PARDISO: Total","",TAU_USER);
        	TAU_PROFILE_START(t_pardiso_total);
    	#endif

  		if(*neq==0)
		{
			#ifdef PROFILING_ON
        		TAU_PROFILE_STOP(t_pardiso_total);
    		#endif
			return;
		}
  		pardiso_factor(ad,au,adb,aub,sigma,icol,irow, 
			neq,nzs,symmetryflag,inputformat,jq,nzs3);

  		pardiso_solve(b,neq,symmetryflag,nrhs);

  		pardiso_cleanup(neq,symmetryflag);

		#ifdef PROFILING_ON
        	TAU_PROFILE_STOP(t_pardiso_total);
    	#endif

  		return;
	}

#endif

