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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "main.h"

void readinput(char *jobnamec, char **inpcp, ITG *nline, ITG *nset,
	       ITG *ipoinp, ITG **inpp, ITG **ipoinpcp, ITG *ithermal,
               ITG *nuel){

  /*   reads and stores the input deck in inpcp; determines the
       number of sets  */

  /* Array to manage file pointers for input files, allows handling multiple files (e.g., when files are nested) */
  FILE *f1[10];

  char buff[1320]="", fninp[132]="", includefn[132]="", *inpc=NULL,
       textpart[2112]="",*set=NULL;

  ITG i,j,k,n,in=0,nlinemax=100000,irestartread,irestartstep,
      icntrl,nload,nforc,nboun,nk,ne,nmpc,nalset,nmat,ntmat,npmat,
      norien,nam,nprint,mi[3],ntrans,ncs,namtot,ncmat,memmpc,ne1d,
      ne2d,nflow,*meminset=NULL,*rmeminset=NULL, *inp=NULL,ntie,
      nener,nstate,nentries=17,ifreeinp,ikey,lincludefn,nslavs,
      nbody,ncharmax=1000000,*ipoinpc=NULL,ichangefriction=0,nkon,
      ifile,mcs,initialtemperature=0,nprop,mortar,ifacecount,
      nintpoint,infree[4],iheading=0,ichangesurfacebehavior=0,
      nef; 

  /* initialization */

  /* nentries is the number of different keyword cards for which
     the input deck order is important, cf keystart.f */

      
  /* Allocate memory for storing input commands and data dynamically */
  NNEW(inpc,char,ncharmax);
  NNEW(ipoinpc,ITG,nlinemax+1);
  NNEW(inp,ITG,3*nlinemax);
  
  // Initialize the line counter
  *nline=0;
  
   /* Initialize input order array */
  for(int i=0;i<2*nentries;i++)
  {
    ipoinp[i]=0;
  }
  
  ifreeinp=1;
  ikey=0;

  /* Initialize the input file handling */
  printf("Parsing input file... \n");
  strcpy(fninp,jobnamec);// Copy job name to file name buffer
  strcat(fninp,".inp");  // Append '.inp' extension to form the full input file name
  
  /* Open the input file and check for errors */
  if((f1[in]=fopen(fninp,"r"))==NULL)
  {
    printf("*ERROR in readinput: cannot open file %s\n",fninp);
    exit(0); // Exit program if file cannot be opened
  }

/* Parse each line in the .inp file and store the order in which the 
    input is to be read in the fields ipoinp and inp  */


/* Read from the input file until the end */
  do
  {
    if(fgets(buff,1320,f1[in])==NULL) // Read a line from the file
    {
	    fclose(f1[in]); // Close the current file if nothing to be read
	    if(in!=0)
      {
	      in--; // If there are nested files, move to the previous file
	      continue; // Continue reading from the previous file
	    }
	  else
      {
        break;
      } // If no more files, exit the loop
    }

    /* check for heading lines that start from "*": should not be changed */
    if(iheading==1)
    {
	    if((buff[0]=='*')&&(buff[1]!='*'))
      {
	      iheading=0;
	    }
    }

    /* storing the significant characters */
    /* get rid of blanks  */
    k=0;
    i=-1;
    if(iheading==0)
    {
	    do
      {
	      i++;
	      if((buff[i]=='\0')||(buff[i]=='\n')||(buff[i]=='\r')||(k==1320)) break;
	      if((buff[i]==' ')||(buff[i]=='\t')) continue;
	      buff[k]=buff[i];
	      k++;
	    }
      while(1);
    }
    else
    {
	    do
      {
	      i++;
	      if((buff[i]=='\0')||(buff[i]=='\n')||(buff[i]=='\r')||(k==1320)) break;
	      buff[k]=buff[i];
	      k++;
	    }
      while(1);
    }
	
    /* check for blank lines and comments */
    if(k==0) continue;
      if(strcmp1(&buff[0],"**")==0) continue;

      /* changing to uppercase except filenames */
      if(iheading==0)
      {
	      j=0;
	      ifile=0;
	      do
        {
	        if(j>=6)
          {
		        if(strcmp1(&buff[j-6],"INPUT=")==0) ifile=1;
	        }
	        if(j>=7)
          {
		        if(strcmp1(&buff[j-7],"OUTPUT=")==0) ifile=1;
	        }
	        if(j>=9)
          {
		        if(strcmp1(&buff[j-9],"FILENAME=")==0) ifile=1;
	        }
	        if(ifile==1)
          {
		        do
            {
		          if(strcmp1(&buff[j],",")!=0)
              {
			          j++;
		          }
              else
              {
			          ifile=0;
			          break;
		          }
		        }
            while(j<k);
	        }
          else
          {
		        buff[j]=toupper(buff[j]);
	        }
	        j++;
	      }
        while(j<k);
      }
      /* check for a *HEADING card */
      if(strcmp1(&buff[0],"*HEADING")==0)
      {
	      iheading=1;
      }
      /* check for a *KINEMATIC or *DISTRIBUTING card and change
         the asterisk into a C (is a "dependent" card of the 
         *COUPLING card */

      if((strcmp1(&buff[0],"*KINEMATIC")==0)||
         ((strcmp1(&buff[0],"*DISTRIBUTING")==0)&&
          (strcmp1(&buff[0],"*DISTRIBUTINGCOUPLING")!=0)))
        {
	        buff[0]='C';
        }
	  
      /* check for include statements */
	  
      if(strcmp1(&buff[0],"*INCLUDE")==0)
      {
        printf("Reading INPUT card...");
	      lincludefn=k;
	      FORTRAN(includefilename,(buff,includefn,&lincludefn));
        includefn[lincludefn]='\0';
	      in++;

	      if(in>9)
        {
	        printf("*ERROR in readinput: include statements can \n not be cascaded over more than 9 levels\n");
	      }

	      if((f1[in]=fopen(includefn,"r"))==NULL)
        {
	        printf("*ERROR in readinput: cannot open file %s\n",includefn);
	        exit(0);
	      }

        if((f1[in]=fopen(includefn,"r")) != NULL)
        {
	        printf("Include file:  %s\n",includefn);
	      }

        continue;
      }

      /* adding a line */
	    (*nline)++;

      if(*nline>nlinemax)
      {
	      nlinemax=(ITG)(1.1*nlinemax);
	      RENEW(ipoinpc,ITG,nlinemax+1);
	      RENEW(inp,ITG,3*nlinemax);
      }

      /* checking the total number of characters */

      if(ipoinpc[*nline-1]+k>ncharmax)
      {
	      ncharmax=(ITG)(1.1*ncharmax);
	      RENEW(inpc,char,ncharmax);
      }
	  
      /* copying into inpc */

      for(j=0;j<k;j++)
      {
	      inpc[ipoinpc[*nline-1]+j]=buff[j];
      }

      ipoinpc[*nline]=ipoinpc[*nline-1]+k;
      /* counting sets */
      
      if(strcmp1(&buff[0],"*AMPLITUDE")==0)
      {
        printf("ERROR: AMPLITUDE NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"AMPLITUDE",
                          nline,&ikey));                         
      }
      else if(strcmp1(&buff[0],"*CHANGEFRICTION")==0)
      {
        printf("ERROR: CHANGEFRICTION NOT SUPPORTED! \n");
        exit(0); 
	      ichangefriction=1;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CHANGESURFACEBEHAVIOR")==0)
      {
        printf("ERROR: CHANGESURFACEBEHAVIOR NOT SUPPORTED! \n");
        exit(0); 
	      ichangesurfacebehavior=1;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONDUCTIVITY")==0)
      {
        printf("ERROR: CONDUCTIVITY NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONTACTDAMPING")==0)
      {
        printf("ERROR: CONTACTDAMPING NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CONTACTPAIR")==0)
      {
        printf("ERROR: CONTACTPAIR NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"CONTACTPAIR",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*COUPLING")==0)
      {
        printf("ERROR: COUPLING NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"COUPLING",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CREEP")==0)
      {
        printf("ERROR: CREEP NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*CYCLICHARDENING")==0)
      {
        printf("ERROR: CYCLICHARDENING NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DAMPING")==0)
      {
        printf("ERROR: DAMPING NOT SUPPORTED! \n");
        exit(0); 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELASTIC")==0)
      { 
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DEFORMATIONPLASTICITY")==0)
      {
        printf("ERROR: DEFORMATIONPLASTICITY NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DENSITY")==0)
      {
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*DEPVAR")==0)
      {
        printf("ERROR: DEPVAR NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELASTIC")==0)
      {
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ELECTRICALCONDUCTIVITY")==0)
      {
        printf("ERROR: ELECTRICALCONDUCTIVITY NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*ELEMENT")==0)&&
              (strcmp1(&buff[0],"*ELEMENTOUTPUT")!=0))
              {
                (*nset)++;
                FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ELEMENT",
                          nline,&ikey));
              }
      else if(strcmp1(&buff[0],"*ELSET")==0)
      {
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ELSET",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*EXPANSION")==0)
      {
        printf("ERROR: EXPANSION NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*FLUIDCONSTANTS")==0)
      {
        printf("ERROR: FLUIDCONSTANTS NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*FRICTION")==0)&&(ichangefriction==0))
      {
        printf("ERROR: FRICTION NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*GAPCONDUCTANCE")==0)
      {
        printf("ERROR: GAPCONDUCTANCE NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*GAPHEATGENERATION")==0)
      {
        printf("ERROR: GAPHEATGENERATION NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*HYPERELASTIC")==0)
      {
        printf("ERROR: HYPERELASTIC NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*HYPERFOAM")==0)
      {
        printf("ERROR: HYPERFOAM NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*INITIALCONDITIONS")==0)
      {
        printf("ERROR:INNITIALCONDITIONS NOT SUPPORTED! \n");
        exit(0);
	      FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INITIALCONDITIONS",
			    nline,&ikey));

	      FORTRAN(splitline,(buff,textpart,&n));
	      for(i=0;i<n;i++)
        {
	        if(strcmp1(&textpart[(long long)132*i],"TYPE=TEMPERATURE")==0)
          {
		        initialtemperature=1;
	        }
        }
      }
      else if(strcmp1(&buff[0],"*MAGNETICPERMEABILITY")==0)
      {
        printf("ERROR:MAGNETICPERMEABILITY NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*MATERIAL")==0)
      {
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if((strcmp1(&buff[0],"*NODE")==0)&&
	      (strcmp1(&buff[0],"*NODEPRINT")!=0)&&
	      (strcmp1(&buff[0],"*NODEOUTPUT")!=0)&&
	      (strcmp1(&buff[0],"*NODEFILE")!=0))
        {
          (*nset)++;
          FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"NODE",
                          nline,&ikey));
        }
      else if(strcmp1(&buff[0],"*NSET")==0)
      {
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"NSET",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*ORIENTATION")==0)
      {
        printf("ERROR: ORIENTATION NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"ORIENTATION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*PLASTIC")==0)
      {
        printf("ERROR: PLASTIC NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*RESTART")==0)
      {
        printf("ERROR: RESTART NOT SUPPORTED! \n");
        exit(0);
	      irestartread=0;
	      irestartstep=0;
	      strcpy1(&buff[k]," ",1);
	      FORTRAN(splitline,(buff,textpart,&n));

	      for(i=0;i<n;i++)
        {
	        if(strcmp1(&textpart[(long long)132*i],"READ")==0)
          {
		        irestartread=1;
	        }
	        if(strcmp1(&textpart[(long long)132*i],"STEP")==0)
          {
		        irestartstep=atoi(&textpart[(long long)132*i+5]);
	        }
        }
          if(irestartread==1)
          {
            icntrl=0;
            FORTRAN(restartshort,(nset,&nload,&nbody,&nforc,&nboun,&nk,
              &ne,&nmpc,&nalset,&nmat,&ntmat,&npmat,&norien,&nam,
              &nprint,mi,&ntrans,&ncs,&namtot,&ncmat,&memmpc,
              &ne1d,&ne2d,&nflow,set,meminset,rmeminset,jobnamec,
	            &irestartstep,&icntrl,ithermal,&nener,&nstate,&ntie,
	            &nslavs,&nkon,&mcs,&nprop,&mortar,&ifacecount,&nintpoint,
	            infree,&nef));

            FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"RESTART,READ",
                              nline,&ikey));
	        }

          else
          {
            FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                              nline,&ikey));
          }
      } // Restart end
      else if(strcmp1(&buff[0],"*SPECIFICGASCONSTANT")==0)
      {
        printf("ERROR: SPECIFICGASCONSTANT NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SPECIFICHEAT")==0)
      {
        printf("ERROR: SPECIFICHEAT NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SUBMODEL")==0)
      { 
        printf("ERROR: SUBMODEL NOT SUPPORTED! \n");
        exit(0);

	      (*nset)+=2;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SURFACEINTERACTION")==0)
      {
        printf("ERROR: SURFACEINTEGRATION NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*SURFACEBEHAVIOR")==0)
      {
        printf("ERROR: SURFACEBEHAVIOR NOT SUPPORTED! \n");
        exit(0);
	      if(ichangesurfacebehavior==0)
        {
	        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"INTERACTION",
                          nline,&ikey));
	      }
        else
        {
	        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
				  nline,&ikey));
	      }
      }
      else if(strcmp1(&buff[0],"*SURFACE")==0)
      {
        (*nset)++;
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"SURFACE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*TIE")==0)
      {
        printf("ERROR: TIE NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"TIE",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*TRANSFORM")==0)
      {
        printf("ERROR: TRANSFORM NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"TRANSFORM",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*USERELEMENT")==0)
      {
        printf("ERROR: USERELEMENT NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"USERELEMENT",
                          nline,&ikey));
	    (*nuel)++;
      }
      else if(strcmp1(&buff[0],"*USERMATERIAL")==0)
      {
        printf("ERROR: USERMATERIAL NOT SUPPORTED! \n");
        exit(0);
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"MATERIAL",
                          nline,&ikey));
      }
      else if(strcmp1(&buff[0],"*")==0)
      {
        FORTRAN(keystart,(&ifreeinp,ipoinp,inp,"REST",
                          nline,&ikey));

        /* checking whether the calculation is mechanical,
           thermal or thermomechanical: needed to know
           which mpc's to apply to 2-D elements */

	      if((strcmp1(&buff[0],"*STATIC")==0)||
	      (strcmp1(&buff[0],"*VISCO")==0)||
	      (strcmp1(&buff[0],"*DYNAMIC")==0))
        {
	        if(ithermal[1]==0)
          {
		        if(initialtemperature==1)ithermal[1]=1;
	        }
          else if(ithermal[1]==2)
          {
		        ithermal[1]=3;
	        }
	      }
        else if(strcmp1(&buff[0],"*HEATTRANSFER")==0)
        {
          printf("ERROR: HEATTRANSFER NOT SUPPORTED! \n");
          exit(0);
	        if(ithermal[1]<2) ithermal[1]=ithermal[1]+2;
	      }
        else if(strcmp1(&buff[0],"*COUPLEDTEMPERATURE-DISPLACEMENT")==0)
        {
          printf("ERROR: COUPLEDTEMPERATURE-DISPLACEMENT NOT SUPPORTED! \n");
          exit(0);
	        ithermal[1]=3;
	      }
        else if(strcmp1(&buff[0],"*UNCOUPLEDTEMPERATURE-DISPLACEMENT")==0)
        {
          printf("ERROR: UNCOUPLEDTEMPERATURE-DISPLACEMENT NOT SUPPORTED! \n");
          exit(0);
	        ithermal[1]=3;
	      }
        else if(strcmp1(&buff[0],"*ELECTROMAGNETICS")==0)
        {
          printf("ERROR: ELECTROMAGNETICS NOT SUPPORTED! \n");
          exit(0);
	        ithermal[1]=3;
	      }
      }
    }
    while(1);
    inp[3*ipoinp[2*ikey-1]-2]=*nline;
    RENEW(inpc,char,(long long)132**nline);
    RENEW(inp,ITG,3*ipoinp[2*ikey-1]);
    *inpcp=inpc;
    *ipoinpcp=ipoinpc;
    *inpp=inp;

  
//  FORTRAN(writeinput,(inpc,ipoinp,inp,nline,&ipoinp[2*ikey-1],ipoinpc));

  return;

}






