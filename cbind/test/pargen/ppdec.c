/*----------------------------------------------------------------------------------*/
/*              Parallel Sparse BLAS  v2.2					    */
/*    (C) Copyright 2007 Salvatore Filippone    University of Rome Tor Vergata      */
/* 										    */
/*  Redistribution and use in source and binary forms, with or without		    */
/*  modification, are permitted provided that the following conditions		    */
/*  are met:									    */
/*    1. Redistributions of source code must retain the above copyright		    */
/*       notice, this list of conditions and the following disclaimer.		    */
/*    2. Redistributions in binary form must reproduce the above copyright	    */
/*       notice, this list of conditions, and the following disclaimer in the	    */
/*       documentation and/or other materials provided with the distribution.	    */
/*    3. The name of the PSBLAS group or the names of its contributors may	    */
/*       not be used to endorse or promote products derived from this		    */
/*       software without specific written permission.				    */
/* 										    */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS		    */
/*  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED	    */
/*  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR	    */
/*  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS  */
/*  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR	    */
/*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF	    */
/*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS	    */
/*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN	    */
/*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)	    */
/*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE	    */
/*  POSSIBILITY OF SUCH DAMAGE.							    */
/* 										    */
/*  										    */
/* File: ppdec.c   								    */
/*										    */
/* Program: ppdec 								    */
/* This sample program shows how to build and solve a sparse linear		    */
/*										    */
/* The program  solves a linear system based on the partial differential	    */
/* equation 									    */
/*										    */
/* 										    */
/*										    */
/* The equation generated is							    */
/*										    */
/*   b1 d d (u)  b2  d d (u)    a1 d (u))  a2 d (u)))   			    */
/* -   ------   -    ------  +    ----- +  ------     + a3 u = 0		    */
/*      dx dx         dy dy         dx        dy        			    */
/*										    */
/* 										    */
/* with  Dirichlet boundary conditions on the unit cube 			    */
/*										    */
/*    0<=x,y,z<=1								    */
/* 										    */
/* The equation is discretized with finite differences and uniform stepsize;	    */
/* the resulting  discrete  equation is						    */
/*										    */
/* ( u(x,y,z)(2b1+2b2+a1+a2)+u(x-1,y)(-b1-a1)+u(x,y-1)(-b2-a2)+			    */
/*  -u(x+1,y)b1-u(x,y+1)b2)*(1/h**2)						    */
/*										    */
/* Example adapted from: C.T.Kelley						    */
/*    Iterative Methods for Linear and Nonlinear Equations			    */
/*    SIAM 1995									    */
/*										    */
/*										    */
/* In this sample program the index space of the discretized			    */
/* computational domain is first numbered sequentially in a standard way, 	    */
/* then the corresponding vector is distributed according to an HPF BLOCK	    */
/* distribution directive.							    */
/*										    */
/* Boundary conditions are set in a very simple way, by adding 			    */
/* equations of the form							    */
/*										    */
/*   u(x,y) = rhs(x,y)								    */
/*                                                                                  */
/*----------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"

#define LINEBUFSIZE 1024
#define NBMAX       20

double  a1(double x, double y, double  z)
{
  return(1.0);
}
double a2(double x, double y, double  z)
{
  return(20.0*y);
}
double a3(double x, double y, double  z)
{
  return(1.0);
}
double  a4(double x, double y, double  z)
{
  return(1.0);
}
double  b1(double x, double y, double  z)
{
  return(1.0);
}
double b2(double x, double y, double  z)
{
  return(1.0);
}
double b3(double x, double y, double  z)
{
  return(1.0);
}

int matgen(int ictxt, int ng,int idim,int vg[],psb_c_dspmat *ah,psb_c_descriptor *cdh,
	   psb_c_dvector *xh, psb_c_dvector *bh, psb_c_dvector *rh)
{
  int iam, np;
  int x, y, z, el,glob_row,i,info,ret;
  double gx, gy, gz, deltah;
  double val[10*NBMAX], zt[NBMAX];
  int irow[10*NBMAX], icol[10*NBMAX];

  info = 0;
  psb_c_info(ictxt,&iam,&np);
  deltah = (double) 1.0/(idim-1);
  psb_c_set_index_base(1);
  for (glob_row=1; glob_row<=ng; glob_row++) {

    /* Check if I have to do something about this entry */
    if (vg[glob_row-1] == iam) {
      el=0;
      if ( (glob_row%(idim*idim)) == 0) {
	x = glob_row/(idim*idim);
      } else {
	x = glob_row/(idim*idim)+1;
      }
      if (((glob_row-(x-1)*idim*idim)%idim) == 0) {
	y = (glob_row-(x-1)*idim*idim)/idim;
      } else {
	y = (glob_row-(x-1)*idim*idim)/idim+1;
      }
      z = glob_row-(x-1)*idim*idim-(y-1)*idim;
      gx=x*deltah;
      gy=y*deltah;
      gz=z*deltah;
      zt[0] = 0.0;
      /*  internal point: build discretization */
      /*  term depending on   (x-1,y,z)        */
	
      if (x==1) {
	val[el] = -b1(gx,gy,gz)-a1(gx,gy,gz);
	val[el] /= deltah*deltah;
	zt[0] = exp(-gy*gy-gz*gz)*(-val[el]);
      } else {
	val[el]=-b1(gx,gy,gz) -a1(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x-2)*idim*idim+(y-1)*idim+(z);
	el=el+1;
      }
      /*  term depending on     (x,y-1,z) */ 
      if (y==1) { 
	val[el]=-b2(gx,gy,gz)-a2(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	zt[0] = exp(-gy*gy-gz*gz)*exp(-gx)*(-val[el]);
      } else {
	val[el]=-b2(gx,gy,gz)-a2(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x-1)*idim*idim+(y-2)*idim+(z);
	el=el+1;
      }
      /* term depending on     (x,y,z-1)*/
      if (z==1) { 
	val[el]=-b3(gx,gy,gz)-a3(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	zt[0] = exp(-gy*gy-gz*gz)*exp(-gx)*(-val[el]);
      } else {
	val[el]=-b3(gx,gy,gz)-a3(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x-1)*idim*idim+(y-1)*idim+(z-1);
	el=el+1;
      }
      /* term depending on     (x,y,z)*/
      val[el]=2*b1(gx,gy,gz)+2*b2(gx,gy,gz)+2*b3(gx,gy,gz)      
	+a1(gx,gy,gz)+a2(gx,gy,gz)+a3(gx,gy,gz);
      val[el] = val[el]/(deltah*deltah);
      icol[el]=(x-1)*idim*idim+(y-1)*idim+(z);
      el=el+1;
      /*  term depending on     (x,y,z+1) */
      if (z==idim) { 
	val[el]=-b1(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	zt[0] = exp(-gy*gy-gz*gz)*exp(-gx)*(-val[el]);
      } else {
	val[el]=-b1(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x-1)*idim*idim+(y-1)*idim+(z+1);
	el=el+1;
      }
      /* term depending on     (x,y+1,z) */
      if (y==idim) { 
	val[el]=-b2(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	zt[0] = exp(-gy*gy-gz*gz)*exp(-gx)*(-val[el]);
      } else {
	val[el]=-b2(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x-1)*idim*idim+(y)*idim+(z);
	el=el+1;
      }
      /*  term depending on     (x+1,y,z) */
      if (x<idim) { 
	val[el]=-b3(gx,gy,gz);
	val[el] = val[el]/(deltah*deltah);
	icol[el]=(x)*idim*idim+(y-1)*idim+(z);
	el=el+1;
      }
      if ((ret=psb_c_dspins(el,irow,icol,val,ah,cdh))!=0) 
	fprintf(stderr,"From psb_c_dspins: %d\n",ret); 
      psb_c_dgeins(1,irow,zt,bh,cdh);
      zt[0]=0.0;
      psb_c_dgeins(1,irow,zt,xh,cdh);
    }
  }
  
  if ((info=psb_c_cdasb(cdh))!=0)  return(info);
  
  if ((info=psb_c_dspasb(ah,cdh))!=0)  return(info); 
  
  if ((info=psb_c_dgeasb(xh,cdh))!=0)  return(info);
  if ((info=psb_c_dgeasb(bh,cdh))!=0)  return(info);
  if ((info=psb_c_dgeasb(rh,cdh))!=0)  return(info);
  return(info);

}
  
int main(int argc, char *argv[])
{
  int ictxt, iam, np;
  char methd[40], ptype[20], afmt[8], buffer[LINEBUFSIZE+1];
  int nparms;
  int idim,info,istop,itmax,itrace,irst,i,iter,ret;
  psb_c_dprec *ph;
  psb_c_dspmat *ah;
  psb_c_dvector *bh, *xh, *rh; 
  int *vg, ng, nb,nlr; 
  double t1,t2,eps,err;
  double *xv, *bv, *rv;
  double one=1.0, zero=0.0, res2;
  /* psb_c_SolverOptions options; */
  psb_c_descriptor *cdh;
  
  ictxt = psb_c_init();
  psb_c_info(ictxt,&iam,&np);
  fprintf(stdout,"Initialization: am %d of %d\n",iam,np);

  fflush(stdout);
  psb_c_barrier(ictxt);
  if (iam == 0) {
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d ",&nparms);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%s",methd);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%s",ptype);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%s",afmt);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d",&idim);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d",&istop);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d",&itmax);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d",&itrace);
    fgets(buffer,LINEBUFSIZE,stdin);
    sscanf(buffer,"%d",&irst);
  } 
  /* Now broadcast the values, and check they're OK */
  psb_c_ibcast(ictxt,1,&nparms,0);
  psb_c_hbcast(ictxt,methd,0);
  psb_c_hbcast(ictxt,ptype,0);
  psb_c_hbcast(ictxt,afmt,0);
  psb_c_ibcast(ictxt,1,&idim,0);
  psb_c_ibcast(ictxt,1,&istop,0);
  psb_c_ibcast(ictxt,1,&itmax,0);
  psb_c_ibcast(ictxt,1,&itrace,0);
  psb_c_ibcast(ictxt,1,&irst,0);
  
  fprintf(stderr,"%d Check on received: methd %s ptype %s afmt %s\n",
	  iam,methd,ptype,afmt);

  psb_c_barrier(ictxt);

  cdh=psb_c_new_descriptor();

  /* Simple minded BLOCK data distribution */ 
  ng = idim*idim*idim;
  nb = (ng+np-1)/np;
  if ((vg=malloc(ng*sizeof(int)))==NULL) {
    fprintf(stderr,"On %d: malloc failure\n",iam);
    psb_c_abort(ictxt);
  }
  for (i=0; i<ng; i++) {
    vg[i] = i/nb; 
  }
  if ((info=psb_c_cdall_vg(ng,vg,ictxt,cdh))!=0) {
    fprintf(stderr,"From cdall: %d\nBailing out\n",info);
    psb_c_abort(ictxt);
  }

  bh  = psb_c_new_dvector();
  xh  = psb_c_new_dvector();
  rh  = psb_c_new_dvector();
  ah  = psb_c_new_dspmat();
  fprintf(stderr,"From psb_c_new_dspmat: %p\n",ah); 

  /* Allocate mem space for sparse matrix and vectors */ 
  ret=psb_c_dspall(ah,cdh); 
  fprintf(stderr,"From psb_c_dspall: %d\n",ret); 
  psb_c_dgeall(bh,cdh);
  psb_c_dgeall(xh,cdh);
  psb_c_dgeall(rh,cdh);
  

  
  /* Matrix generation  */
  if (matgen(ictxt,ng,idim,vg,ah,cdh,xh,bh,rh) != 0) {
    fprintf(stderr,"Error during matrix build loop\n");
    psb_c_abort(ictxt);
  }    
  psb_c_dmat_name_print(ah,"cbindmat.mtx");
  psb_c_barrier(ictxt);
  /* Set up the preconditioner */ 
  ph  = psb_c_new_dprec();
  psb_c_dprecinit(ph,ptype);
  ret=psb_c_dprecbld(ah,cdh,ph);
  fprintf(stderr,"From psb_c_dprecbld: %d\n",ret); 
#if 0 
  /* Set up the solver options */ 
  psb_c_DefaultSolverOptions(&options);
  options.eps    = 1.e-9;
  options.itmax  = itmax;
  options.irst   = irst;
  options.itrace = 1;
  options.istop  = istop;
  psb_c_seterraction_ret();
  t1=psb_c_wtime();
  ret=psb_c_dkrylov(methd,ah,ph,bh,xh,cdh,&options);
  t2=psb_c_wtime();
  iter = options.iter;
  err  = options.err;
  fprintf(stderr,"From krylov: %d %lf, %d %d\n",iter,err,ret,psb_c_get_errstatus());
  if (psb_c_get_errstatus() != 0) {
    psb_c_print_errmsg();
  }
  fprintf(stderr,"After cleanup %d\n",psb_c_get_errstatus());
#endif   
  /* Check 2-norm of residual on exit */ 
  psb_c_dgeaxpby(one,bh,zero,rh,cdh); 
  psb_c_dspmm(-one,ah,xh,one,rh,cdh); 
  res2=psb_c_dgenrm2(rh,cdh);

  if (iam==0) {
    fprintf(stdout,"Time: %lf\n",(t2-t1));
    fprintf(stdout,"Iter: %d\n",iter);
    fprintf(stdout,"Err: %lg\n",err);
    fprintf(stdout,"||r||_2: %lg\n",res2);
    
  }

#if 0
  bv = psb_c_dvect_get_cpy(bh);
  nlr=psb_c_cd_get_local_rows(cdh);
  for (i=0;i<nlr; i++)
    fprintf(stdout,"RHS: %d %d %lf\n",iam,i,bv[i]);


  xv = psb_c_dvect_get_cpy(xh);
  nlr=psb_c_cd_get_local_rows(cdh);
  for (i=0;i<nlr; i++)
    fprintf(stdout,"SOL: %d %d %lf\n",iam,i,xv[i]);


  rv = psb_c_dvect_get_cpy(rh);
  nlr=psb_c_cd_get_local_rows(cdh);
  for (i=0;i<nlr; i++)
    fprintf(stdout,"RES: %d %d %lf\n",iam,i,rv[i]);

#endif
  
  /* Clean up memory */ 
  if ((info=psb_c_dgefree(xh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(ictxt);
  }
  if ((info=psb_c_dgefree(bh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(ictxt);
  }
  if ((info=psb_c_dgefree(rh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(ictxt);
  }

  if ((info=psb_c_cdfree(cdh))!=0) {
    fprintf(stderr,"From cdfree: %d\nBailing out\n",info);
    psb_c_abort(ictxt);
  }
  fprintf(stderr,"pointer  from cdfree: %p\n",cdh->descriptor);
 
  /* Clean up object handles  */   
  free(ph);
  free(xh);
  free(bh);
  free(ah);
  free(cdh);


  fprintf(stderr,"program completed successfully\n");

  psb_c_barrier(ictxt);
  psb_c_exit(ictxt);
}
