/*----------------------------------------------------------------------------------*/
/*              Parallel Sparse BLAS  v 3.5.0					    */
/*    (C) Copyright 2017 Salvatore Filippone                    */
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
/* distribution directive. The discretization ensures there are IDIM                */
/* *internal* points in each direction.                                             */
/*                                                                                  */
/*----------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"
#include "psb_krylov_cbind.h"

#define LINEBUFSIZE 1024
#define NBMAX       20
#define DUMPMATRIX  0

double  a1(double x, double y, double  z)
{
  return(1.0/80.0);
}
double a2(double x, double y, double  z)
{
  return(1.0/80.0);
}
double a3(double x, double y, double  z)
{
  return(1.0/80.0);
}
double  c(double x, double y, double  z)
{
  return(0.0);
}
double  b1(double x, double y, double  z)
{
  return(0.0/sqrt(3.0));
}
double b2(double x, double y, double  z)
{
  return(0.0/sqrt(3.0));
}
double b3(double x, double y, double  z)
{
  return(0.0/sqrt(3.0));
}

double g(double x, double y, double z)
{
  if (x == 1.0) {
    return(1.0);
  } else if (x == 0.0) {
    return( exp(-y*y-z*z));
  } else {
    return(0.0);
  }
}

psb_i_t matgen(psb_c_ctxt cctxt, psb_i_t nl, psb_i_t idim, psb_l_t vl[],
	       psb_c_dspmat *ah,psb_c_descriptor *cdh,
	       psb_c_dvector *xh, psb_c_dvector *bh, psb_c_dvector *rh)
{
  psb_i_t iam, np;
  psb_l_t ix, iy, iz, el,glob_row;
  psb_i_t i, k, info,ret;
  double x, y, z, deltah, sqdeltah, deltah2;
  double val[10*NBMAX], zt[NBMAX];
  psb_l_t irow[10*NBMAX], icol[10*NBMAX];

  info = 0;
  psb_c_info(cctxt,&iam,&np);
  deltah = (double) 1.0/(idim+1);
  sqdeltah = deltah*deltah;
  deltah2  = 2.0* deltah;
  psb_c_set_index_base(0);
  for (i=0; i<nl;  i++) {
    glob_row=vl[i];
    //if ((i%100000 == 0)||(i<10)) fprintf(stderr,"%d: generation loop at %d %ld \n",iam,i,glob_row);
    el=0;
    ix = glob_row/(idim*idim);
    iy = (glob_row-ix*idim*idim)/idim;
    iz = glob_row-ix*idim*idim-iy*idim;
    x=(ix+1)*deltah;
    y=(iy+1)*deltah;
    z=(iz+1)*deltah;
    zt[0] = 0.0;
    /*  internal point: build discretization */
    /*  term depending on   (x-1,y,z)        */
    val[el] = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2;	
    if (ix==0) {
      zt[0] += g(0.0,y,z)*(-val[el]);
    } else {
      icol[el]=(ix-1)*idim*idim+(iy)*idim+(iz);
      el=el+1;
    }
    /*  term depending on     (x,y-1,z) */
    val[el]  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2;	      
    if (iy==0) { 
      zt[0] += g(x,0.0,z)*(-val[el]);
    } else {
      icol[el]=(ix)*idim*idim+(iy-1)*idim+(iz);
      el=el+1;
    }
    /* term depending on     (x,y,z-1)*/
    val[el]=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2;
    if (iz==0) { 
      zt[0] += g(x,y,0.0)*(-val[el]);
    } else {
      icol[el]=(ix)*idim*idim+(iy)*idim+(iz-1);
      el=el+1;
    }
    /* term depending on     (x,y,z)*/
    val[el]=2.0*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah + c(x,y,z);
    icol[el]=(ix)*idim*idim+(iy)*idim+(iz);
    el=el+1;
    /*  term depending on     (x,y,z+1) */
    val[el] = -a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2;
    if (iz==idim-1) { 
      zt[0] += g(x,y,1.0)*(-val[el]);
    } else {
      icol[el]=(ix)*idim*idim+(iy)*idim+(iz+1);
      el=el+1;
    }
    /* term depending on     (x,y+1,z) */
    val[el] = -a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2;
    if (iy==idim-1) { 
      zt[0] += g(x,1.0,z)*(-val[el]);
    } else {
      icol[el]=(ix)*idim*idim+(iy+1)*idim+(iz);
      el=el+1;
    }
    /*  term depending on     (x+1,y,z) */
    val[el] = -a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2;
    if (ix==idim-1) { 
      zt[0] += g(1.0,y,z)*(-val[el]);
    } else {
      icol[el]=(ix+1)*idim*idim+(iy)*idim+(iz);
      el=el+1;
    }
    for (k=0; k<el; k++) irow[k]=glob_row;
    if ((ret=psb_c_dspins(el,irow,icol,val,ah,cdh))!=0) 
      fprintf(stderr,"From psb_c_dspins: %d\n",ret);
    irow[0] = glob_row; 
    psb_c_dgeins(1,irow,zt,bh,cdh);
    zt[0]=0.0;
    psb_c_dgeins(1,irow,zt,xh,cdh);   
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
  psb_c_ctxt *cctxt;
  psb_i_t iam, np;
  char methd[40], ptype[20], afmt[8], buffer[LINEBUFSIZE+1];
  psb_i_t nparms;
  psb_i_t idim,info,istop,itmax,itrace,irst,iter,ret;
  psb_c_dprec *ph;
  psb_c_dspmat *ah;
  psb_c_dvector *bh, *xh, *rh; 
  psb_i_t nb,nlr, nl;
  psb_l_t i,ng, *vl, k;
  double t1,t2,eps,err;
  double *xv, *bv, *rv;
  double one=1.0, zero=0.0, res2;
  psb_c_SolverOptions options; 
  psb_c_descriptor *cdh;
  FILE *vectfile;

  cctxt = psb_c_new_ctxt();
  psb_c_init(cctxt);
  psb_c_info(*cctxt,&iam,&np);

  psb_c_barrier(*cctxt);
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
  psb_c_ibcast(*cctxt,1,&nparms,0);
  psb_c_hbcast(*cctxt,methd,0);
  psb_c_hbcast(*cctxt,ptype,0);
  psb_c_hbcast(*cctxt,afmt,0);
  psb_c_ibcast(*cctxt,1,&idim,0);
  psb_c_ibcast(*cctxt,1,&istop,0);
  psb_c_ibcast(*cctxt,1,&itmax,0);
  psb_c_ibcast(*cctxt,1,&itrace,0);
  psb_c_ibcast(*cctxt,1,&irst,0);
  
  fflush(stderr);
  psb_c_barrier(*cctxt);

  cdh=psb_c_new_descriptor();
  psb_c_set_index_base(0);

  /* Simple minded BLOCK data distribution */ 
  ng = ((psb_l_t) idim)*idim*idim;  
  nb = (ng+np-1)/np;
  nl = nb;
  if ( (ng -iam*nb) < nl) nl = ng -iam*nb; 
    if ((vl=malloc(nb*sizeof(psb_l_t)))==NULL) {
    fprintf(stderr,"On %d: malloc failure\n",iam);
    psb_c_abort(*cctxt);
  }
  i = ((psb_l_t)iam) * nb;
  for (k=0; k<nl; k++)
    vl[k] = i+k; 

  if ((info=psb_c_cdall_vl(nl,vl,*cctxt,cdh))!=0) {
    fprintf(stderr,"From cdall: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }

  bh  = psb_c_new_dvector();
  xh  = psb_c_new_dvector();
  rh  = psb_c_new_dvector();
  ah  = psb_c_new_dspmat();
  //fprintf(stderr,"From psb_c_new_dspmat: %p\n",ah); 

  /* Allocate mem space for sparse matrix and vectors */ 
  ret=psb_c_dspall(ah,cdh); 
  //fprintf(stderr,"From psb_c_dspall: %d\n",ret); 
  psb_c_dgeall(bh,cdh);
  psb_c_dgeall(xh,cdh);
  psb_c_dgeall(rh,cdh);
  

  
  /* Matrix generation  */
  if (matgen(*cctxt,nl,idim,vl,ah,cdh,xh,bh,rh) != 0) {
    fprintf(stderr,"Error during matrix build loop\n");
    psb_c_abort(*cctxt);
  }    
  psb_c_barrier(*cctxt);
  /* Set up the preconditioner */ 
  ph  = psb_c_new_dprec();
  psb_c_dprecinit(*cctxt,ph,ptype);
  ret=psb_c_dprecbld(ah,cdh,ph);
  //fprintf(stderr,"From psb_c_dprecbld: %d\n",ret); 

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
  //fprintf(stderr,"From krylov: %d %lf, %d %d\n",iter,err,ret,psb_c_get_errstatus());
  if (psb_c_get_errstatus() != 0) {
    psb_c_print_errmsg();
  }
  //fprintf(stderr,"After cleanup %d\n",psb_c_get_errstatus());
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

#if DUMPATRIX
  psb_c_dmat_name_print(ah,"cbindmat.mtx");
  nlr = psb_c_cd_get_local_rows(cdh);
  bv  = psb_c_dvect_get_cpy(bh);
  vectfile=fopen("cbindb.mtx","w");
  for (i=0;i<nlr; i++)
    fprintf(vectfile,"%lf\n",bv[i]);
  fclose(vectfile);


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
  if ((info=psb_c_dprecfree(ph))!=0) {
    fprintf(stderr,"From dprecfree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }
  if ((info=psb_c_dgefree(xh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }
  if ((info=psb_c_dgefree(bh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }
  if ((info=psb_c_dgefree(rh,cdh))!=0) {
    fprintf(stderr,"From dgefree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }
  if ((info=psb_c_dspfree(ah,cdh))!=0) {
    fprintf(stderr,"From dspfree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }

  if ((info=psb_c_cdfree(cdh))!=0) {
    fprintf(stderr,"From cdfree: %d\nBailing out\n",info);
    psb_c_abort(*cctxt);
  }
  //fprintf(stderr,"pointer  from cdfree: %p\n",cdh->descriptor);
 
  /* Clean up object handles  */   
  free(ph);
  free(xh);
  free(bh);
  free(rh);
  free(ah);
  free(cdh);

  /* further cleanup */
  free(vl);
  //if (iam == 0) fprintf(stderr,"program completed successfully\n");

  psb_c_barrier(*cctxt);
  psb_c_exit(*cctxt);
  free(cctxt);
}
