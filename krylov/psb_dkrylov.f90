!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
! File: psb_krylov_mod.f90
!  Interfaces for Krylov subspace iterative methods.
!
!
! Subroutine: psb_dkrylov
! 
!    Front-end for the Krylov subspace iterations, realversion
!    
! Arguments:
!
!    methd  -  character                    The specific method; can take the values:
!                                           CG
!                                           CGS
!                                           BICG
!                                           BICGSTAB
!                                           BICGSTABL
!                                           RGMRES
!                                           
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
!    b      -  real,dimension(:)         Input: vector containing the
!                                           right hand side B
!    x      -  real,dimension(:)         Input/Output: vector containing the
!                                           initial guess and final solution X.
!    eps    -  real                         Input: Stopping tolerance; the iteration is
!                                           stopped when the error
!                                           estimate |err| <= eps
!                                           
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                     Output: Return code
!
!    itmax  -  integer(optional)            Input: maximum number of iterations to be
!                                           performed.
!    iter   -  integer(optional)            Output: how many iterations have been
!                                           performed.
!    err    -  real   (optional)            Output: error estimate on exit
!    itrace -  integer(optional)            Input: print an informational message
!                                           with the error estimate every itrace
!                                           iterations
!    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
!                                           BICGSTAB(L) methods
!    istop  -  integer(optional)            Input: stopping criterion, or how
!                                           to estimate the error. 
!                                           1: err =  |r|/(|a||x|+|b|)
!                                           2: err =  |r|/|b|
!                                           where r is the (preconditioned, recursive
!                                           estimate of) residual 
! 
Subroutine psb_dkrylov_vect(method,a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop,cond)

  use psb_base_mod
  use psb_prec_mod,only : psb_dprec_type
  use psb_krylov_mod, psb_protect_name => psb_dkrylov_vect

  character(len=*)                     :: method
  Type(psb_dspmat_type), Intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_dprec_type), intent(inout) :: prec 
  type(psb_d_vect_type), Intent(inout) :: b
  type(psb_d_vect_type), Intent(inout) :: x
  Real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_dpk_), Optional, Intent(out) :: err,cond


  abstract interface
    subroutine psb_dkryl_vect(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      import :: psb_ipk_, psb_dpk_, psb_desc_type, &
           & psb_dspmat_type, psb_dprec_type, psb_d_vect_type
      type(psb_dspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      type(psb_d_vect_type), Intent(inout) :: b
      type(psb_d_vect_type), Intent(inout) :: x
      real(psb_dpk_), intent(in)           :: eps
      class(psb_dprec_type), intent(inout) :: prec
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), optional, intent(in)        :: itmax, itrace,istop
      integer(psb_ipk_), optional, intent(out)       :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_dkryl_vect
    Subroutine psb_dkryl_rest_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      import :: psb_ipk_, psb_dpk_, psb_desc_type, &
           & psb_dspmat_type, psb_dprec_type, psb_d_vect_type
      Type(psb_dspmat_type), Intent(in)    :: a
      Type(psb_desc_type), Intent(in)      :: desc_a
      class(psb_dprec_type), intent(inout) :: prec
      type(psb_d_vect_type), Intent(inout) :: b
      type(psb_d_vect_type), Intent(inout) :: x
      Real(psb_dpk_), Intent(in)           :: eps
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
      integer(psb_ipk_), Optional, Intent(out)       :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_dkryl_rest_vect
    Subroutine psb_dkryl_cond_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,istop,cond)
      import :: psb_ipk_, psb_dpk_, psb_desc_type, &
           & psb_dspmat_type, psb_dprec_type, psb_d_vect_type
      Type(psb_dspmat_type), Intent(in)    :: a
      Type(psb_desc_type), Intent(in)      :: desc_a
      class(psb_dprec_type), intent(inout) :: prec
      type(psb_d_vect_type), Intent(inout) :: b
      type(psb_d_vect_type), Intent(inout) :: x
      Real(psb_dpk_), Intent(in)           :: eps
      integer(psb_ipk_), intent(out)                 :: info
      integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace,istop
      integer(psb_ipk_), Optional, Intent(out)       :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err, cond
    end subroutine psb_dkryl_cond_vect
  end interface

  procedure(psb_dkryl_vect) :: psb_dbicg_vect, psb_dcgstab_vect,&
       & psb_dcgs_vect
  procedure(psb_dkryl_rest_vect) :: psb_drgmres_vect, psb_dcgstabl_vect, psb_dgcr_vect
  procedure(psb_dkryl_cond_vect) :: psb_dcg_vect, psb_dfcg_vect

  logical           :: do_alloc_wrk
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: me,np,err_act, itrace_
  character(len=20)             :: name

  info = psb_success_
  name = 'psb_krylov'
  call psb_erractionsave(err_act)

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)

  ! Default return for COND
  if (present(cond)) cond = dzero

  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = -1
  end if

  do_alloc_wrk = .not.prec%is_allocated_wrk()
  if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v,desc=desc_a)

  select case(psb_toupper(method))
  case('CG') 
    call  psb_dcg_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop,cond=cond)
  case('FCG') 
    call  psb_dfcg_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop,cond=cond)
  case('GCR') 
    call  psb_dgcr_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop)
  case('CGS') 
    call  psb_dcgs_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop)
  case('BICG') 
    call  psb_dbicg_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop)
  case('BICGSTAB') 
    call  psb_dcgstab_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop)
  case('RGMRES','GMRES')
    call  psb_drgmres_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,irst=irst,istop=istop)
  case('BICGSTABL')
    call  psb_dcgstabl_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,irst=irst,istop=istop)
  case default
    if (me == 0) write(psb_err_unit,*) trim(name),&
         & ': Warning: Unknown method  ',method,&
         & ', defaulting to BiCGSTAB'
    call  psb_dcgstab_vect(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace=itrace_,istop=istop)
  end select

  if ((info==psb_success_).and.do_alloc_wrk) call prec%free_wrk(info)
  
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err=trim(method))
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_dkrylov_vect

