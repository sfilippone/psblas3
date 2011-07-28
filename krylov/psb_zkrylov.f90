!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!
! File: psb_krylov_mod.f90
!  Interfaces for Krylov subspace iterative methods.
  !
  ! Subroutine: psb_zkrylov
  ! 
  !    Front-end for the Krylov subspace iterations, complexversion
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
  !    a      -  type(psb_zspmat_type)      Input: sparse matrix containing A.
  !    prec   -  class(psb_zprec_type)       Input: preconditioner
  !    b      -  complex,dimension(:)         Input: vector containing the
  !                                           right hand side B
  !    x      -  complex,dimension(:)         Input/Output: vector containing the
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
Subroutine psb_zkrylov(method,a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod,only : psb_sprec_type, psb_dprec_type, psb_cprec_type, psb_zprec_type
  use psb_krylov_mod, psb_protect_name => psb_zkrylov 

  character(len=*)                   :: method
  Type(psb_zspmat_type), Intent(in)  :: a
  Type(psb_desc_type), Intent(in)    :: desc_a
  class(psb_zprec_type), intent(in)   :: prec 
  complex(psb_dpk_), Intent(in)    :: b(:)
  complex(psb_dpk_), Intent(inout) :: x(:)
  Real(psb_dpk_), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
  Integer, Optional, Intent(out)     :: iter
  Real(psb_dpk_), Optional, Intent(out) :: err

  interface 
    subroutine psb_zcg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)    :: b(:)
      complex(psb_dpk_), intent(inout) :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      class(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcg
    subroutine psb_zbicg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)      :: b(:)
      complex(psb_dpk_), intent(inout)   :: x(:)
      real(psb_dpk_), intent(in)         :: eps
      class(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zbicg
    subroutine psb_zcgstab(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)       :: b(:)
      complex(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      class(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcgstab
    Subroutine psb_zcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      Type(psb_zspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      class(psb_zprec_type), intent(in)   :: prec 
      complex(psb_dpk_), Intent(in)    :: b(:)
      complex(psb_dpk_), Intent(inout) :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_zcgstabl
    Subroutine psb_zrgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      Type(psb_zspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      class(psb_zprec_type), intent(in)   :: prec 
      complex(psb_dpk_), Intent(in)    :: b(:)
      complex(psb_dpk_), Intent(inout) :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_zrgmres
    subroutine psb_zcgs(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)       :: b(:)
      complex(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      class(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcgs
  end interface

  integer                            :: ictxt,me,np,err_act
  character(len=20)             :: name

  info = psb_success_
  name = 'psb_krylov'
  call psb_erractionsave(err_act)


  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)


  select case(psb_toupper(method))
  case('CG') 
    call  psb_zcg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
  case('CGS') 
    call  psb_zcgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
  case('BICG') 
    call  psb_zbicg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
  case('BICGSTAB') 
    call  psb_zcgstab(a,prec,b,x,eps,desc_a,info,&
         & itmax,iter,err,itrace,istop)
  case('RGMRES')
    call  psb_zrgmres(a,prec,b,x,eps,desc_a,info,&
         & itmax,iter,err,itrace,irst,istop)
  case('BICGSTABL')
    call  psb_zcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
  case default
    if (me == 0) write(psb_err_unit,*) trim(name),': Warning: Unknown method  ',method,&
         & ', defaulting to BiCGSTAB'
    call  psb_zcgstab(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
  end select

  if(info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if

end subroutine psb_zkrylov

