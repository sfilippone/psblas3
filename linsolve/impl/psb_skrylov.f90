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
!         software without specific prior written permission.
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
! Subroutine: psb_skrylov
! 
!    Front-end for the Krylov subspace iterations, realversion
!    
! Arguments:
!
!    methd  -  character                    The specific method; can take the values:
!                                           CG
!                                           FCG
!                                           CGS
!                                           BICG
!                                           BICGSTAB
!                                           BICGSTABL
!                                           RGMRES
!                                           
!    a      -  type(psb_sspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_sprec_type)       Input: preconditioner
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
subroutine psb_skrylov_vect(method, a, prec, b, x, eps, desc_a, info, &
                        & itmax, iter, err, itrace, irst, istop, cond, &
                        & steps, base_type, eigext, Gram_solver, FGS_sweeps)

  use psb_base_mod
  use psb_prec_mod, only : psb_sprec_type
  use psb_linsolve_mod, psb_protect_name => psb_skrylov_vect

  character(len=*)                     :: method
  type(psb_sspmat_type), intent(in)    :: a
  type(psb_desc_type), intent(in)      :: desc_a
  class(psb_sprec_type), intent(inout) :: prec 
  type(psb_s_vect_type), intent(inout) :: b
  type(psb_s_vect_type), intent(inout) :: x
  real(psb_spk_), intent(in)           :: eps
  integer(psb_ipk_), intent(out)       :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop, steps
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_spk_), optional, intent(out)     :: err, cond
  character, optional, intent(in)           :: base_type
  real(psb_spk_), optional, intent(in)      :: eigext(2)
  character(len=3), optional, intent(in)    :: Gram_solver
  integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps

  abstract interface
    subroutine psb_skryl_vect(a, prec, b, x, eps, desc_a, info, &
                          & itmax, iter, err, itrace, istop)
      import :: psb_ipk_, psb_spk_, psb_desc_type, &
           & psb_sspmat_type, psb_sprec_type, psb_s_vect_type
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      type(psb_s_vect_type), intent(inout) :: b
      type(psb_s_vect_type), intent(inout) :: x
      real(psb_spk_), intent(in)           :: eps
      class(psb_sprec_type), intent(inout) :: prec
      integer(psb_ipk_), intent(out)           :: info
      integer(psb_ipk_), optional, intent(in)  :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out) :: iter
      real(psb_spk_), optional, intent(out)    :: err
    end subroutine psb_skryl_vect

    subroutine psb_skryl_rest_vect(a, prec, b, x, eps, desc_a, info, &
                                & itmax, iter, err, itrace, irst, istop)
      import :: psb_ipk_, psb_spk_, psb_desc_type, &
           & psb_sspmat_type, psb_sprec_type, psb_s_vect_type
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      class(psb_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type), intent(inout) :: b
      type(psb_s_vect_type), intent(inout) :: x
      real(psb_spk_), intent(in)           :: eps
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err
    end subroutine psb_skryl_rest_vect

    subroutine psb_skryl_cond_vect(a, prec, b, x, eps, desc_a, info, &
         & itmax, iter, err, itrace, istop, cond)
      import :: psb_ipk_, psb_spk_, psb_desc_type, &
           & psb_sspmat_type, psb_sprec_type, psb_s_vect_type
      type(psb_sspmat_type), intent(in)    :: a
      type(psb_desc_type), intent(in)      :: desc_a
      class(psb_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type), intent(inout) :: b
      type(psb_s_vect_type), intent(inout) :: x
      real(psb_spk_), intent(in)           :: eps
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err, cond
    end subroutine psb_skryl_cond_vect

    subroutine psb_skryl_step_vect(a, prec, b, x, s, eps, desc_a, info, &
                                & itmax, iter, err, itrace, istop, &
                                & base_type, eigext, Gram_solver, FGS_sweeps)
      import :: psb_ipk_, psb_spk_, psb_desc_type, &
           & psb_sspmat_type, psb_sprec_type, psb_s_vect_type
      type(psb_sspmat_type), intent(in)    :: a
      class(psb_sprec_type), intent(inout) :: prec
      type(psb_s_vect_type), intent(inout) :: b
      type(psb_s_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)        :: s
      real(psb_spk_), intent(in)           :: eps
      type(psb_desc_type), intent(in)      :: desc_a
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err
      character, optional, intent(in)           :: base_type
      real(psb_spk_), optional, intent(in)      :: eigext(2)
      character(len=3), optional, intent(in)    :: Gram_solver
      integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps
    end subroutine psb_skryl_step_vect
  end interface

  procedure(psb_skryl_vect)       :: psb_sbicg_vect, psb_scgstab_vect, psb_scgs_vect
  procedure(psb_skryl_rest_vect)  :: psb_srgmres_vect, psb_scgstabl_vect, psb_sgcr_vect
  procedure(psb_skryl_cond_vect)  :: psb_scg_vect, psb_sfcg_vect, psb_sminres_vect
  procedure(psb_skryl_step_vect)  :: psb_sscg_vect, psb_sscg2_vect

  logical             :: do_alloc_wrk
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np, err_act, itrace_, steps_
  character(len=20)   :: name

  info = psb_success_
  name = 'psb_krylov'
  call psb_erractionsave(err_act)

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)

  ! Default return for COND
  if (present(cond)) cond = szero

  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = -1
  end if

  do_alloc_wrk = .not. prec%is_allocated_wrk()
  if (do_alloc_wrk) call prec%allocate_wrk(info, vmold=x%v, desc=desc_a)

  select case(psb_toupper(method))
    case('CG') 
      call psb_scg_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop, cond = cond)
    case('FCG') 
      call psb_sfcg_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop, cond = cond)
    case('GCR') 
      call psb_sgcr_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop)
    case('CGS') 
      call psb_scgs_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop)
    case('BICG') 
      call psb_sbicg_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop)
    case('BICGSTAB') 
      call psb_scgstab_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop)
    case('RGMRES', 'GMRES')
      call psb_srgmres_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, irst = irst, istop = istop)
    case('MINRES', 'PMINRES')
      call psb_sminres_vect(a, prec, b, x, eps, desc_a, info, &
		  & itmax, iter, err, itrace = itrace_, istop = istop)
    case('BICGSTABL')
      call psb_scgstabl_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, irst = irst, istop = istop)
    case('SSTEPCG')
      ! steps (default = 5)
      steps_ = 5
      if(present(steps)) steps_ = steps
      
      call psb_sscg_vect(a, prec, b, x, steps_, eps, desc_a, info, &
          & itmax = itmax, iter = iter, err = err, itrace = itrace_, istop = istop, &
          & base_type = base_type, eigext = eigext, Gram_solver = Gram_solver, FGS_sweeps = FGS_sweeps)
      
    case('SSTEPCG1')
      ! steps (default = 5)
      steps_ = 5
      if(present(steps)) steps_ = steps

      call psb_sscg2_vect(a, prec, b, x, steps_, eps, desc_a, info, &
          & itmax = itmax, iter = iter, err = err, itrace = itrace_, istop = istop, &
          & base_type = base_type, eigext = eigext, Gram_solver = Gram_solver, FGS_sweeps = FGS_sweeps)  
                     
    case default
      if (me == psb_root_) write(psb_err_unit, *) trim(name) , &
          & ': Warning: Unknown method  ', method, ', defaulting to BiCGSTAB'
      
      call psb_scgstab_vect(a, prec, b, x, eps, desc_a, info, &
          & itmax, iter, err, itrace = itrace_, istop = istop)
  end select

  if ((info == psb_success_) .and. do_alloc_wrk) call prec%free_wrk(info)
  
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name, a_err = trim(method))
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt, err_act)
  return
end subroutine psb_skrylov_vect