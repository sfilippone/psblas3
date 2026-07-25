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
! File: psb_linsolve_mod.f90
!  Interfaces for linear solvers.
!
module psb_linsolve_mod
  use psb_const_mod
  public 

  interface psb_krylov
    subroutine psb_skrylov_vect(method, a, prec, b, x, eps, desc_a, info, &
                          & itmax, iter, err, itrace, irst, istop, cond, &
                          & steps, base_type, eigext, Gram_solver, FGS_sweeps)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_sspmat_type, psb_spk_, psb_s_vect_type
      use psb_prec_mod, only : psb_sprec_type
      character(len=*)                      :: method
      type(psb_sspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      class(psb_sprec_type), intent(inout)  :: prec 
      type(psb_s_vect_type), intent(inout)  :: b
      type(psb_s_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in)            :: eps
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop, steps
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err, cond
      character, optional, intent(in)           :: base_type
      real(psb_spk_), optional, intent(in)      :: eigext(2)
      character(len=3), optional, intent(in)    :: Gram_solver
      integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps
    end subroutine psb_skrylov_vect

    subroutine psb_ckrylov_vect(method, a, prec, b, x, eps, desc_a, info, &
                          & itmax, iter, err, itrace, irst, istop, cond, &
                          & steps, base_type, eigext, Gram_solver, FGS_sweeps)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_cspmat_type, psb_spk_, psb_c_vect_type
      use psb_prec_mod, only : psb_cprec_type
      character(len=*)                      :: method
      type(psb_cspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      class(psb_cprec_type), intent(inout)  :: prec 
      type(psb_c_vect_type), intent(inout)  :: b
      type(psb_c_vect_type), intent(inout)  :: x
      real(psb_spk), intent(in)            :: eps
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop, steps
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk), optional, intent(out)     :: err, cond
      character, optional, intent(in)           :: base_type
      real(psb_spk), optional, intent(in)      :: eigext(2)
      character(len=3), optional, intent(in)    :: Gram_solver
      integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps
    end subroutine psb_ckrylov_vect

    subroutine psb_dkrylov_vect(method, a, prec, b, x, eps, desc_a, info, &
                          & itmax, iter, err, itrace, irst, istop, cond, &
                          & steps, base_type, eigext, Gram_solver, FGS_sweeps)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_dspmat_type, psb_dpk_, psb_d_vect_type
      use psb_prec_mod, only : psb_dprec_type
      character(len=*)                      :: method
      type(psb_dspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      class(psb_dprec_type), intent(inout)  :: prec 
      type(psb_d_vect_type), intent(inout)  :: b
      type(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in)            :: eps
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop, steps
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_dpk_), optional, intent(out)     :: err, cond
      character, optional, intent(in)           :: base_type
      real(psb_dpk_), optional, intent(in)      :: eigext(2)
      character(len=3), optional, intent(in)    :: Gram_solver
      integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps
    end subroutine psb_zkrylov_vect

    subroutine psb_zkrylov_vect(method, a, prec, b, x, eps, desc_a, info, &
                          & itmax, iter, err, itrace, irst, istop, cond, &
                          & steps, base_type, eigext, Gram_solver, FGS_sweeps)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_zspmat_type, psb_dpk_, psb_z_vect_type
      use psb_prec_mod, only : psb_zprec_type
      character(len=*)                      :: method
      type(psb_zspmat_type), intent(in)     :: a
      type(psb_desc_type), intent(in)       :: desc_a
      class(psb_zprec_type), intent(inout)  :: prec 
      type(psb_z_vect_type), intent(inout)  :: b
      type(psb_z_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in)            :: eps
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, irst, istop, steps
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_dpk_), optional, intent(out)     :: err, cond
      character, optional, intent(in)           :: base_type
      real(psb_dpk_), optional, intent(in)      :: eigext(2)
      character(len=3), optional, intent(in)    :: Gram_solver
      integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps
    end subroutine psb_zkrylov_vect
  end interface

  interface psb_richardson
    subroutine psb_srichardson_vect(a, prec, b, x, eps, desc_a, info, &
                                  & itmax, iter, err, itrace, istop)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_sspmat_type, psb_spk_, psb_s_vect_type
      use psb_prec_mod, only : psb_sprec_type
      type(psb_sspmat_type), intent(in)     :: a
      class(psb_sprec_type), intent(inout)  :: prec 
      type(psb_s_vect_type), intent(inout)  :: b
      type(psb_s_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in)            :: eps
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err
    end subroutine psb_srichardson_vect

    subroutine psb_crichardson_vect(a, prec, b, x, eps, desc_a, info, &
                                  & itmax, iter, err, itrace, istop)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_cspmat_type, psb_spk_, psb_c_vect_type
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)     :: a
      class(psb_cprec_type), intent(inout)  :: prec 
      type(psb_c_vect_type), intent(inout)  :: b
      type(psb_c_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in)            :: eps
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_spk_), optional, intent(out)     :: err
    end subroutine psb_crichardson_vect

    subroutine psb_drichardson_vect(a, prec, b, x, eps, desc_a, info, &
                                  & itmax, iter, err, itrace, istop)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_dspmat_type, psb_dpk_, psb_d_vect_type
      use psb_prec_mod, only : psb_dprec_type
      type(psb_dspmat_type), intent(in)     :: a
      class(psb_dprec_type), intent(inout)  :: prec 
      type(psb_d_vect_type), intent(inout)  :: b
      type(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in)            :: eps
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_dpk_), optional, intent(out)     :: err
    end subroutine psb_drichardson_vect

    subroutine psb_zrichardson_vect(a, prec, b, x, eps, desc_a, info, &
                                  & itmax, iter, err, itrace, istop)
      use psb_base_mod, only : psb_ipk_, psb_desc_type, psb_zspmat_type, psb_dpk_, psb_z_vect_type
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)     :: a
      class(psb_zprec_type), intent(inout)  :: prec 
      type(psb_z_vect_type), intent(inout)  :: b
      type(psb_z_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in)            :: eps
      type(psb_desc_type), intent(in)       :: desc_a
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), optional, intent(out)  :: iter
      real(psb_dpk_), optional, intent(out)     :: err
    end subroutine psb_zrichardson_vect
  end interface
end module psb_linsolve_mod
