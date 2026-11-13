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
Module psb_c_linsolve_mod

  use psb_const_mod
  public 

  interface psb_krylov
    Subroutine psb_ckrylov_vect(method,a,prec,b,x,eps,desc_a,info,&
         & itmax,iter,err,itrace,irst,istop,cond,s1,s2)
      use psb_base_mod, only  : psb_ipk_, psb_desc_type, psb_cspmat_type, &
           & psb_spk_, psb_c_vect_type
      use psb_prec_mod, only : psb_cprec_type
      character(len=*)                      :: method
      Type(psb_cspmat_type), Intent(in)     :: a
      Type(psb_desc_type), Intent(in)       :: desc_a
      class(psb_cprec_type), intent(inout)  :: prec 
      type(psb_c_vect_type), Intent(inout)  :: b
      type(psb_c_vect_type), Intent(inout)  :: x
      Real(psb_spk_), Intent(in)             :: eps
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), Optional, Intent(in)   :: itmax, itrace, irst,istop
      integer(psb_ipk_), Optional, Intent(out)  :: iter
      Real(psb_spk_), Optional, Intent(out) :: err,cond
      type(psb_c_vect_type), optional  :: s1,s2
    end Subroutine psb_ckrylov_vect
  end interface

  interface psb_richardson
    Subroutine psb_crichardson_vect(a,prec,b,x,eps,desc_a,info,&
         & itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_ipk_, psb_desc_type, psb_cspmat_type, &
           & psb_spk_, psb_c_vect_type
      use psb_prec_mod, only : psb_cprec_type
      Type(psb_cspmat_type), Intent(in)     :: a
      Type(psb_desc_type), Intent(in)       :: desc_a
      class(psb_cprec_type), intent(inout)  :: prec 
      type(psb_c_vect_type), Intent(inout)  :: b
      type(psb_c_vect_type), Intent(inout)  :: x
      Real(psb_spk_), Intent(in)            :: eps
      integer(psb_ipk_), intent(out)            :: info
      integer(psb_ipk_), Optional, Intent(in)   :: itmax, itrace, istop
      integer(psb_ipk_), Optional, Intent(out)  :: iter
      Real(psb_spk_), Optional, Intent(out) :: err
    end Subroutine psb_crichardson_vect
  end interface
  
end module psb_c_linsolve_mod
