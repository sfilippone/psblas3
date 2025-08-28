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
!    Moved here from MLD2P4, original copyright below.
!  
!   
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!    
!    (C) Copyright 2008-2018 
!  
!        Salvatore Filippone  
!        Pasqua D'Ambra   
!        Daniela di Serafino   
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
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
!
! File: psb_d_ilu_fact_mod.f90
!
! Module: psb_d_ilu_fact_mod
!
!  This module defines some interfaces used internally by the implementation of
!  psb_d_ilu_solver, but not visible to the end user. 
!
!
module psb_d_ilu_fact_mod
  use psb_base_mod
  use psb_prec_const_mod
  interface psb_ilu0_fact
    subroutine psb_dilu0_fact(ialg,a,l,u,d,info,blck,upd,shft)
      import psb_dspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)         :: ialg
      integer(psb_ipk_), intent(out)        :: info
      type(psb_dspmat_type),intent(in)      :: a
      type(psb_dspmat_type),intent(inout)   :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      character, intent(in), optional       :: upd
      real(psb_dpk_), intent(inout)      :: d(:)
      real(psb_dpk_), intent(in), optional :: shft
    end subroutine psb_dilu0_fact
  end interface

  interface psb_iluk_fact
    subroutine psb_diluk_fact(fill_in,ialg,a,l,u,d,info,blck,shft)
      import psb_dspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in,ialg
      integer(psb_ipk_), intent(out)       :: info
      type(psb_dspmat_type),intent(in)     :: a
      type(psb_dspmat_type),intent(inout)  :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
      real(psb_dpk_), intent(in), optional :: shft
    end subroutine psb_diluk_fact
  end interface

  interface psb_ilut_fact
    subroutine psb_dilut_fact(fill_in,thres,a,l,u,d,info,blck,iscale,shft)
      import  psb_dspmat_type, psb_dpk_, psb_ipk_
      integer(psb_ipk_), intent(in)        :: fill_in
      real(psb_dpk_), intent(in)           :: thres
      integer(psb_ipk_), intent(out)       :: info
      type(psb_dspmat_type),intent(in)     :: a
      type(psb_dspmat_type),intent(inout)  :: l,u
      real(psb_dpk_), intent(inout)     :: d(:)
      type(psb_dspmat_type),intent(in), optional, target :: blck
      integer(psb_ipk_), intent(in), optional  :: iscale
      real(psb_dpk_), intent(in), optional :: shft
    end subroutine psb_dilut_fact
  end interface

end module psb_d_ilu_fact_mod
