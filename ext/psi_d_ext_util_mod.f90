!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

module psi_d_ext_util_mod

  use psb_base_mod, only : psb_ipk_, psb_dpk_

  interface psi_xtr_dia_from_coo
    subroutine psi_d_xtr_dia_from_coo(nr,nc,nz,ia,ja,val,d,nrd,ncd,data,info,&
         & initdata,rdisp)    
      import  :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_), intent(in)  :: nr, nc, nz, nrd, ncd, ia(:), ja(:), d(:)
      real(psb_dpk_),    intent(in)  :: val(:)
      real(psb_dpk_),    intent(out) :: data(nrd,ncd)
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional  :: initdata
      integer(psb_ipk_), intent(in), optional :: rdisp
      
    end subroutine psi_d_xtr_dia_from_coo
  end interface

  interface psi_xtr_ell_from_coo
    subroutine psi_d_xtr_ell_from_coo(i,nr,mxrwl,iac,jac,&
         & valc,ja,val,irn,diag,ld)
      import  :: psb_ipk_, psb_dpk_
      implicit none 
      integer(psb_ipk_) :: i,nr,mxrwl,ld
      integer(psb_ipk_) :: iac(*),jac(*),ja(ld,*),irn(*),diag(*)
      real(psb_dpk_)    :: valc(*), val(ld,*)
      
    end subroutine psi_d_xtr_ell_from_coo
  end interface psi_xtr_ell_from_coo
      
  interface psi_xtr_coo_from_dia
    subroutine psi_d_xtr_coo_from_dia(nr,nc,ia,ja,val,nz,nrd,ncd,data,offsets,&
         & info,rdisp)
      import :: psb_ipk_, psb_dpk_
      
      implicit none 
      
      integer(psb_ipk_), intent(in)    :: nr,nc, nrd,ncd, offsets(:) 
      integer(psb_ipk_), intent(inout) :: ia(:), ja(:), nz
      real(psb_dpk_),    intent(inout) :: val(:)
      real(psb_dpk_),    intent(in)    :: data(nrd,ncd)
      integer(psb_ipk_), intent(out)   :: info
      integer(psb_ipk_), intent(in), optional :: rdisp
    end subroutine psi_d_xtr_coo_from_dia
  end interface
  
end module psi_d_ext_util_mod
