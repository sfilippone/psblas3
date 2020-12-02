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
!    Moved here from AMG-AINV, original copyright below.
!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
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
!
!
!
module psb_s_ainv_tools_mod

  interface sp_drop
    subroutine psb_s_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
      use psb_base_mod, only : psb_spk_, psb_ipk_
      implicit none
      real(psb_spk_), intent(in)      :: sp_thresh
      integer(psb_ipk_), intent(in)     :: idiag, nzrmax
      integer(psb_ipk_), intent(inout)  :: nz
      integer(psb_ipk_), intent(inout)  :: iz(:)
      real(psb_spk_), intent(inout)    :: valz(:)
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psb_s_sp_drop
  end interface

  interface rwclip
    subroutine psb_s_rwclip(nz,ia,ja,val,imin,imax,jmin,jmax)
      use psb_base_mod, only : psb_spk_, psb_ipk_

      implicit none
      integer(psb_ipk_), intent(inout)  :: nz
      integer(psb_ipk_), intent(inout)  :: ia(*), ja(*)
      real(psb_spk_), intent(inout)    :: val(*)
      integer(psb_ipk_), intent(in)     :: imin,imax,jmin,jmax
    end subroutine psb_s_rwclip
  end interface

  interface sparsify
    subroutine psb_s_sparsify(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,info, &
         & istart,iheap,ikr)
      use psb_base_mod, only : psb_spk_, psb_ipk_, psb_i_heap
      implicit none

      real(psb_spk_), intent(in)            :: sp_thresh
      integer(psb_ipk_), intent(in)           :: idiag, n, nzrmax
      real(psb_spk_), intent(inout)          :: zw(:)
      integer(psb_ipk_), intent(out)          :: nz
      integer(psb_ipk_), intent(out)          :: iz(:)
      real(psb_spk_), intent(out)            :: valz(:)
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), intent(in), optional :: istart
      type(psb_i_heap), optional              :: iheap
      integer(psb_ipk_), optional             :: ikr(:)
    end subroutine psb_s_sparsify
    subroutine psb_s_sparsify_list(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,lhead,listv,ikr,info)
      use psb_base_mod, only : psb_spk_, psb_ipk_
      implicit none

      real(psb_spk_), intent(in)     :: sp_thresh
      integer(psb_ipk_), intent(in)    :: idiag, n, nzrmax
      real(psb_spk_), intent(inout)   :: zw(:)
      integer(psb_ipk_), intent(out)   :: nz
      integer(psb_ipk_), intent(out)   :: iz(:)
      real(psb_spk_), intent(out)      :: valz(:)
      integer(psb_ipk_), intent(out)   :: info
      integer(psb_ipk_), intent(inout) :: lhead, listv(:)
      integer(psb_ipk_)                :: ikr(:)
    end subroutine psb_s_sparsify_list

  end interface

end module psb_s_ainv_tools_mod
