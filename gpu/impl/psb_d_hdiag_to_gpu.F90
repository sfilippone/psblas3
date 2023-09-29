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
  

subroutine psb_d_hdiag_to_gpu(a,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use hdiagdev_mod
  use psb_vectordev_mod
  use psb_d_hdiag_mat_mod, psb_protect_name => psb_d_hdiag_to_gpu
#else 
  use psb_d_hdiag_mat_mod
#endif
  use iso_c_binding
  implicit none 
  class(psb_d_hdiag_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_) :: nr, nc, hacksize, hackcount, allocheight 
#ifdef HAVE_SPGPU
  type(hdiagdev_parms) :: gpu_parms
#endif

  info = 0

#ifdef HAVE_SPGPU
  nr = a%get_nrows()
  nc = a%get_ncols()
  hacksize  = a%hackSize
  hackCount = a%nhacks
  if (.not.allocated(a%hackOffsets)) then 
    info = -1 
    return
  end if
  allocheight = a%hackOffsets(hackCount+1) 
!!$  write(*,*) 'HDIAG TO GPU:',nr,nc,hacksize,hackCount,allocheight,&
!!$       & size(a%hackoffsets),size(a%diaoffsets), size(a%val)
  if (.not.allocated(a%diaOffsets)) then
    info = -2 
    return
  end if
  if (.not.allocated(a%val)) then
    info = -3
    return
  end if

  if (c_associated(a%deviceMat)) then 
     call freeHdiagDevice(a%deviceMat)
  endif

  info = FAllocHdiagDevice(a%deviceMat,nr,nc,&
       & allocheight,hacksize,hackCount,spgpu_type_double)
  if (info == 0) info = &
       & writeHdiagDevice(a%deviceMat,a%val,a%diaOffsets,a%hackOffsets)

#endif

end subroutine psb_d_hdiag_to_gpu
