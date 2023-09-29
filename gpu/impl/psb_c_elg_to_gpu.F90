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
  

subroutine psb_c_elg_to_gpu(a,info,nzrm) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use elldev_mod
  use psb_vectordev_mod
  use psb_c_elg_mat_mod, psb_protect_name => psb_c_elg_to_gpu
#else 
  use psb_c_elg_mat_mod
#endif
  implicit none 
  class(psb_c_elg_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: nzrm

  integer(psb_ipk_)  :: m, nzm, n, pitch,maxrowsize, nzt 
#ifdef HAVE_SPGPU
  type(elldev_parms) :: gpu_parms
#endif

  info = 0

#ifdef HAVE_SPGPU
  if ((.not.allocated(a%val)).or.(.not.allocated(a%ja))) return
  
  m   = a%get_nrows()
  nzm = psb_size(a%val,2)
  n   = a%get_ncols()
  nzt = a%get_nzeros()
  if (present(nzrm)) nzm = max(nzm,nzrm)
  
  gpu_parms = FgetEllDeviceParams(m,nzm,nzt,n,spgpu_type_complex_float,1)
  
  if (c_associated(a%deviceMat)) then 
    pitch      = getEllDevicePitch(a%deviceMat)
    maxrowsize = getEllDeviceMaxRowSize(a%deviceMat)
  else
    pitch      = -1
    maxrowsize = -1
  end if

  if ((pitch /= gpu_parms%pitch).or.(maxrowsize /= gpu_parms%maxRowSize)) then 
    if (c_associated(a%deviceMat)) then 
      call freeEllDevice(a%deviceMat)
    endif
    info       = FallocEllDevice(a%deviceMat,m,nzm,nzt,n,spgpu_type_complex_float,1)
    pitch      = getEllDevicePitch(a%deviceMat)
    maxrowsize = getEllDeviceMaxRowSize(a%deviceMat)
  end if
  if (info == 0) then 
    if ((pitch /= psb_size(a%val,1)).or.(maxrowsize /= psb_size(a%val,2))) then 
      call psb_realloc(pitch,maxrowsize,a%val,info)
      if (info == 0) call psb_realloc(pitch,maxrowsize,a%ja,info)
    end if
  end if
  if (info == 0)  info = &
       & writeEllDevice(a%deviceMat,a%val,a%ja,size(a%ja,1),a%irn,a%idiag)
  call a%set_sync()
#endif

end subroutine psb_c_elg_to_gpu
