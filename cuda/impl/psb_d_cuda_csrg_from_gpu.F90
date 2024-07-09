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
  

subroutine psb_d_cuda_csrg_from_gpu(a,info) 
  
  use psb_base_mod
  use elldev_mod
  use psb_vectordev_mod
  use psb_d_cuda_csrg_mat_mod, psb_protect_name => psb_d_cuda_csrg_from_gpu
  implicit none 
  class(psb_d_cuda_csrg_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)             :: info

  integer(psb_ipk_)  :: m, n, nz

  info = 0

  if (.not.(c_associated(a%deviceMat%mat))) then 
    call a%free()
    return
  end if

  info = CSRGDeviceGetParms(a%deviceMat,m,n,nz)
  if (info /= psb_success_) return
  
  if (info == 0) call psb_realloc(m+1,a%irp,info)
  if (info == 0) call psb_realloc(nz,a%ja,info)
  if (info == 0) call psb_realloc(nz,a%val,info)
  if (info == 0)  info = &
       & CSRGDevice2Host(a%deviceMat,m,n,nz,a%irp,a%ja,a%val)
#if (CUDA_SHORT_VERSION <= 10) || (CUDA_VERSION <  11030)
  a%irp(:) = a%irp(:)+1
  a%ja(:) = a%ja(:)+1
#endif

  call a%set_sync()

end subroutine psb_d_cuda_csrg_from_gpu
