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
  

subroutine psb_z_cuda_hlg_to_gpu(a,info,nzrm) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use hlldev_mod
  use psb_vectordev_mod
  use psb_z_cuda_hlg_mat_mod, psb_protect_name => psb_z_cuda_hlg_to_gpu
#else 
  use psb_z_cuda_hlg_mat_mod
#endif
  use iso_c_binding
  implicit none 
  class(psb_z_cuda_hlg_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: nzrm

  integer(psb_ipk_) :: m, nzm, nza, n, pitch,maxrowsize, allocsize

  info = 0

#ifdef HAVE_SPGPU
  if ((.not.allocated(a%val)).or.(.not.allocated(a%ja))) return
  
  n   = a%get_nrows()
  allocsize = a%get_size()
  nza = a%get_nzeros()
  if (c_associated(a%deviceMat)) then 
     call freehllDevice(a%deviceMat)
  endif
  info       = FallochllDevice(a%deviceMat,a%hksz,n,nza,allocsize,spgpu_type_complex_double,1)
  if (info == 0)  info = &
       & writehllDevice(a%deviceMat,a%val,a%ja,a%hkoffs,a%irn,a%idiag)
!  if (info /= 0) goto 9999
#endif

end subroutine psb_z_cuda_hlg_to_gpu
