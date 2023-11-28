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
  

subroutine psb_s_cuda_cp_hlg_from_fmt(a,b,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use hlldev_mod
  use psb_vectordev_mod
  use psb_s_cuda_hlg_mat_mod, psb_protect_name => psb_s_cuda_cp_hlg_from_fmt
#else 
  use psb_s_cuda_hlg_mat_mod
#endif
  implicit none 

  class(psb_s_cuda_hlg_sparse_mat), intent(inout) :: a
  class(psb_s_base_sparse_mat), intent(in)   :: b
  integer(psb_ipk_), intent(out)             :: info

  info = psb_success_

  select type(b)
  type is (psb_s_coo_sparse_mat)
    call a%cp_from_coo(b,info) 
  class default
    call a%psb_s_hll_sparse_mat%cp_from_fmt(b,info)
#ifdef HAVE_SPGPU
    if (info == 0) call a%to_gpu(info)
#endif
  end select
  if (info /= 0) goto 9999

  return

9999 continue
  info = psb_err_alloc_dealloc_
  return

end subroutine psb_s_cuda_cp_hlg_from_fmt
