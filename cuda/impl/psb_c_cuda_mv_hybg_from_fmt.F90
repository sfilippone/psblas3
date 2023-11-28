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
#if CUDA_SHORT_VERSION <= 10 
  
subroutine psb_c_cuda_mv_hybg_from_fmt(a,b,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use cusparse_mod
  use psb_c_cuda_hybg_mat_mod, psb_protect_name => psb_c_cuda_mv_hybg_from_fmt
#else 
  use psb_c_cuda_hybg_mat_mod
#endif
  implicit none 

  class(psb_c_cuda_hybg_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)              :: info

  !locals
  info = psb_success_

  select type(b)
  type is (psb_c_coo_sparse_mat)
    call a%mv_from_coo(b,info) 
  class default
    call a%psb_c_csr_sparse_mat%mv_from_fmt(b,info) 
    if (info /= 0) return
#ifdef HAVE_SPGPU
    call a%to_gpu(info)
#endif
  end select
end subroutine psb_c_cuda_mv_hybg_from_fmt
#endif
