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
  

module psb_cuda_mod
  use psb_const_mod
  use psb_cuda_env_mod

  use psb_i_cuda_vect_mod
  use psb_s_cuda_vect_mod
  use psb_d_cuda_vect_mod
  use psb_c_cuda_vect_mod
  use psb_z_cuda_vect_mod

  use psb_i_cuda_multivect_mod
  use psb_s_cuda_multivect_mod
  use psb_d_cuda_multivect_mod
  use psb_c_cuda_multivect_mod
  use psb_z_cuda_multivect_mod

  use psb_d_ell_mat_mod
  use psb_d_cuda_elg_mat_mod
  use psb_s_ell_mat_mod
  use psb_s_cuda_elg_mat_mod
  use psb_z_ell_mat_mod
  use psb_z_cuda_elg_mat_mod
  use psb_c_ell_mat_mod
  use psb_c_cuda_elg_mat_mod

  use psb_s_hll_mat_mod
  use psb_s_cuda_hlg_mat_mod
  use psb_d_hll_mat_mod
  use psb_d_cuda_hlg_mat_mod
  use psb_c_hll_mat_mod
  use psb_c_cuda_hlg_mat_mod
  use psb_z_hll_mat_mod
  use psb_z_cuda_hlg_mat_mod
  
  use psb_s_cuda_csrg_mat_mod
  use psb_d_cuda_csrg_mat_mod
  use psb_c_cuda_csrg_mat_mod
  use psb_z_cuda_csrg_mat_mod
#if CUDA_SHORT_VERSION <= 10 
  use psb_s_cuda_hybg_mat_mod
  use psb_d_cuda_hybg_mat_mod
  use psb_c_cuda_hybg_mat_mod
  use psb_z_cuda_hybg_mat_mod
#endif
  use psb_d_cuda_diag_mat_mod
  use psb_d_cuda_hdiag_mat_mod

  use psb_s_cuda_dnsg_mat_mod
  use psb_d_cuda_dnsg_mat_mod
  use psb_c_cuda_dnsg_mat_mod
  use psb_z_cuda_dnsg_mat_mod
  
  use psb_s_cuda_hdiag_mat_mod
  ! use psb_s_cuda_diag_mat_mod

end module psb_cuda_mod

