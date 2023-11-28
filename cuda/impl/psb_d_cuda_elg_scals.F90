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
  

subroutine psb_d_cuda_elg_scals(d,a,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use elldev_mod
  use psb_vectordev_mod
  use psb_d_cuda_elg_mat_mod, psb_protect_name => psb_d_cuda_elg_scals 
#else 
  use psb_d_cuda_elg_mat_mod
#endif
  implicit none 
  class(psb_d_cuda_elg_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)   :: info

  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_dev()) call a%sync()
  if (a%is_unit()) then 
    call a%make_nonunit()
  end if

  a%val(:,:) = a%val(:,:) * d

#ifdef HAVE_SPGPU
  call a%to_gpu(info)
  if (info /= 0) goto 9999
#endif

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cuda_elg_scals
