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
  

subroutine psb_z_cuda_hlg_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use hlldev_mod
  use psb_vectordev_mod
  use psb_z_cuda_hlg_mat_mod, psb_protect_name => psb_z_cuda_hlg_inner_vect_sv
#else 
  use psb_z_cuda_hlg_mat_mod
#endif
  use psb_z_cuda_vect_mod
  implicit none 
  class(psb_z_cuda_hlg_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)                 :: alpha, beta
  class(psb_z_base_vect_type), intent(inout) :: x, y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)            :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='z_base_inner_vect_sv'
  logical, parameter :: debug=.false.
  complex(psb_dpk_), allocatable      :: rx(:), ry(:)

  call psb_get_erraction(err_act)
  info = psb_success_

  
  call x%sync()
  call y%sync()
  if (a%is_dev()) call a%sync()
  call a%psb_z_hll_sparse_mat%inner_spsm(alpha,x,beta,y,info,trans)
  call y%set_host()

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_z_cuda_hlg_inner_vect_sv
