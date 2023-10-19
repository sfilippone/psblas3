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
  
subroutine psb_z_csrg_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use elldev_mod
  use psb_vectordev_mod
  use psb_z_csrg_mat_mod, psb_protect_name => psb_z_csrg_inner_vect_sv
#else 
  use psb_z_csrg_mat_mod
#endif
  use psb_z_gpu_vect_mod
  implicit none 
  class(psb_z_csrg_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)                :: alpha, beta
  class(psb_z_base_vect_type), intent(inout) :: x, y
  integer(psb_ipk_), intent(out)           :: info
  character, optional, intent(in)          :: trans

  complex(psb_dpk_), allocatable      :: rx(:), ry(:)
  logical        :: tra
  character      :: trans_
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='z_csrg_inner_vect_sv'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_success_


  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')

#ifdef HAVE_SPGPU  
  if (tra.or.(beta/=dzero)) then 
    call x%sync()
    call y%sync()
    call a%psb_z_csr_sparse_mat%inner_spsm(alpha,x,beta,y,info,trans)
    call y%set_host()
  else
    select type (xx => x) 
    type is (psb_z_vect_gpu)
      select type(yy => y) 
      type is (psb_z_vect_gpu)
        if (xx%is_host()) call xx%sync()
        if (beta /= dzero) then 
          if (yy%is_host()) call yy%sync()
        end if
        info = spsvCSRGDevice(a%deviceMat,alpha,xx%deviceVect,&
             & beta,yy%deviceVect)
        if (info /= 0) then 
          call psb_errpush(psb_err_from_subroutine_ai_,name,&
               & a_err='spsvCSRGDevice',i_err=(/info,izero,izero,izero,izero/))
          info = psb_err_from_subroutine_ai_
          goto 9999
        end if
        call yy%set_dev()
      class default
        rx = xx%get_vect()
        ry = y%get_vect()
        call a%psb_z_csr_sparse_mat%inner_spsm(alpha,rx,beta,ry,info)
        call y%bld(ry)
      end select
    class default
      rx = x%get_vect()
      ry = y%get_vect()
      call a%psb_z_csr_sparse_mat%inner_spsm(alpha,rx,beta,ry,info)
      call y%bld(ry)
    end select
  end if
#else
  call x%sync()
  call y%sync()
  call a%psb_z_csr_sparse_mat%inner_spsm(alpha,x,beta,y,info,trans)
  call y%set_host()
#endif
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='csrg_vect_sv')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_z_csrg_inner_vect_sv
