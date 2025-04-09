!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
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
!
! File: psb_dvmlt.f90
!!$
!!$subroutine psb_dvmlt(x,y,desc_a,info)
!!$  use psb_base_mod, psb_protect_name => psb_dvmlt
!!$  implicit none
!!$  type(psb_d_vect_type), intent (inout) ::  x
!!$  type(psb_d_vect_type), intent (inout) ::  y
!!$  type(psb_desc_type), intent (in)      :: desc_a
!!$  integer(psb_ipk_), intent(out)                  :: info
!!$
!!$  ! locals
!!$  integer(psb_ipk_) :: ctxt, np, me,&
!!$       & err_act, iix, jjx, iiy, jjy
!!$  integer(psb_lpk_) :: ix, ijx, iy, ijy, m
!!$  character(len=20)        :: name, ch_err
!!$
!!$  name='psb_dgevmlt'
!!$  if (psb_errstatus_fatal()) return
!!$  info=psb_success_
!!$  call psb_erractionsave(err_act)
!!$
!!$  ctxt=desc_a%get_context()
!!$
!!$  call psb_info(ctxt, me, np)
!!$  if (np == -ione) then
!!$    info = psb_err_context_error_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$  if (.not.allocated(x%v)) then
!!$    info = psb_err_invalid_vect_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$  if (.not.allocated(y%v)) then
!!$    info = psb_err_invalid_vect_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$
!!$  ix = ione
!!$  iy = ione
!!$
!!$  m = desc_a%get_global_rows()
!!$
!!$  ! check vector correctness
!!$  call psb_chkvect(m,lone,x%get_nrows(),ix,lone,desc_a,info,iix,jjx)
!!$  if(info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    ch_err='psb_chkvect 1'
!!$    call psb_errpush(info,name,a_err=ch_err)
!!$    goto 9999
!!$  end if
!!$  call psb_chkvect(m,lone,y%get_nrows(),iy,lone,desc_a,info,iiy,jjy)
!!$  if(info /= psb_success_) then
!!$    info=psb_err_from_subroutine_
!!$    ch_err='psb_chkvect 2'
!!$    call psb_errpush(info,name,a_err=ch_err)
!!$    goto 9999
!!$  end if
!!$
!!$  if ((iix /= ione).or.(iiy /= ione)) then
!!$    info=psb_err_ix_n1_iy_n1_unsupported_
!!$    call psb_errpush(info,name)
!!$  end if
!!$
!!$  if(desc_a%get_local_rows() > 0) then
!!$    call y%base_mlt_v(desc_a%get_local_rows(),&
!!$         & alpha,x,beta,info)
!!$  end if
!!$
!!$  call psb_erractionrestore(err_act)
!!$  return
!!$
!!$9999 call psb_error_handler(ctxt,err_act)
!!$
!!$  return
!!$
!!$end subroutine psb_dvmlt

!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
!    
! File: psb_dvmlt.f90
!
! Function: psb_ddot_multivect
!    psb_dvmlta computes the inner product of two distributed multivectors,
!
!    dot(:,:) := ( X(:,:) )**C * ( Y(:,:) )
!
!
! Arguments:
!    x      -  type(psb_d_multivect_type) The input vector containing the entries of sub( X ).
!    y      -  type(psb_d_multivect_type) The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).  The communication descriptor.
!    info   -  integer.              Return code
!    global -  logical(optional)     Whether to perform the global sum, default: .true.
!
!  Note: from a functional point of view, X and Y are input, but here
!        they are declared INOUT because of the sync() methods. 
!
!
subroutine psb_smlt_multivect(x, y, res,desc_a,info,global) 
  use psb_desc_mod
  use psb_s_base_mat_mod
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_s_multivect_mod
  use psb_s_psblas_mod, psb_protect_name => psb_dmlt_multivect
  implicit none 
  real(psb_spk_), dimension(:,:), allocatable, intent(inout) :: res
  type(psb_s_multivect_type), intent(inout) :: x, y
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  logical, intent(in), optional        :: global

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, idx, ndm,&
       & err_act, iix, jjx, iiy, jjy, i, nr
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m, n, nx, ny
  logical :: global_
  character(len=20)      :: name, ch_err

  name='psb_dmlt_multivect'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (np == -ione) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(y%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(global)) then
    global_ = global
  else
    global_ = .true.
  end if

  ix = ione
  ijx = ione

  iy = ione
  ijy = ione

  m = desc_a%get_global_rows()
  nx = x%get_ncols()
  ny = y%get_ncols()

  ! check vector correctness
  call psb_chkvect(m,nx,x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ny,y%get_nrows(),iy,ijy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((iix /= ione).or.(iiy /= ione)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (x%get_ncols() /= y%get_ncols()) then
    info=psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  else
    if (allocated(res)) then
      if ((size(res,1) /= x%get_ncols()).or.(size(res,2) /= y%get_ncols())) then
        deallocate(res,stat=info)
      end if
    end if
  end if

  nr = desc_a%get_local_rows() 
  if(nr > 0) then
    call x%mlt(nr,y,res,info) 
    ! FIXME
    ! adjust dot_local because overlapped elements are computed more than once
    if (size(desc_a%ovrlap_elem,1)>0) then
      if (x%is_dev()) call x%sync()
      if (y%is_dev()) call y%sync()
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        ! FIXME: case of AS-type descriptors
        ! Since I'm coputing res via a dgemm on the whole vector, I need to adjust
        ! the result by removing the contribution of the overlapped elements
        ! specifically: res(:,:) = res(:,:) - x%v%v(idx,:)^T*y%v%v(idx,:)
        ! using dgemm to compute the matrix-matrix product of the form R = R - X'*Y
        ! where R is the result, X' is the transpose of the matrix x%v%v(idx,:) 
        ! and Y is the matrix y%v%v(idx,:)
        ! call dgemm('T','N',size(x%v%v(idx,:),1),size(y%v%v(idx,:),1),&
        !      & size(x%v%v(idx,:),1),-done,x%v%v(idx,:),size(x%v%v(idx,:),1),&
        !      & y%v%v(idx,:),size(y%v%v(idx,:),1),done,res,y%get_ncols())
        info = psb_err_internal_error_
        ch_err='over_elem_unsup'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end do
    end if
  else
    res = dzero
  end if

  ! compute global sum
  if (global_) call psb_sum(ctxt, res)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_smlt_multivect
