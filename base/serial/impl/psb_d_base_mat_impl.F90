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
! == ==================================
!
!
!
! Data management
!
!
!
!
!
! == ==================================

subroutine psb_d_base_cp_to_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_cp_to_coo

subroutine psb_d_base_cp_from_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_cp_from_coo


subroutine psb_d_base_cp_to_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%cp_to_coo(b,info)
  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_to_fmt

subroutine psb_d_base_cp_from_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%cp_from_coo(b,info)
  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_from_fmt


subroutine psb_d_base_mv_to_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_coo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to coo')
    goto 9999
  end if

  call a%free()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_to_coo

subroutine psb_d_base_mv_from_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_from_coo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='from coo')
    goto 9999
  end if

  call b%free()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_from_coo


subroutine psb_d_base_mv_to_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%mv_to_coo(b,info)
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

  return

end subroutine psb_d_base_mv_to_fmt

subroutine psb_d_base_mv_from_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%mv_from_coo(b,info)
  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
  return

end subroutine psb_d_base_mv_from_fmt

subroutine  psb_d_base_clean_zeros(a, info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_clean_zeros
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)              :: info
  !
  type(psb_d_coo_sparse_mat) :: tmpcoo

  call a%mv_to_coo(tmpcoo,info)
  if (info == 0) call tmpcoo%clean_zeros(info)
  if (info == 0) call a%mv_from_coo(tmpcoo,info)

end subroutine psb_d_base_clean_zeros


subroutine psb_d_base_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csput_a
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csput'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_csput_a

subroutine psb_d_base_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csput_v
  use psb_d_base_vect_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_vect_type), intent(inout)  :: val
  class(psb_i_base_vect_type), intent(inout)  :: ia, ja
  integer(psb_ipk_), intent(in)               :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_)  :: err_act, nzin, nzout
  character(len=20)  :: name='csput_v'
  integer :: jmin_, jmax_
  logical :: append_, rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(val%v).and.allocated(ia%v).and.allocated(ja%v)) then
    if (a%is_dev())   call a%sync()
    if (val%is_dev()) call val%sync()
    if (ia%is_dev())  call ia%sync()
    if (ja%is_dev())  call ja%sync()
    call a%csput(nz,ia%v,ja%v,val%v,imin,imax,jmin,jmax,info)
  else
    info = psb_err_invalid_mat_state_
  endif
  if (info /= 0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_csput_v

subroutine psb_d_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csgetrow
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale,chksz
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_csgetrow

!
! Here we have the base implementation of getblk and clip:
! this is just based on the getrow.
! If performance is critical it can be overridden.
!
subroutine psb_d_base_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale,chksz)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csgetblk
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale,chksz
  integer(psb_ipk_)  :: err_act, nzin, nzout
  character(len=20)  :: name='csget'
  integer(psb_ipk_)  :: jmin_, jmax_
  logical :: append_, rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(append)) then
    append_ = append
  else
    append_ = .false.
  endif
  if (append_) then
    nzin = a%get_nzeros()
  else
    nzin = 0
  endif
  if (present(rscale)) then
    rscale_=rscale
  else
    rscale_=.false.
  end if
  if (present(cscale)) then
    cscale_=cscale
  else
    cscale_=.false.
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif

  if (append_.and.(rscale_.or.cscale_)) then
    write(psb_err_unit,*) &
         & 'd_csgetblk: WARNING: dubious input: append_ and rscale_|cscale_'
  end if

  if (rscale_) then
    call b%set_nrows(imax-imin+1)
  else
    call b%set_nrows(max(min(imax,a%get_nrows()),b%get_nrows()))
  end if

  if (cscale_) then
    call b%set_ncols(jmax_-jmin_+1)
  else
    call b%set_ncols(max(min(jmax_,a%get_ncols()),b%get_ncols()))
  end if

  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
       & nzin=nzin, rscale=rscale, cscale=cscale, chksz=chksz)

  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_csgetblk


subroutine psb_d_base_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csclip
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_)  :: err_act, nzin, nzout, imin_, imax_, jmin_, jmax_, mb,nb
  character(len=20)  :: name='csget'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  nzin = 0
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = a%get_nrows() ! Should this be imax_ ??
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = a%get_ncols()  ! Should this be jmax_ ??
  endif
  call b%allocate(mb,nb)
  call a%csget(imin_,imax_,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin_, jmax=jmax_, append=.false., &
       & nzin=nzin, rscale=rscale_, cscale=cscale_)
  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_csclip


!
! Here we have the base implementation of tril and triu
! this is just based on the getrow.
! If performance is critical it can be overridden.
!
subroutine psb_d_base_tril(a,l,info,&
     & diag,imin,imax,jmin,jmax,rscale,cscale,u)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_tril
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(out) :: l
  integer(psb_ipk_),intent(out)            :: info
  integer(psb_ipk_), intent(in), optional  :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional            :: rscale,cscale
  class(psb_d_coo_sparse_mat), optional, intent(out) :: u

  integer(psb_ipk_) :: err_act, nzin, nzout, i, j, k, ibk
  integer(psb_ipk_) :: imin_, imax_, jmin_, jmax_, mb,nb, diag_, nzlin, nzuin, nz
  integer(psb_ipk_), allocatable :: ia(:), ja(:)
  real(psb_dpk_), allocatable    :: val(:)
  character(len=20)  :: name='tril'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.
  integer(psb_ipk_), parameter :: nbk=8

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(diag)) then
    diag_ = diag
  else
    diag_ = 0
  end if
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = imax_
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = jmax_
  endif


  nz = a%get_nzeros()
  call l%allocate(mb,nb,nz)

  if (present(u)) then
    nzlin = l%get_nzeros() ! At this point it should be 0
    call u%allocate(mb,nb,nz)
    nzuin = u%get_nzeros() ! At this point it should be 0
    call psb_realloc(max(mb,nb),ia,info)
    call psb_realloc(max(mb,nb),ja,info)
    call psb_realloc(max(mb,nb),val,info)
    do i=imin_,imax_, nbk
      ibk = min(nbk,imax_-i+1)
      call a%csget(i,i+ibk-1,nzout,ia,ja,val,info,&
           & jmin=jmin_, jmax=jmax_)
      do k=1, nzout
        if ((ja(k)-ia(k))<=diag_) then
          nzlin = nzlin + 1
          l%ia(nzlin)  = ia(k)
          l%ja(nzlin)  = ja(k)
          l%val(nzlin) = val(k)
        else
          nzuin = nzuin + 1
          u%ia(nzuin)  = ia(k)
          u%ja(nzuin)  = ja(k)
          u%val(nzuin) = val(k)
        end if
      end do
    end do

    call l%set_nzeros(nzlin)
    call u%set_nzeros(nzuin)
    call u%fix(info)
    nzout = u%get_nzeros()
    if (rscale_) &
         & u%ia(1:nzout) = u%ia(1:nzout) - imin_ + 1
    if (cscale_) &
         & u%ja(1:nzout) = u%ja(1:nzout) - jmin_ + 1
    if ((diag_ >= -1).and.(imin_ == jmin_)) then
      call u%set_triangle(.true.)
      call u%set_lower(.false.)
    end if
  else
    nzin = l%get_nzeros() ! At this point it should be 0
    do i=imin_,imax_
      k = min(jmax_,i+diag_)
      call a%csget(i,i,nzout,l%ia,l%ja,l%val,info,&
           & jmin=jmin_, jmax=k, append=.true., &
           & nzin=nzin)
      if (info /= psb_success_) goto 9999
      call l%set_nzeros(nzin+nzout)
      nzin = nzin+nzout
    end do
  end if
  call l%fix(info)
  nzout = l%get_nzeros()
  if (rscale_) &
       & l%ia(1:nzout) = l%ia(1:nzout) - imin_ + 1
  if (cscale_) &
       & l%ja(1:nzout) = l%ja(1:nzout) - jmin_ + 1

  if ((diag_ <= 0).and.(imin_ == jmin_)) then
    call l%set_triangle(.true.)
    call l%set_lower(.true.)
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_tril

subroutine psb_d_base_triu(a,u,info,&
     & diag,imin,imax,jmin,jmax,rscale,cscale,l)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_triu
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(out) :: u
  integer(psb_ipk_),intent(out)            :: info
  integer(psb_ipk_), intent(in), optional  :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional            :: rscale,cscale
  class(psb_d_coo_sparse_mat), optional, intent(out) :: l

  integer(psb_ipk_) :: err_act, nzin, nzout, i, j, k, ibk
  integer(psb_ipk_) :: imin_, imax_, jmin_, jmax_, mb,nb, diag_, nzlin, nzuin, nz
  integer(psb_ipk_), allocatable :: ia(:), ja(:)
  real(psb_dpk_), allocatable    :: val(:)
  character(len=20)  :: name='triu'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.
  integer(psb_ipk_), parameter :: nbk=8

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(diag)) then
    diag_ = diag
  else
    diag_ = 0
  end if
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = imax_
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = jmax_
  endif


  nz = a%get_nzeros()
  call u%allocate(mb,nb,nz)

  if (present(l)) then
    nzuin = u%get_nzeros() ! At this point it should be 0
    call l%allocate(mb,nb,nz)
    nzlin = l%get_nzeros() ! At this point it should be 0
    call psb_realloc(max(mb,nb),ia,info)
    call psb_realloc(max(mb,nb),ja,info)
    call psb_realloc(max(mb,nb),val,info)
    do i=imin_,imax_, nbk
      ibk = min(nbk,imax_-i+1)
      call a%csget(i,i+ibk-1,nzout,ia,ja,val,info,&
           & jmin=jmin_, jmax=jmax_)
      do k=1, nzout
        if ((ja(k)-ia(k))<diag_) then
          nzlin = nzlin + 1
          l%ia(nzlin)  = ia(k)
          l%ja(nzlin)  = ja(k)
          l%val(nzlin) = val(k)
        else
          nzuin = nzuin + 1
          u%ia(nzuin)  = ia(k)
          u%ja(nzuin)  = ja(k)
          u%val(nzuin) = val(k)
        end if
      end do
    end do
    call u%set_nzeros(nzuin)
    call l%set_nzeros(nzlin)
    call l%fix(info)
    nzout = l%get_nzeros()
    if (rscale_) &
         & l%ia(1:nzout) = l%ia(1:nzout) - imin_ + 1
    if (cscale_) &
         & l%ja(1:nzout) = l%ja(1:nzout) - jmin_ + 1
    if ((diag_ <=1).and.(imin_ == jmin_)) then
      call l%set_triangle(.true.)
      call l%set_lower(.true.)
    end if
  else
    nzin = u%get_nzeros()
    do i=imin_,imax_
      k = max(jmin_,i+diag_)
      call a%csget(i,i,nzout,u%ia,u%ja,u%val,info,&
           & jmin=k, jmax=jmax_, append=.true., &
           & nzin=nzin)
      if (info /= psb_success_) goto 9999
      call u%set_nzeros(nzin+nzout)
      nzin = nzin+nzout
    end do
  end if
  call u%fix(info)
  nzout = u%get_nzeros()
  if (rscale_) &
       & u%ia(1:nzout) = u%ia(1:nzout) - imin_ + 1
  if (cscale_) &
       & u%ja(1:nzout) = u%ja(1:nzout) - jmin_ + 1

  if ((diag_ >= 0).and.(imin_ == jmin_)) then
    call u%set_triangle(.true.)
    call u%set_upper(.true.)
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_triu



subroutine psb_d_base_clone(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_clone
  use psb_error_mod
  implicit none

  class(psb_d_base_sparse_mat), intent(inout)              :: a
  class(psb_d_base_sparse_mat), allocatable, intent(inout) :: b
  integer(psb_ipk_), intent(out) :: info

  info = 0
  if (allocated(b)) then
    call b%free()
    deallocate(b, stat=info)
  end if
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
    return
  end if

  ! Do not use SOURCE allocation: this makes sure that
  ! memory allocated elsewhere is treated properly.
  allocate(b,mold=a,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call b%cp_from_fmt(a, info)

end subroutine psb_d_base_clone

subroutine psb_d_base_make_nonunit(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_make_nonunit
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  type(psb_d_coo_sparse_mat) :: tmp

  integer(psb_ipk_) :: i, j, m, n, nz, mnm, info

  if (a%is_unit()) then
    call a%mv_to_coo(tmp,info)
    if (info /= 0) return
    m = tmp%get_nrows()
    n = tmp%get_ncols()
    mnm = min(m,n)
    nz = tmp%get_nzeros()
    call tmp%reallocate(nz+mnm)
    do i=1, mnm
      tmp%val(nz+i) = done
      tmp%ia(nz+i)  = i
      tmp%ja(nz+i)  = i
    end do
    call tmp%set_nzeros(nz+mnm)
    call tmp%set_unit(.false.)
    call tmp%fix(info)
    if (info /= 0) &
         & call a%mv_from_coo(tmp,info)
  end if

end subroutine psb_d_base_make_nonunit

subroutine psb_d_base_mold(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mold
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in)                 :: a
  class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='base_mold'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_mold

subroutine psb_d_base_transp_2mat(a,b)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_transp_2mat
  use psb_error_mod
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_base_sparse_mat), intent(out)    :: b

  type(psb_d_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='d_base_transp'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_d_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transp()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=(/ione/))
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_transp_2mat

subroutine psb_d_base_transc_2mat(a,b)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_transc_2mat
  implicit none

  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_base_sparse_mat), intent(out)    :: b

  type(psb_d_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='d_base_transc'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_d_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transc()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=(/ione/))
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_d_base_transc_2mat

subroutine psb_d_base_transp_1mat(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_transp_1mat
  use psb_error_mod
  implicit none

  class(psb_d_base_sparse_mat), intent(inout) :: a

  type(psb_d_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='d_base_transp'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transp()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_missing_override_method_
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_transp_1mat

subroutine psb_d_base_transc_1mat(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_transc_1mat
  implicit none

  class(psb_d_base_sparse_mat), intent(inout) :: a

  type(psb_d_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='d_base_transc'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transc()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_missing_override_method_
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_transc_1mat


! == ==================================
!
!
!
! Computational routines
!
!
!
!
!
!
! == ==================================

subroutine psb_d_base_csmm(alpha,a,x,beta,y,info,trans)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csmm
  use psb_error_mod

  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_base_csmm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_csmm


subroutine psb_d_base_csmv(alpha,a,x,beta,y,info,trans)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csmv
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_base_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)


end subroutine psb_d_base_csmv


subroutine psb_d_base_inner_cssm(alpha,a,x,beta,y,info,trans)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_inner_cssm
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_base_inner_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_inner_cssm


subroutine psb_d_base_inner_cssv(alpha,a,x,beta,y,info,trans)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_inner_cssv
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_base_inner_cssv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_inner_cssv


subroutine psb_d_base_cssm(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cssm
  use psb_error_mod
  use psb_string_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)

  real(psb_dpk_), allocatable :: tmp(:,:)
  integer(psb_ipk_) :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  character(len=20)  :: name='d_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = min(size(x,2), size(y,2))
  if (size(x,1) < nac) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/3_psb_ipk_,nac/))
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,nar/))
    goto 9999
  end if

  if (.not. (a%is_triangle())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(d)) then
    if (present(scale)) then
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then
      if (size(d,1) < nac) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nac/))
        goto 9999
      end if

      allocate(tmp(nac,nc),stat=info)
      if (info /= psb_success_) info = psb_err_alloc_dealloc_
      if (info == psb_success_) then
        do i=1, nac
          tmp(i,1:nc) = d(i)*x(i,1:nc)
        end do
      end if
      if (info == psb_success_)&
           & call a%inner_spsm(alpha,tmp,beta,y,info,trans)

      if (info == psb_success_) then
        deallocate(tmp,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then

      if (size(d,1) < nar) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nar/))
        goto 9999
      end if

      allocate(tmp(nar,nc),stat=info)
      if (info /= psb_success_) info = psb_err_alloc_dealloc_
      if (info == psb_success_)&
           & call a%inner_spsm(done,x,dzero,tmp,info,trans)

      if (info == psb_success_)then
        do i=1, nar
          tmp(i,1:nc) = d(i)*tmp(i,1:nc)
        end do
      end if
      if (info == psb_success_)&
           & call psb_geaxpby(nar,nc,alpha,tmp,beta,y,info)

      if (info == psb_success_) then
        deallocate(tmp,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else
      info = 31
      call psb_errpush(info,name,i_err=(/8_psb_ipk_,izero/),a_err=scale_)
      goto 9999
    end if
  else
    ! Scale is ignored in this case
    call a%inner_spsm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cssm


subroutine psb_d_base_cssv(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cssv
  use psb_error_mod
  use psb_string_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)

  real(psb_dpk_), allocatable :: tmp(:)
  integer(psb_ipk_)  :: err_act, nar,nac,nc, i
  character(len=1)   :: scale_
  character(len=20)  :: name='d_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = 1
  if (size(x,1) < nac) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/3_psb_ipk_,nac/))
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,nar/))
    goto 9999
  end if

  if (.not. (a%is_triangle())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(d)) then
    if (present(scale)) then
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then
      if (size(d,1) < nac) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nac/))
        goto 9999
      end if

      allocate(tmp(nac),stat=info)
      if (info /= psb_success_) info = psb_err_alloc_dealloc_
      if (info == psb_success_) call inner_vscal(nac,d,x,tmp)
      if (info == psb_success_)&
           & call a%inner_spsm(alpha,tmp,beta,y,info,trans)

      if (info == psb_success_) then
        deallocate(tmp,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then
      if (size(d,1) < nar) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nar/))
        goto 9999
      end if

      if (beta == dzero) then
        call a%inner_spsm(alpha,x,dzero,y,info,trans)
        if (info == psb_success_)  call inner_vscal1(nar,d,y)
      else
        allocate(tmp(nar),stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
        if (info == psb_success_)&
             & call a%inner_spsm(alpha,x,dzero,tmp,info,trans)

        if (info == psb_success_)  call inner_vscal1(nar,d,tmp)
        if (info == psb_success_)&
             & call psb_geaxpby(nar,done,tmp,beta,y,info)
        if (info == psb_success_) then
          deallocate(tmp,stat=info)
          if (info /= psb_success_) info = psb_err_alloc_dealloc_
        end if
      end if

    else
      info = 31
      call psb_errpush(info,name,i_err=(/8_psb_ipk_,izero/),a_err=scale_)
      goto 9999
    end if
  else
    ! Scale is ignored in this case
    call a%inner_spsm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='inner_spsm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
contains
  subroutine inner_vscal(n,d,x,y)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    real(psb_dpk_), intent(in)  :: d(*),x(*)
    real(psb_dpk_), intent(out) :: y(*)
    integer(psb_ipk_) :: i

    do i=1,n
      y(i) = d(i)*x(i)
    end do
  end subroutine inner_vscal


  subroutine inner_vscal1(n,d,x)
    implicit none
    integer(psb_ipk_), intent(in)         :: n
    real(psb_dpk_), intent(in)  :: d(*)
    real(psb_dpk_), intent(inout) :: x(*)
    integer(psb_ipk_) :: i

    do i=1,n
      x(i) = d(i)*x(i)
    end do
  end subroutine inner_vscal1

end subroutine psb_d_base_cssv


subroutine psb_d_base_scals(d,a,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_scals
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_scals'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_scals

subroutine psb_d_base_scalplusidentity(d,a,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_scalplusidentity
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_scalplusidentity'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_scalplusidentity

subroutine psb_d_base_scal(d,a,info,side)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_scal
  use psb_error_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer(psb_ipk_), intent(out)            :: info
  character, intent(in), optional :: side

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_scal'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_scal

function psb_d_base_maxval(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_maxval

  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='maxval'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  res = dzero
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end function psb_d_base_maxval


function psb_d_base_csnmi(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csnmi

  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='csnmi'
  real(psb_dpk_), allocatable  :: vt(:)

  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  res = dzero
  call psb_realloc(a%get_nrows(),vt,info)
  if (info /= 0) then
    info  = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%arwsum(vt)
  res = maxval(vt)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end function psb_d_base_csnmi

function psb_d_base_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_csnm1

  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='csnm1'
  real(psb_dpk_), allocatable  :: vt(:)

  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  res = dzero
  call psb_realloc(a%get_ncols(),vt,info)
  if (info /= 0) then
    info  = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%aclsum(vt)
  res = maxval(vt)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end function psb_d_base_csnm1

subroutine psb_d_base_rowsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_rowsum
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_rowsum

subroutine psb_d_base_arwsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_arwsum
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_arwsum

subroutine psb_d_base_colsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_colsum
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_colsum

subroutine psb_d_base_aclsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_aclsum
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_aclsum

subroutine psb_d_base_get_diag(a,d,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_get_diag

  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_d_base_get_diag

subroutine psb_d_base_spaxpby(alpha,a,beta,b,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_spaxpby

  real(psb_dpk_), intent(in)                   :: alpha
  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)                   :: beta
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                :: info

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='spaxpby'
  logical, parameter           :: debug=.false.
  type(psb_d_coo_sparse_mat) :: acoo

  call psb_erractionsave(err_act)
  if((a%get_ncols() /= b%get_ncols()).or.(a%get_nrows() /= b%get_nrows())) then
    info  = psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  call a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_coo')
    goto 9999
  end if

  call acoo%spaxpby(alpha,beta,b,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='spaxby')
    goto 9999
  end if

  call acoo%mv_to_fmt(a,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_fmt')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_d_base_spaxpby

function psb_d_base_cmpval(a,val,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cmpval

  class(psb_d_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)                   :: val
  real(psb_dpk_), intent(in)                  :: tol
  integer(psb_ipk_), intent(out)                :: info
  logical                                       :: res

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='cmpval'
  logical, parameter           :: debug=.false.
  type(psb_d_coo_sparse_mat) :: acoo

  call a%cp_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cp_to_coo')
    goto 9999
  end if

  res = acoo%spcmp(val,tol,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cmpval')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end function psb_d_base_cmpval

function psb_d_base_cmpmat(a,b,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cmpmat

  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  real(psb_dpk_), intent(in)                  :: tol
  integer(psb_ipk_), intent(out)                :: info
  logical                                       :: res

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='cmpmat'
  logical, parameter           :: debug=.false.
  type(psb_d_coo_sparse_mat) :: acoo

  call a%cp_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cp_to_coo')
    goto 9999
  end if

  ! Fix the indexes
  call acoo%fix(info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='fix')
    goto 9999
  end if

  res = acoo%spcmp(b,tol,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cmpmat')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end function psb_d_base_cmpmat

! == ==================================
!
!
!
! Computational routines for d_VECT
! variables. If the actual data type is
! a "normal" one, these are sufficient.
!
!
!
!
! == ==================================



subroutine psb_d_base_vect_mv(alpha,a,x,beta,y,info,trans)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_vect_mv
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_d_base_vect_type), intent(inout) :: x
  class(psb_d_base_vect_type), intent(inout) :: y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans

  ! For the time being we just throw everything back
  ! onto the normal routines.
  call x%sync()
  call y%sync()
  call a%spmm(alpha,x%v,beta,y%v,info,trans)
  call y%set_host()
end subroutine psb_d_base_vect_mv

subroutine psb_d_base_vect_cssv(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_vect_cssv
  use psb_d_base_vect_mod
  use psb_error_mod
  use psb_string_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_d_base_vect_type), intent(inout) :: x,y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans, scale
  class(psb_d_base_vect_type), intent(inout),optional  :: d

  real(psb_dpk_), allocatable :: tmp(:)
  class(psb_d_base_vect_type), allocatable :: tmpv
  integer(psb_ipk_)  :: err_act, nar,nac,nc, i
  character(len=1)   :: scale_
  character(len=20)  :: name='d_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = 1
  if (x%get_nrows() < nac) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/3_psb_ipk_,nac/))
    goto 9999
  end if
  if (y%get_nrows() < nar) then
    info = psb_err_input_asize_small_i_
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,nar/))
    goto 9999
  end if

  if (.not. (a%is_triangle())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call x%sync()
  call y%sync()
  if (present(d)) then
    call d%sync()
    if (present(scale)) then
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then
      if (d%get_nrows() < nac) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nac/))
        goto 9999
      end if
      allocate(tmpv, mold=y,stat=info)
      if (info /= psb_success_) info = psb_err_alloc_dealloc_
      if (info == psb_success_) call tmpv%mlt(done,d%v(1:nac),x,dzero,info)
      if (info == psb_success_)&
           & call a%inner_spsm(alpha,tmpv,beta,y,info,trans)

      if (info == psb_success_) then
        call tmpv%free(info)
        if (info == psb_success_) deallocate(tmpv,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then
      if (d%get_nrows() < nar) then
        info = psb_err_input_asize_small_i_
        call psb_errpush(info,name,i_err=(/9_psb_ipk_,nar/))
        goto 9999
      end if

      if (beta == dzero) then
        call a%inner_spsm(alpha,x,dzero,y,info,trans)
        if (info == psb_success_)  call y%mlt(d%v(1:nar),info)

      else
        allocate(tmpv, mold=y,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
        if (info == psb_success_)&
             & call a%inner_spsm(alpha,x,dzero,tmpv,info,trans)

        if (info == psb_success_)  call tmpv%mlt(d%v(1:nar),info)
        if (info == psb_success_)&
             & call y%axpby(nar,done,tmpv,beta,info)
        if (info == psb_success_) then
          call tmpv%free(info)
          if (info == psb_success_) deallocate(tmpv,stat=info)
          if (info /= psb_success_) info = psb_err_alloc_dealloc_
        end if
      end if

    else
      info = 31
      call psb_errpush(info,name,i_err=(/8_psb_ipk_,izero/),a_err=scale_)
      goto 9999
    end if
  else
    ! Scale is ignored in this case
    call a%inner_spsm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='inner_spsm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_vect_cssv


subroutine psb_d_base_inner_vect_sv(alpha,a,x,beta,y,info,trans)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_inner_vect_sv
  use psb_error_mod
  use psb_string_mod
  use psb_d_base_vect_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_d_base_vect_type), intent(inout) :: x, y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_base_inner_vect_sv'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%inner_spsm(alpha,x%v,beta,y%v,info,trans)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='inner_spsm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_inner_vect_sv


subroutine psb_d_base_cp_to_lcoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_to_lcoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_lcoo'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
   !
  info = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_coo(tmp,info)
  if (info == psb_success_) call tmp%mv_to_lcoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_to_lcoo

subroutine psb_d_base_cp_from_lcoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_from_lcoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  call tmp%cp_from_lcoo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_from_lcoo

subroutine psb_d_base_cp_to_lfmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_to_lfmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(in) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_lfmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%cp_to_lcoo(b,info)
  class default
    call a%cp_to_coo(icoo,info)
    if (info == psb_success_) call icoo%mv_to_lcoo(lcoo,info)
    if (info == psb_success_) call b%mv_from_coo(lcoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_to_lfmt

subroutine psb_d_base_cp_from_lfmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_cp_from_lfmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_lfmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo
  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%cp_from_lcoo(b,info)
  class default
    call b%cp_to_coo(lcoo,info)
    if (info == psb_success_) call lcoo%mv_to_icoo(icoo,info)
    if (info == psb_success_) call a%mv_from_coo(icoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_cp_from_lfmt


subroutine psb_d_base_mv_to_lcoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_to_lcoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_lcoo'
  logical, parameter :: debug=.false.


  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_lcoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to coo')
    goto 9999
  end if

  call a%free()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_to_lcoo

subroutine psb_d_base_mv_from_lcoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_from_lcoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_lcoo'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_from_lcoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='from coo')
    goto 9999
  end if

  call b%free()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_from_lcoo


subroutine psb_d_base_mv_to_lfmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_to_lfmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_lfmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%mv_to_lcoo(b,info)
  class default
    call a%mv_to_coo(icoo,info)
    if (info == psb_success_) call icoo%mv_to_lcoo(lcoo,info)
    if (info == psb_success_) call b%mv_from_coo(lcoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_to_lfmt

subroutine psb_d_base_mv_from_lfmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_d_base_mv_from_lfmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_d_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_lfmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo
  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%mv_from_lcoo(b,info)
  class default
    call b%mv_to_coo(lcoo,info)
    if (info == psb_success_) call lcoo%mv_to_icoo(icoo,info)
    if (info == psb_success_) call a%mv_from_coo(icoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_base_mv_from_lfmt

!
!
! ld implementation
!
!
! == ==================================
!
!
!
! Data management
!
!
!
!
!
! == ==================================

subroutine psb_ld_base_cp_to_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_cp_to_coo

subroutine psb_ld_base_cp_from_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_cp_from_coo


subroutine psb_ld_base_cp_to_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%cp_to_coo(b,info)
  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_to_fmt

subroutine psb_ld_base_cp_from_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%cp_from_coo(b,info)
  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_from_fmt


subroutine psb_ld_base_mv_to_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_coo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to coo')
    goto 9999
  end if

  call a%free()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_to_coo

subroutine psb_ld_base_mv_from_coo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_from_coo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='from coo')
    goto 9999
  end if

  call b%free()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_from_coo


subroutine psb_ld_base_mv_to_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%mv_to_coo(b,info)
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

  return

end subroutine psb_ld_base_mv_to_fmt

subroutine psb_ld_base_mv_from_fmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  !
  info = psb_success_
  select type(b)
  type is (psb_ld_coo_sparse_mat)
    call a%mv_from_coo(b,info)
  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
  return

end subroutine psb_ld_base_mv_from_fmt

subroutine  psb_ld_base_clean_zeros(a, info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_clean_zeros
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)              :: info
  !
  type(psb_ld_coo_sparse_mat) :: tmpcoo

  call a%mv_to_coo(tmpcoo,info)
  if (info == 0) call tmpcoo%clean_zeros(info)
  if (info == 0) call a%mv_from_coo(tmpcoo,info)

end subroutine psb_ld_base_clean_zeros


subroutine psb_ld_base_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csput_a
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_lpk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csput'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_csput_a

subroutine psb_ld_base_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_error_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csput_v
  use psb_d_base_vect_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_vect_type), intent(inout)  :: val
  class(psb_l_base_vect_type), intent(inout)  :: ia, ja
  integer(psb_lpk_), intent(in)               :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_lpk_)  :: nzin, nzout
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csput_v'
  integer :: jmin_, jmax_
  logical :: append_, rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (allocated(val%v).and.allocated(ia%v).and.allocated(ja%v)) then
    if (a%is_dev())   call a%sync()
    if (val%is_dev()) call val%sync()
    if (ia%is_dev())  call ia%sync()
    if (ja%is_dev())  call ja%sync()
    call a%csput_a(nz,ia%v,ja%v,val%v,imin,imax,jmin,jmax,info)
  else
    info = psb_err_invalid_mat_state_
  endif
  if (info /= 0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_csput_v

subroutine psb_ld_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csgetrow
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_lpk_), intent(out)                 :: nz
  integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_lpk_), intent(in), optional        :: iren(:)
  integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_csgetrow



!
! Here we have the base implementation of getblk and clip:
! this is just based on the getrow.
! If performance is critical it can be overridden.
!
subroutine psb_ld_base_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csgetblk
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(inout) :: b
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_lpk_), intent(in), optional        :: iren(:)
  integer(psb_lpk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_)  :: err_act
  integer(psb_lpk_)  :: nzin, nzout
  character(len=20)  :: name='csget'
  integer(psb_lpk_)  :: jmin_, jmax_
  logical :: append_, rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(append)) then
    append_ = append
  else
    append_ = .false.
  endif
  if (append_) then
    nzin = a%get_nzeros()
  else
    nzin = 0
  endif
  if (present(rscale)) then
    rscale_=rscale
  else
    rscale_=.false.
  end if
  if (present(cscale)) then
    cscale_=cscale
  else
    cscale_=.false.
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif

  if (append_.and.(rscale_.or.cscale_)) then
    write(psb_err_unit,*) &
         & 'ld_csgetblk: WARNING: dubious input: append_ and rscale_|cscale_'
  end if

  if (rscale_) then
    call b%set_nrows(imax-imin+1)
  else
    call b%set_nrows(max(min(imax,a%get_nrows()),b%get_nrows()))
  end if

  if (cscale_) then
    call b%set_ncols(jmax_-jmin_+1)
  else
    call b%set_ncols(max(min(jmax_,a%get_ncols()),b%get_ncols()))
  end if

  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
       & nzin=nzin, rscale=rscale, cscale=cscale)

  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_csgetblk


subroutine psb_ld_base_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csclip
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_lpk_)  :: nzin, nzout, imin_, imax_, jmin_, jmax_, mb,nb
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='csget'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  nzin = 0
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = a%get_nrows() ! Should this be imax_ ??
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = a%get_ncols()  ! Should this be jmax_ ??
  endif
  call b%allocate(mb,nb)
  call a%csget(imin_,imax_,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin_, jmax=jmax_, append=.false., &
       & nzin=nzin, rscale=rscale_, cscale=cscale_)
  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_csclip


!
! Here we have the base implementation of tril and triu
! this is just based on the getrow.
! If performance is critical it can be overridden.
!
subroutine psb_ld_base_tril(a,l,info,&
     & diag,imin,imax,jmin,jmax,rscale,cscale,u)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_tril
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(out) :: l
  integer(psb_ipk_),intent(out)            :: info
  integer(psb_lpk_), intent(in), optional  :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional            :: rscale,cscale
  class(psb_ld_coo_sparse_mat), optional, intent(out) :: u

  integer(psb_ipk_) :: err_act
  integer(psb_lpk_) :: nzin, nzout, i, j, k, ibk
  integer(psb_lpk_) :: imin_, imax_, jmin_, jmax_, mb,nb, diag_, nzlin, nzuin, nz
  integer(psb_lpk_), allocatable :: ia(:), ja(:)
  real(psb_dpk_), allocatable    :: val(:)
  character(len=20)  :: name='tril'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.
  integer(psb_lpk_), parameter :: nbk=8

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(diag)) then
    diag_ = diag
  else
    diag_ = 0
  end if
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = imax_
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = jmax_
  endif


  nz = a%get_nzeros()
  call l%allocate(mb,nb,nz)

  if (present(u)) then
    nzlin = l%get_nzeros() ! At this point it should be 0
    call u%allocate(mb,nb,nz)
    nzuin = u%get_nzeros() ! At this point it should be 0
    call psb_realloc(max(mb,nb),ia,info)
    call psb_realloc(max(mb,nb),ja,info)
    call psb_realloc(max(mb,nb),val,info)
    do i=imin_,imax_, nbk
      ibk = min(nbk,imax_-i+1)
      call a%csget(i,i+ibk-1,nzout,ia,ja,val,info,&
           & jmin=jmin_, jmax=jmax_)
      do k=1, nzout
        if ((ja(k)-ia(k))<=diag_) then
          nzlin = nzlin + 1
          l%ia(nzlin)  = ia(k)
          l%ja(nzlin)  = ja(k)
          l%val(nzlin) = val(k)
        else
          nzuin = nzuin + 1
          u%ia(nzuin)  = ia(k)
          u%ja(nzuin)  = ja(k)
          u%val(nzuin) = val(k)
        end if
      end do
    end do

    call l%set_nzeros(nzlin)
    call u%set_nzeros(nzuin)
    call u%fix(info)
    nzout = u%get_nzeros()
    if (rscale_) &
         & u%ia(1:nzout) = u%ia(1:nzout) - imin_ + 1
    if (cscale_) &
         & u%ja(1:nzout) = u%ja(1:nzout) - jmin_ + 1
    if ((diag_ >= -1).and.(imin_ == jmin_)) then
      call u%set_triangle(.true.)
      call u%set_lower(.false.)
    end if
  else
    nzin = l%get_nzeros() ! At this point it should be 0
    do i=imin_,imax_
      k = min(jmax_,i+diag_)
      call a%csget(i,i,nzout,l%ia,l%ja,l%val,info,&
           & jmin=jmin_, jmax=k, append=.true., &
           & nzin=nzin)
      if (info /= psb_success_) goto 9999
      call l%set_nzeros(nzin+nzout)
      nzin = nzin+nzout
    end do
  end if
  call l%fix(info)
  nzout = l%get_nzeros()
  if (rscale_) &
       & l%ia(1:nzout) = l%ia(1:nzout) - imin_ + 1
  if (cscale_) &
       & l%ja(1:nzout) = l%ja(1:nzout) - jmin_ + 1

  if ((diag_ <= 0).and.(imin_ == jmin_)) then
    call l%set_triangle(.true.)
    call l%set_lower(.true.)
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_tril

subroutine psb_ld_base_triu(a,u,info,&
     & diag,imin,imax,jmin,jmax,rscale,cscale,l)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_triu
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_ld_coo_sparse_mat), intent(out) :: u
  integer(psb_ipk_),intent(out)            :: info
  integer(psb_lpk_), intent(in), optional  :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional            :: rscale,cscale
  class(psb_ld_coo_sparse_mat), optional, intent(out) :: l

  integer(psb_ipk_) :: err_act
  integer(psb_lpk_) :: nzin, nzout, i, j, k, ibk
  integer(psb_lpk_) :: imin_, imax_, jmin_, jmax_, mb,nb, diag_, nzlin, nzuin, nz
  integer(psb_lpk_), allocatable :: ia(:), ja(:)
  real(psb_dpk_), allocatable    :: val(:)
  character(len=20)  :: name='triu'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.
  integer(psb_lpk_), parameter :: nbk=8

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(diag)) then
    diag_ = diag
  else
    diag_ = 0
  end if
  if (present(imin)) then
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then
    mb = imax_ - imin_ +1
  else
    mb = imax_
  endif
  if (cscale_) then
    nb = jmax_ - jmin_ +1
  else
    nb = jmax_
  endif


  nz = a%get_nzeros()
  call u%allocate(mb,nb,nz)

  if (present(l)) then
    nzuin = u%get_nzeros() ! At this point it should be 0
    call l%allocate(mb,nb,nz)
    nzlin = l%get_nzeros() ! At this point it should be 0
    call psb_realloc(max(mb,nb),ia,info)
    call psb_realloc(max(mb,nb),ja,info)
    call psb_realloc(max(mb,nb),val,info)
    do i=imin_,imax_, nbk
      ibk = min(nbk,imax_-i+1)
      call a%csget(i,i+ibk-1,nzout,ia,ja,val,info,&
           & jmin=jmin_, jmax=jmax_)
      do k=1, nzout
        if ((ja(k)-ia(k))<diag_) then
          nzlin = nzlin + 1
          l%ia(nzlin)  = ia(k)
          l%ja(nzlin)  = ja(k)
          l%val(nzlin) = val(k)
        else
          nzuin = nzuin + 1
          u%ia(nzuin)  = ia(k)
          u%ja(nzuin)  = ja(k)
          u%val(nzuin) = val(k)
        end if
      end do
    end do
    call u%set_nzeros(nzuin)
    call l%set_nzeros(nzlin)
    call l%fix(info)
    nzout = l%get_nzeros()
    if (rscale_) &
         & l%ia(1:nzout) = l%ia(1:nzout) - imin_ + 1
    if (cscale_) &
         & l%ja(1:nzout) = l%ja(1:nzout) - jmin_ + 1
    if ((diag_ <=1).and.(imin_ == jmin_)) then
      call l%set_triangle(.true.)
      call l%set_lower(.true.)
    end if
  else
    nzin = u%get_nzeros()
    do i=imin_,imax_
      k = max(jmin_,i+diag_)
      call a%csget(i,i,nzout,u%ia,u%ja,u%val,info,&
           & jmin=k, jmax=jmax_, append=.true., &
           & nzin=nzin)
      if (info /= psb_success_) goto 9999
      call u%set_nzeros(nzin+nzout)
      nzin = nzin+nzout
    end do
  end if
  call u%fix(info)
  nzout = u%get_nzeros()
  if (rscale_) &
       & u%ia(1:nzout) = u%ia(1:nzout) - imin_ + 1
  if (cscale_) &
       & u%ja(1:nzout) = u%ja(1:nzout) - jmin_ + 1

  if ((diag_ >= 0).and.(imin_ == jmin_)) then
    call u%set_triangle(.true.)
    call u%set_upper(.true.)
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_triu



subroutine psb_ld_base_clone(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_clone
  use psb_error_mod
  implicit none

  class(psb_ld_base_sparse_mat), intent(inout)              :: a
  class(psb_ld_base_sparse_mat), allocatable, intent(inout) :: b
  integer(psb_ipk_), intent(out) :: info

  info = 0
  if (allocated(b)) then
    call b%free()
    deallocate(b, stat=info)
  end if
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
    return
  end if

  ! Do not use SOURCE allocation: this makes sure that
  ! memory allocated elsewhere is treated properly.
  allocate(b,mold=a,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call b%cp_from_fmt(a, info)

end subroutine psb_ld_base_clone

subroutine psb_ld_base_make_nonunit(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_make_nonunit
  use psb_error_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  type(psb_ld_coo_sparse_mat) :: tmp

  integer(psb_ipk_) :: info
  integer(psb_lpk_) :: i, j, m, n, nz, mnm

  if (a%is_unit()) then
    call a%mv_to_coo(tmp,info)
    if (info /= 0) return
    m = tmp%get_nrows()
    n = tmp%get_ncols()
    mnm = min(m,n)
    nz = tmp%get_nzeros()
    call tmp%reallocate(nz+mnm)
    do i=1, mnm
      tmp%val(nz+i) = done
      tmp%ia(nz+i)  = i
      tmp%ja(nz+i)  = i
    end do
    call tmp%set_nzeros(nz+mnm)
    call tmp%set_unit(.false.)
    call tmp%fix(info)
    if (info /= 0) &
         & call a%mv_from_coo(tmp,info)
  end if

end subroutine psb_ld_base_make_nonunit

subroutine psb_ld_base_mold(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mold
  use psb_error_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(in)                 :: a
  class(psb_ld_base_sparse_mat), intent(inout), allocatable :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='base_mold'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_mold

subroutine psb_ld_base_transp_2mat(a,b)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_transp_2mat
  use psb_error_mod
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_lbase_sparse_mat), intent(out)    :: b

  type(psb_ld_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='ld_base_transp'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_ld_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transp()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=(/ione/))
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_transp_2mat

subroutine psb_ld_base_transc_2mat(a,b)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_transc_2mat
  implicit none

  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_lbase_sparse_mat), intent(out)    :: b

  type(psb_ld_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='ld_base_transc'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_ld_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transc()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=(/ione/))
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_ld_base_transc_2mat

subroutine psb_ld_base_transp_1mat(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_transp_1mat
  use psb_error_mod
  implicit none

  class(psb_ld_base_sparse_mat), intent(inout) :: a

  type(psb_ld_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='ld_base_transp'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transp()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_missing_override_method_
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_transp_1mat

subroutine psb_ld_base_transc_1mat(a)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_transc_1mat
  implicit none

  class(psb_ld_base_sparse_mat), intent(inout) :: a

  type(psb_ld_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  character(len=*), parameter :: name='ld_base_transc'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transc()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_missing_override_method_
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_transc_1mat

subroutine psb_ld_base_scals(d,a,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_scals
  use psb_error_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='ld_scals'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_scals

subroutine psb_ld_base_scalplusidentity(d,a,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_scalplusidentity
  use psb_error_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='ld_scalplusidentity'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_scalplusidentity

subroutine psb_ld_base_scal(d,a,info,side)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_scal
  use psb_error_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer(psb_ipk_), intent(out)            :: info
  character, intent(in), optional :: side

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='ld_scal'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_scal

function psb_ld_base_maxval(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_maxval

  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='maxval'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  res = dzero
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end function psb_ld_base_maxval

function psb_ld_base_csnmi(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csnmi

  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='csnmi'
  real(psb_dpk_), allocatable  :: vt(:)

  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  res = dzero
  call psb_realloc(a%get_nrows(),vt,info)
  if (info /= 0) then
    info  = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%arwsum(vt)
  res = maxval(vt)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_ld_base_csnmi

function psb_ld_base_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_csnm1

  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='csnm1'
  real(psb_dpk_), allocatable  :: vt(:)

  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  res = dzero
  call psb_realloc(a%get_ncols(),vt,info)
  if (info /= 0) then
    info  = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%aclsum(vt)
  res = maxval(vt)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_ld_base_csnm1

subroutine psb_ld_base_rowsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_rowsum
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_rowsum

subroutine psb_ld_base_arwsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_arwsum
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_arwsum

subroutine psb_ld_base_colsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_colsum
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_colsum

subroutine psb_ld_base_aclsum(d,a)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_aclsum
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_)  :: err_act, info
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_aclsum

subroutine psb_ld_base_spaxpby(alpha,a,beta,b,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_spaxpby

  real(psb_dpk_), intent(in)                   :: alpha
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)                   :: beta
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                :: info

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='spaxpby'
  logical, parameter           :: debug=.false.
  type(psb_ld_coo_sparse_mat) :: acoo

  call psb_erractionsave(err_act)
  if((a%get_ncols() /= b%get_ncols()).or.(a%get_nrows() /= b%get_nrows())) then
    info  = psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  end if

  call a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_coo')
    goto 9999
  end if

  call acoo%spaxpby(alpha,beta,b,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='spaxby')
    goto 9999
  end if

  call acoo%mv_to_fmt(a,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_fmt')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_ld_base_spaxpby

function psb_ld_base_cmpval(a,val,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cmpval

  class(psb_ld_base_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)                   :: val
  real(psb_dpk_), intent(in)                  :: tol
  integer(psb_ipk_), intent(out)                :: info
  logical                                       :: res

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='cmpval'
  logical, parameter           :: debug=.false.
  type(psb_ld_coo_sparse_mat) :: acoo

  call a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_coo')
    goto 9999
  end if

  res = acoo%spcmp(val,tol,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cmpval')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end function psb_ld_base_cmpval

function psb_ld_base_cmpmat(a,b,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cmpmat

  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  real(psb_dpk_), intent(in)                  :: tol
  integer(psb_ipk_), intent(out)                :: info
  logical                                       :: res

  ! Auxiliary
  integer(psb_ipk_)            :: err_act
  character(len=20)            :: name='cmpmat'
  logical, parameter           :: debug=.false.
  type(psb_ld_coo_sparse_mat) :: acoo

  call a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='mv_to_coo')
    goto 9999
  end if

  ! Fix the indexes
  call acoo%fix(info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='fix')
    goto 9999
  end if

  res = acoo%spcmp(b,tol,info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='cmpmat')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end function psb_ld_base_cmpmat

subroutine psb_ld_base_get_diag(a,d,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_get_diag

  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  call psb_error_handler(err_act)

end subroutine psb_ld_base_get_diag



subroutine psb_ld_base_cp_to_icoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_to_icoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
   !
  info = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_coo(tmp,info)
  if (info == psb_success_) call tmp%mv_to_icoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_to_icoo

subroutine psb_ld_base_cp_from_icoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_from_icoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_icoo'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
   !
  info = psb_success_
  call psb_erractionsave(err_act)

  call tmp%cp_from_icoo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_from_icoo


subroutine psb_ld_base_cp_to_ifmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_to_ifmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(in) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_ifmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%cp_to_icoo(b,info)
  class default
    call a%cp_to_coo(lcoo,info)
    if (info == psb_success_) call lcoo%mv_to_icoo(icoo,info)
    if (info == psb_success_) call b%mv_from_coo(icoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_to_ifmt

subroutine psb_ld_base_cp_from_ifmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_cp_from_ifmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_ifmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo
  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%cp_from_icoo(b,info)
  class default
    call b%cp_to_coo(icoo,info)
    if (info == psb_success_) call icoo%mv_to_lcoo(lcoo,info)
    if (info == psb_success_) call a%mv_from_coo(lcoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_cp_from_ifmt


subroutine psb_ld_base_mv_to_icoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_to_icoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_icoo'
  logical, parameter :: debug=.false.


  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_icoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to coo')
    goto 9999
  end if

  call a%free()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_to_icoo

subroutine psb_ld_base_mv_from_icoo(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_from_icoo
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_icoo'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_from_icoo(b,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='from coo')
    goto 9999
  end if

  call b%free()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_from_icoo


subroutine psb_ld_base_mv_to_ifmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_to_ifmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='to_ifmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo

  !
  ! Default implementation
  !
  info = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%mv_to_icoo(b,info)
  class default
    call a%mv_to_coo(lcoo,info)
    if (info == psb_success_) call lcoo%mv_to_icoo(icoo,info)
    if (info == psb_success_) call b%mv_from_coo(icoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_to_ifmt

subroutine psb_ld_base_mv_from_ifmt(a,b,info)
  use psb_d_base_mat_mod, psb_protect_name => psb_ld_base_mv_from_ifmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none
  class(psb_ld_base_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='from_ifmt'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat)  :: icoo
  type(psb_ld_coo_sparse_mat) :: lcoo
  !
  ! Default implementation
  !
  info  = psb_success_
  call psb_erractionsave(err_act)

  select type(b)
  type is (psb_d_coo_sparse_mat)
    call a%mv_from_icoo(b,info)
  class default
    call b%mv_to_coo(icoo,info)
    if (info == psb_success_) call icoo%mv_to_lcoo(lcoo,info)
    if (info == psb_success_) call a%mv_from_coo(lcoo,info)
  end select

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_base_mv_from_ifmt
