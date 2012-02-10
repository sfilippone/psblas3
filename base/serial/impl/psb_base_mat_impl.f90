function psb_base_get_nz_row(idx,a) result(res)
  use psb_error_mod
  use psb_base_mat_mod, psb_protect_name => psb_base_get_nz_row
  implicit none 
  integer(psb_ipk_), intent(in)                    :: idx
  class(psb_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_) :: res

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='base_get_nz_row'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  res = -1
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end function psb_base_get_nz_row

function psb_base_get_nzeros(a) result(res)
  use psb_base_mat_mod, psb_protect_name => psb_base_get_nzeros
  use psb_error_mod
  implicit none 
  class(psb_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_) :: res

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='base_get_nzeros'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  res = -1
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end function psb_base_get_nzeros

function psb_base_get_size(a) result(res)
  use psb_base_mat_mod, psb_protect_name => psb_base_get_size
  use psb_error_mod
  implicit none 
  class(psb_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_) :: res

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='get_size'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  res = -1
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end function psb_base_get_size

subroutine psb_base_reinit(a,clear)
  use psb_base_mat_mod, psb_protect_name => psb_base_reinit
  use psb_error_mod
  implicit none 

  class(psb_base_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reinit'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  info = psb_err_missing_override_method_
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_reinit

subroutine psb_base_sparse_print(iout,a,iv,head,ivr,ivc)
  use psb_base_mat_mod, psb_protect_name => psb_base_sparse_print
  use psb_error_mod
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_base_sparse_mat), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  info = psb_err_missing_override_method_
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_sparse_print

subroutine psb_base_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_base_mat_mod, psb_protect_name => psb_base_csgetptn
  implicit none

  class(psb_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_csgetptn

subroutine psb_base_get_neigh(a,idx,neigh,n,info,lev)
  use psb_base_mat_mod, psb_protect_name => psb_base_get_neigh
  use psb_error_mod
  use psb_realloc_mod
  use psb_sort_mod
  implicit none 
  class(psb_base_sparse_mat), intent(in) :: a   
  integer(psb_ipk_), intent(in)                :: idx 
  integer(psb_ipk_), intent(out)               :: n   
  integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), optional, intent(in)      :: lev 

  integer(psb_ipk_) :: lev_, i, nl, ifl,ill,&
       &  n1, err_act, nn, nidx,ntl,ma
  integer(psb_ipk_), allocatable :: ia(:), ja(:)
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if(present(lev)) then 
    lev_ = lev
  else
    lev_=1
  end if
  ! Turns out we can write get_neigh at this
  ! level 
  n = 0
  ma = a%get_nrows()
  call a%csget(idx,idx,n,ia,ja,info)
  if (info == psb_success_) call psb_realloc(n,neigh,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if
  neigh(1:n) = ja(1:n)
  ifl = 1
  ill = n
  do nl = 2, lev_ 
    n1 = ill - ifl + 1
    call psb_ensure_size(ill+n1*n1,neigh,info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
    ntl = 0
    do i=ifl,ill
      nidx=neigh(i)
      if ((nidx /= idx).and.(nidx > 0).and.(nidx <= ma)) then
        call a%csget(nidx,nidx,nn,ia,ja,info)
        if (info == psb_success_) call psb_ensure_size(ill+ntl+nn,neigh,info)
        if (info /= psb_success_) then 
          call psb_errpush(psb_err_alloc_dealloc_,name)
          goto 9999
        end if
        neigh(ill+ntl+1:ill+ntl+nn)=ja(1:nn)
        ntl = ntl+nn
      end if
    end do
    call psb_msort_unique(neigh(ill+1:ill+ntl),nn)
    ifl = ill + 1
    ill = ill + nn
  end do
  call psb_msort_unique(neigh(1:ill),nn,dir=psb_sort_up_)
  n = nn

  call psb_erractionrestore(err_act)
  return

9999 continue

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_get_neigh

subroutine  psb_base_allocate_mnnz(m,n,a,nz) 
  use psb_base_mat_mod, psb_protect_name => psb_base_allocate_mnnz
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in) :: m,n
  class(psb_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(in), optional  :: nz
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_allocate_mnnz

subroutine  psb_base_reallocate_nz(nz,a) 
  use psb_base_mat_mod, psb_protect_name => psb_base_reallocate_nz
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_reallocate_nz

subroutine  psb_base_free(a) 
  use psb_base_mat_mod, psb_protect_name => psb_base_free
  use psb_error_mod
  implicit none 
  class(psb_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='free'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_free

subroutine  psb_base_trim(a) 
  use psb_base_mat_mod, psb_protect_name => psb_base_trim
  use psb_error_mod
  implicit none 
  class(psb_base_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  call psb_errpush(psb_err_missing_override_method_,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_base_trim

