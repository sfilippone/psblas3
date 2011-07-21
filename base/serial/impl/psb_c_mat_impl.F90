! == ===================================
!
!
!
! Setters 
!
!
!
!
!
!
! == ===================================


subroutine  psb_c_set_nrows(m,a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_nrows
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  integer, intent(in) :: m
  Integer :: err_act, info
  character(len=20)  :: name='set_nrows'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_nrows(m)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_nrows


subroutine  psb_c_set_ncols(n,a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_ncols
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call a%a%set_ncols(n)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_ncols



subroutine  psb_c_set_state(n,a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_state
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call a%a%set_state(n)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_state



subroutine  psb_c_set_dupl(n,a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_dupl
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_dupl(n)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_dupl


subroutine  psb_c_set_null(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_null
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_null()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_null


subroutine  psb_c_set_bld(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_bld
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_bld()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_bld


subroutine  psb_c_set_upd(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_upd
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upd()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if


end subroutine psb_c_set_upd


subroutine  psb_c_set_asb(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_asb
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_asb()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_asb


subroutine psb_c_set_sorted(a,val) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_sorted
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_sorted(val)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_sorted


subroutine psb_c_set_triangle(a,val) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_triangle
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_triangle(val)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_triangle


subroutine psb_c_set_unit(a,val) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_unit
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_unit(val)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_unit


subroutine psb_c_set_lower(a,val) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_lower
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_lower(val)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_lower


subroutine psb_c_set_upper(a,val) 
  use psb_c_mat_mod, psb_protect_name => psb_c_set_upper
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upper(val)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_set_upper



! == ===================================
!
!
!
! Data management
!
!
!
!
!
! == ===================================  


subroutine psb_c_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_c_sparse_print
  use psb_error_mod
  implicit none 

  integer, intent(in)               :: iout
  class(psb_cspmat_type), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%print(iout,iv,eirs,eics,head,ivr,ivc)

  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_sparse_print


subroutine psb_c_n_sparse_print(fname,a,iv,eirs,eics,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_c_n_sparse_print
  use psb_error_mod
  implicit none 

  character(len=*), intent(in)  :: fname   
  class(psb_cspmat_type), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act, info, iout
  logical :: isopen
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  iout = max(psb_inp_unit,psb_err_unit,psb_out_unit) + 1
  do 
    inquire(unit=iout, opened=isopen)
    if (.not.isopen) exit
    iout = iout + 1
    if (iout > 99) exit
  end do
  if (iout > 99) then 
    write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
    return
  end if
  open(iout,file=fname,iostat=info)
  if (info == psb_success_) then 
    call a%a%print(iout,iv,eirs,eics,head,ivr,ivc)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_n_sparse_print


subroutine psb_c_get_neigh(a,idx,neigh,n,info,lev)
  use psb_c_mat_mod, psb_protect_name => psb_c_get_neigh
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(in) :: a   
  integer, intent(in)                :: idx 
  integer, intent(out)               :: n   
  integer, allocatable, intent(out)  :: neigh(:)
  integer, intent(out)               :: info
  integer, optional, intent(in)      :: lev 

  Integer :: err_act
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_neigh(idx,neigh,n,info,lev)

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_get_neigh



subroutine psb_c_csall(nr,nc,a,info,nz) 
  use psb_c_mat_mod, psb_protect_name => psb_c_csall
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(out) :: a
  integer, intent(in)             :: nr,nc
  integer, intent(out)            :: info
  integer, intent(in), optional   :: nz

  Integer :: err_act 
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  info = psb_success_
  allocate(psb_c_coo_sparse_mat :: a%a, stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  call a%a%allocate(nr,nc,nz)
  call a%set_bld() 

  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_csall


subroutine  psb_c_reallocate_nz(nz,a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_reallocate_nz
  use psb_error_mod
  implicit none 
  integer, intent(in) :: nz
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%reallocate(nz)

  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_reallocate_nz


subroutine  psb_c_free(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_free
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a

  if (allocated(a%a)) then 
    call a%a%free()
    deallocate(a%a) 
  endif

end subroutine psb_c_free


subroutine  psb_c_trim(a) 
  use psb_c_mat_mod, psb_protect_name => psb_c_trim
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%trim()

  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_trim



subroutine psb_c_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_c_mat_mod, psb_protect_name => psb_c_csput
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)

  Integer :: err_act
  character(len=20)  :: name='csput'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.a%is_bld()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info,gtl) 
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_csput


subroutine psb_c_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetptn
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_csgetptn


subroutine psb_c_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetrow
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_csgetrow




subroutine psb_c_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetblk
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  class(psb_cspmat_type), intent(out) :: b
  integer, intent(in)                  :: imin,imax
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    

  if (info == psb_success_) call a%a%csget(imin,imax,acoo,info,&
       & jmin,jmax,iren,append,rscale,cscale)
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_csgetblk




subroutine psb_c_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csclip
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  class(psb_cspmat_type), intent(out) :: b
  integer,intent(out)                  :: info
  integer, intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == psb_success_) call a%a%csclip(acoo,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_csclip


subroutine psb_c_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_b_csclip
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  type(psb_c_coo_sparse_mat), intent(out) :: b
  integer,intent(out)                  :: info
  integer, intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csclip(b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_b_csclip




subroutine psb_c_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv
  implicit none 
  class(psb_cspmat_type), intent(in)    :: a
  class(psb_cspmat_type), intent(out)   :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_c_base_sparse_mat), intent(in), optional :: mold


  class(psb_c_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    call b%set_dupl(dupl)
  else if (a%is_bld()) then 
    ! Does this make sense at all?? Who knows..
    call b%set_dupl(psb_dupl_def_)
  end if

  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

#if defined(HAVE_MOLD)
    allocate(altmp, mold=mold,stat=info) 
#else
    call mold%mold(altmp,info)
#endif

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_c_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_c_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_c_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psb_c_csr_sparse_mat :: altmp, stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(psb_err_unit,*) 'Converting from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%cp_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,b%a)
  call b%set_asb() 
  call b%trim()
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_cscnv



subroutine psb_c_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv_ip
  implicit none 

  class(psb_cspmat_type), intent(inout) :: a
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psb_c_base_sparse_mat), intent(in), optional :: mold


  class(psb_c_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv_ip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    call a%set_dupl(dupl)
  else if (a%is_bld()) then 
    call a%set_dupl(psb_dupl_def_)
  end if

  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

#if defined(HAVE_MOLD)
    allocate(altmp, mold=mold,stat=info) 
#else
    call mold%mold(altmp,info)
#endif

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_c_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_c_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_c_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psb_c_csr_sparse_mat :: altmp, stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(psb_err_unit,*) 'Converting in-place from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%mv_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,a%a)
  call a%set_asb() 
  call a%trim()
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_cscnv_ip



subroutine psb_c_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv_base
  implicit none 
  class(psb_cspmat_type), intent(in)       :: a
  class(psb_c_base_sparse_mat), intent(out) :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl


  type(psb_c_coo_sparse_mat)  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cp_to_coo(altmp,info )
  if ((info == psb_success_).and.present(dupl)) then 
    call altmp%set_dupl(dupl)
  end if
  call altmp%fix(info)
  if (info == psb_success_) call altmp%trim()
  if (info == psb_success_) call altmp%set_asb() 
  if (info == psb_success_) call b%mv_from_coo(altmp,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_cscnv_base



subroutine psb_c_clip_d(a,b,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_clip_d
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  class(psb_cspmat_type), intent(out) :: b
  integer,intent(out)                  :: info

  Integer :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  integer :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == psb_success_) call a%a%cp_to_coo(acoo,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  nz = acoo%get_nzeros()
  j = 0
  do i=1, nz
    if (acoo%ia(i) /= acoo%ja(i)) then 
      j = j + 1 
      acoo%ia(j)  = acoo%ia(i)
      acoo%ja(j)  = acoo%ja(i)
      acoo%val(j) = acoo%val(i)
    end if
  end do
  call acoo%set_nzeros(j)
  call acoo%trim()
  call b%mv_from(acoo)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_clip_d



subroutine psb_c_clip_d_ip(a,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_clip_d_ip
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  integer,intent(out)                  :: info

  Integer :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  integer :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == psb_success_) call a%a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  nz = acoo%get_nzeros()
  j = 0
  do i=1, nz
    if (acoo%ia(i) /= acoo%ja(i)) then 
      j = j + 1 
      acoo%ia(j)  = acoo%ia(i)
      acoo%ja(j)  = acoo%ja(i)
      acoo%val(j) = acoo%val(i)
    end if
  end do
  call acoo%set_nzeros(j)
  call acoo%trim()
  call a%mv_from(acoo)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_clip_d_ip


subroutine psb_c_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_from
  implicit none 
  class(psb_cspmat_type), intent(out) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer :: info

#if defined(HAVE_MOLD)
  allocate(a%a,mold=b, stat=info)
#else
  call b%mold(a%a,info)
#endif
  call a%a%mv_from_fmt(b,info)
  call b%free()

  return
end subroutine psb_c_mv_from


subroutine psb_c_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_from
  implicit none 
  class(psb_cspmat_type), intent(out) :: a
  class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
  Integer :: err_act, info
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

#if defined(HAVE_MOLD)
  allocate(a%a,mold=b,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
#else
  call b%mold(a%a,info)
#endif
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call a%a%cp_from_fmt(b, info)    
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine psb_c_cp_from


subroutine psb_c_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_to
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(out) :: b
  integer :: info

  call b%mv_from_fmt(a%a,info)

  return
end subroutine psb_c_mv_to


subroutine psb_c_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_to
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  class(psb_c_base_sparse_mat), intent(out) :: b
  integer :: info

  call b%cp_from_fmt(a%a,info)

  return
end subroutine psb_c_cp_to



subroutine psb_cspmat_type_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_cspmat_type_move
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(out)   :: b
  integer, intent(out)                   :: info

  Integer :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = psb_success_
  call move_alloc(a%a,b%a)

  return
end subroutine psb_cspmat_type_move


subroutine psb_cspmat_type_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_cspmat_type_clone
  implicit none 
  class(psb_cspmat_type), intent(in)  :: a
  class(psb_cspmat_type), intent(out) :: b
  integer, intent(out)                 :: info

  Integer :: err_act
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

#if defined(HAVE_MOLD)
  allocate(b%a,mold=a%a,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
#else
  call a%a%mold(b%a,info)
#endif
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call b%a%cp_from_fmt(a%a, info)    
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_cspmat_type_clone



subroutine psb_c_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transp_1mat
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a

  Integer :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transp()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_transp_1mat



subroutine psb_c_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transp_2mat
  implicit none 
  class(psb_cspmat_type), intent(out) :: a
  class(psb_cspmat_type), intent(in)  :: b

  Integer :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (b%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

#if defined(HAVE_MOLD)
  allocate(a%a,mold=b%a,stat=info)
#else
  call b%a%mold(a%a,info)
#endif
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transp(b%a)    

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_transp_2mat


subroutine psb_c_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transc_1mat
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a

  Integer :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transc()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_transc_1mat



subroutine psb_c_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transc_2mat
  implicit none 
  class(psb_cspmat_type), intent(out) :: a
  class(psb_cspmat_type), intent(in)  :: b

  Integer :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (b%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

#if defined(HAVE_MOLD)
  allocate(a%a,mold=b%a,stat=info)
#else
  call b%a%mold(a%a,info)
#endif
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transc(b%a)    

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_transc_2mat




subroutine psb_c_reinit(a,clear)
  use psb_c_mat_mod, psb_protect_name => psb_c_reinit
  use psb_error_mod
  implicit none 

  class(psb_cspmat_type), intent(inout) :: a   
  logical, intent(in), optional :: clear
  Integer :: err_act, info
  character(len=20)  :: name='reinit'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%reinit(clear)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_c_reinit




! == ===================================
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
! == ===================================


subroutine psb_c_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csmm
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans
  Integer :: err_act
  character(len=20)  :: name='psb_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csmm(alpha,x,beta,y,info,trans) 
  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_csmm


subroutine psb_c_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csmv
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans
  Integer :: err_act
  character(len=20)  :: name='psb_csmv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csmm(alpha,x,beta,y,info,trans) 
  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_csmv


subroutine psb_c_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cssm
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_spk_), intent(in), optional :: d(:)
  Integer :: err_act
  character(len=20)  :: name='psb_cssm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cssm(alpha,x,beta,y,info,trans,scale,d) 
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_cssm


subroutine psb_c_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cssv
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_spk_), intent(in), optional :: d(:)
  Integer :: err_act
  character(len=20)  :: name='psb_cssv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cssm(alpha,x,beta,y,info,trans,scale,d) 

  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_cssv



function psb_c_csnmi(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_c_csnmi
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_)         :: res

  Integer :: err_act, info
  character(len=20)  :: name='csnmi'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%csnmi()
  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end function psb_c_csnmi


subroutine psb_c_get_diag(a,d,info)
  use psb_c_mat_mod, psb_protect_name => psb_c_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(out)          :: d(:)
  integer, intent(out)                 :: info

  Integer :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_diag(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_get_diag


subroutine psb_c_scal(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_scal
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d(:)
  integer, intent(out)                    :: info

  Integer :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_scal


subroutine psb_c_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_scals
  implicit none 
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d
  integer, intent(out)                    :: info

  Integer :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_c_scals



