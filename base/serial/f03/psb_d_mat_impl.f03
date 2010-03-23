!=====================================
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
!=====================================


subroutine  psb_d_set_nrows(m,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_nrows
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  integer, intent(in) :: m
  Integer :: err_act, info
  character(len=20)  :: name='set_nrows'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_nrows


subroutine  psb_d_set_ncols(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_ncols
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_ncols



subroutine  psb_d_set_state(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_state
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_state



subroutine  psb_d_set_dupl(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_dupl
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  integer, intent(in) :: n
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_dupl


subroutine  psb_d_set_null(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_null
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_null


subroutine  psb_d_set_bld(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_bld
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_bld


subroutine  psb_d_set_upd(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_upd
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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


end subroutine psb_d_set_upd


subroutine  psb_d_set_asb(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_asb
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_asb


subroutine psb_d_set_sorted(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_sorted
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_sorted


subroutine psb_d_set_triangle(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_triangle
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_triangle


subroutine psb_d_set_unit(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_unit
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_unit


subroutine psb_d_set_lower(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_lower
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_lower


subroutine psb_d_set_upper(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_upper
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  logical, intent(in), optional :: val
  Integer :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_set_upper



!=====================================
!
!
!
! Data management
!
!
!
!
!
!=====================================  


subroutine psb_d_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_d_mat_mod, psb_protect_name => psb_d_sparse_print
  use psb_error_mod
  implicit none 

  integer, intent(in)               :: iout
  class(psb_d_sparse_mat), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = 0
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_sparse_print




subroutine psb_d_get_neigh(a,idx,neigh,n,info,lev)
  use psb_d_mat_mod, psb_protect_name => psb_d_get_neigh
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a   
  integer, intent(in)                :: idx 
  integer, intent(out)               :: n   
  integer, allocatable, intent(out)  :: neigh(:)
  integer, intent(out)               :: info
  integer, optional, intent(in)      :: lev 

  Integer :: err_act
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_neigh(idx,neigh,n,info,lev)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_get_neigh



subroutine psb_d_csall(nr,nc,a,info,nz) 
  use psb_d_mat_mod, psb_protect_name => psb_d_csall
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(out) :: a
  integer, intent(in)             :: nr,nc
  integer, intent(out)            :: info
  integer, intent(in), optional   :: nz

  Integer :: err_act 
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  info = 0
  allocate(psb_d_coo_sparse_mat :: a%a, stat=info)
  if (info /= 0) then 
    info = 4000 
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

end subroutine psb_d_csall


subroutine  psb_d_reallocate_nz(nz,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_reallocate_nz
  use psb_error_mod
  implicit none 
  integer, intent(in) :: nz
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_reallocate_nz


subroutine  psb_d_free(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_free
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='free'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%free()
  deallocate(a%a) 
  return

9999 continue

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_free


subroutine  psb_d_trim(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_trim
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end subroutine psb_d_trim



subroutine psb_d_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_d_mat_mod, psb_protect_name => psb_d_csput
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)

  Integer :: err_act
  character(len=20)  :: name='csput'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.a%is_bld()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info,gtl) 
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_csput


subroutine psb_d_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetptn
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
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

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_csgetptn


subroutine psb_d_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetrow
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_csgetrow




subroutine psb_d_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetblk
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
  class(psb_d_sparse_mat), intent(out) :: b
  integer, intent(in)                  :: imin,imax
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo


  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    

  if (info == 0) call a%a%csget(imin,imax,acoo,info,&
       & jmin,jmax,iren,append,rscale,cscale)
  if (info == 0) call move_alloc(acoo,b%a)
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_csgetblk




subroutine psb_d_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csclip
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
  class(psb_d_sparse_mat), intent(out) :: b
  integer,intent(out)                  :: info
  integer, intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == 0) call a%a%csclip(acoo,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info == 0) call move_alloc(acoo,b%a)
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_csclip


subroutine psb_d_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_b_csclip
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
  type(psb_d_coo_sparse_mat), intent(out) :: b
  integer,intent(out)                  :: info
  integer, intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csclip(b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_b_csclip




subroutine psb_d_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv
  implicit none 
  class(psb_d_sparse_mat), intent(in)    :: a
  class(psb_d_sparse_mat), intent(out)   :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_d_base_sparse_mat), intent(in), optional :: mold


  class(psb_d_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = 1121
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
    info = 583
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, source=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_d_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_d_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = 136 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
  end if

  if (info /= 0) then 
    info = 4000
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(0,*) 'Converting from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%cp_from_fmt(a%a, info)

  if (info /= 0) then
    info = 4010
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

end subroutine psb_d_cscnv



subroutine psb_d_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv_ip
  implicit none 

  class(psb_d_sparse_mat), intent(inout) :: a
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psb_d_base_sparse_mat), intent(in), optional :: mold


  class(psb_d_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv_ip'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    call a%set_dupl(dupl)
  else if (a%is_bld()) then 
    call a%set_dupl(psb_dupl_def_)
  end if

  if (count( (/present(mold),present(type) /)) > 1) then
    info = 583
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, source=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_d_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_d_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = 136 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
  end if

  if (info /= 0) then 
    info = 4000
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(0,*) 'Converting in-place from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%mv_from_fmt(a%a, info)

  if (info /= 0) then
    info = 4010
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

end subroutine psb_d_cscnv_ip



subroutine psb_d_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv_base
  implicit none 
  class(psb_d_sparse_mat), intent(in)       :: a
  class(psb_d_base_sparse_mat), intent(out) :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl


  type(psb_d_coo_sparse_mat)  :: altmp
  Integer :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cp_to_coo(altmp,info )
  if ((info == 0).and.present(dupl)) then 
    call altmp%set_dupl(dupl)
  end if
  call altmp%fix(info)
  if (info == 0) call altmp%trim()
  if (info == 0) call altmp%set_asb() 
  if (info == 0) call b%mv_from_coo(altmp,info)

  if (info /= 0) then
    info = 4010
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

end subroutine psb_d_cscnv_base



subroutine psb_d_clip_d(a,b,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_clip_d
  implicit none

  class(psb_d_sparse_mat), intent(in) :: a
  class(psb_d_sparse_mat), intent(out) :: b
  integer,intent(out)                  :: info

  Integer :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo
  integer :: i, j, nz

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == 0) call a%a%cp_to_coo(acoo,info)
  if (info /= 0) then 
    info = 4000
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

end subroutine psb_d_clip_d



subroutine psb_d_clip_d_ip(a,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_clip_d_ip
  implicit none

  class(psb_d_sparse_mat), intent(inout) :: a
  integer,intent(out)                  :: info

  Integer :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo
  integer :: i, j, nz

  info = 0
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == 0) call a%a%mv_to_coo(acoo,info)
  if (info /= 0) then
    info = 4000
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

end subroutine psb_d_clip_d_ip


subroutine psb_d_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_from
  implicit none 
  class(psb_d_sparse_mat), intent(out) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer :: info

  allocate(a%a,source=b, stat=info)
  call a%a%mv_from_fmt(b,info)

  return
end subroutine psb_d_mv_from


subroutine psb_d_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_from
  implicit none 
  class(psb_d_sparse_mat), intent(out) :: a
  class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
  Integer :: err_act, info
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0

  allocate(a%a,source=b,stat=info)
  if (info /= 0) info = 4000
  if (info == 0) call a%a%cp_from_fmt(b, info)    
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
end subroutine psb_d_cp_from


subroutine psb_d_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_to
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(out) :: b
  integer :: info

  call b%mv_from_fmt(a%a,info)

  return
end subroutine psb_d_mv_to


subroutine psb_d_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_to
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  class(psb_d_base_sparse_mat), intent(out) :: b
  integer :: info

  call b%cp_from_fmt(a%a,info)

  return
end subroutine psb_d_cp_to



subroutine psb_d_sparse_mat_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_sparse_mat_move
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  class(psb_d_sparse_mat), intent(out)   :: b
  integer, intent(out)                   :: info

  Integer :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = 0
  call move_alloc(a%a,b%a)

  return
end subroutine psb_d_sparse_mat_move


subroutine psb_d_sparse_mat_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_sparse_mat_clone
  implicit none 
  class(psb_d_sparse_mat), intent(in)  :: a
  class(psb_d_sparse_mat), intent(out) :: b
  integer, intent(out)                 :: info

  Integer :: err_act
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0

  allocate(b%a,source=a%a,stat=info)
  if (info /= 0) info = 4000
  if (info == 0) call b%a%cp_from_fmt(a%a, info)    
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psb_d_sparse_mat_clone



subroutine psb_d_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transp_1mat
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a

  Integer :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
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

end subroutine psb_d_transp_1mat



subroutine psb_d_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transp_2mat
  implicit none 
  class(psb_d_sparse_mat), intent(out) :: a
  class(psb_d_sparse_mat), intent(in)  :: b

  Integer :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (b%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(a%a,source=b%a,stat=info)
  if (info /= 0) then 
    info = 4000
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

end subroutine psb_d_transp_2mat


subroutine psb_d_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transc_1mat
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a

  Integer :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
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

end subroutine psb_d_transc_1mat



subroutine psb_d_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transc_2mat
  implicit none 
  class(psb_d_sparse_mat), intent(out) :: a
  class(psb_d_sparse_mat), intent(in)  :: b

  Integer :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (b%is_null()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(a%a,source=b%a,stat=info)
  if (info /= 0) then 
    info = 4000
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

end subroutine psb_d_transc_2mat




subroutine psb_d_reinit(a,clear)
  use psb_d_mat_mod, psb_protect_name => psb_d_reinit
  use psb_error_mod
  implicit none 

  class(psb_d_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear
  Integer :: err_act, info
  character(len=20)  :: name='reinit'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = 1121
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

end subroutine psb_d_reinit




!=====================================
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
!=====================================


subroutine psb_d_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csmm
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans
  Integer :: err_act
  character(len=20)  :: name='psb_csmm'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csmm(alpha,x,beta,y,info,trans) 
  if (info /= 0) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csmm


subroutine psb_d_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csmv
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans
  Integer :: err_act
  character(len=20)  :: name='psb_csmv'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csmm(alpha,x,beta,y,info,trans) 
  if (info /= 0) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csmv


subroutine psb_d_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cssm
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)
  Integer :: err_act
  character(len=20)  :: name='psb_cssm'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cssm(alpha,x,beta,y,info,trans,scale,d) 
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_cssm


subroutine psb_d_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cssv
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)
  Integer :: err_act
  character(len=20)  :: name='psb_cssv'
  logical, parameter :: debug=.false.

  info = 0 
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cssm(alpha,x,beta,y,info,trans,scale,d) 

  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_cssv



function psb_d_csnmi(a) result(res)
  use psb_d_mat_mod, psb_protect_name => psb_d_csnmi
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  Integer :: err_act, info
  character(len=20)  :: name='csnmi'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
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

end function psb_d_csnmi


subroutine psb_d_get_diag(a,d,info)
  use psb_d_mat_mod, psb_protect_name => psb_d_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_d_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)          :: d(:)
  integer, intent(out)                 :: info

  Integer :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_diag(d,info)
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_get_diag


subroutine psb_d_scal(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_scal
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d(:)
  integer, intent(out)                    :: info

  Integer :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_scal


subroutine psb_d_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_scals
  implicit none 
  class(psb_d_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d
  integer, intent(out)                    :: info

  Integer :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_scals



