subroutine psbn_d_csall(nr,nc,a,info,nz) 
  use psbn_d_base_mat_mod
  use psb_realloc_mod
  use psb_sort_mod
  use psbn_d_mat_mod, psb_protect_name => psbn_d_csall
  implicit none 
  type(psbn_d_sparse_mat), intent(out) :: a
  integer, intent(in)             :: nr,nc
  integer, intent(out)            :: info
  integer, intent(in), optional   :: nz
  

  info = 0
  call a%allocate(nr,nc,nz)
  call a%set_bld() 
  return

end subroutine psbn_d_csall


subroutine psbn_d_csput(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
  use psbn_d_base_mat_mod
  use psb_error_mod
  use psbn_d_mat_mod, psb_protect_name => psbn_d_csput
  implicit none 
  type(psbn_d_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)

  Integer :: err_act
  character(len=20)  :: name='psbn_csput'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.a%is_bld()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csput(nz,val,ia,ja,imin,imax,jmin,jmax,info,gtl) 
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  
end subroutine psbn_d_csput



subroutine psbn_d_spcnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psbn_d_mat_mod, psb_protect_name => psbn_d_spcnv
  implicit none 
  type(psbn_d_sparse_mat), intent(in)    :: a
  type(psbn_d_sparse_mat), intent(out)   :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold
  

  class(psbn_d_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='psbn_cscnv'
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
    call b%set_dupl(psbn_dupl_def_)
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
      allocate(psbn_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psbn_d_coo_sparse_mat :: altmp, stat=info) 
    case default
      info = 136 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psbn_d_csr_sparse_mat :: altmp, stat=info) 
  end if
  
  if (info /= 0) then 
    info = 4000
    call psb_errpush(info,name)
    goto 9999
  end if

  call altmp%cp_from_fmt(a%a, info)
  
  if (info /= 0) then
    info = 4010
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,b%a)
  call b%set_asb() 

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psbn_d_spcnv

subroutine psbn_d_spcnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psbn_d_mat_mod, psb_protect_name => psbn_d_spcnv_ip
  implicit none 

  type(psbn_d_sparse_mat), intent(inout) :: a
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold


  class(psbn_d_base_sparse_mat), allocatable  :: altmp
  Integer :: err_act
  character(len=20)  :: name='psbn_cscnv'
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
    call a%set_dupl(psbn_dupl_def_)
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
      allocate(psbn_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psbn_d_coo_sparse_mat :: altmp, stat=info) 
    case default
      info = 136 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(psbn_d_csr_sparse_mat :: altmp, stat=info) 
  end if
  
  if (info /= 0) then 
    info = 4000
    call psb_errpush(info,name)
    goto 9999
  end if

  call altmp%mv_from_fmt(a%a, info)
  
  if (info /= 0) then
    info = 4010
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,a%a)
  call a%set_asb() 

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

end subroutine psbn_d_spcnv_ip
