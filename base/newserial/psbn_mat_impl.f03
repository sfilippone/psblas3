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
  call a%set_state(psbn_spmat_bld_) 
  return

end subroutine psbn_d_csall


subroutine psbn_d_csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
  use psbn_d_base_mat_mod
  use psb_error_mod
  use psbn_d_mat_mod, psb_protect_name => psbn_d_csins
  implicit none 
  type(psbn_d_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)

  Integer :: err_act
  character(len=20)  :: name='psbn_csins'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (.not.a%is_bld()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csins(nz,val,ia,ja,imin,imax,jmin,jmax,info,gtl) 
  if (info /= 0) goto 9999 

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  
end subroutine psbn_d_csins



subroutine psbn_d_spcnv(a,b,info,type,mold,upd,dupl)
  use psbn_d_mat_mod, psb_protect_name => psbn_d_spcnv
  use psb_realloc_mod
  use psb_sort_mod
  implicit none 
  type(psbn_d_sparse_mat), intent(in)    :: a
  type(psbn_d_sparse_mat), intent(out)   :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold
  

  write(0,*) 'TO BE IMPLEMENTED '

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
  class(psbn_d_base_sparse_mat), pointer :: aslct
  type(psbn_d_csr_sparse_mat)    :: csrtmp
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
  else
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
    allocate(altmp, source=csrtmp,stat=info) 
  end if

  select type ( aa => a%a ) 
  class is (psbn_d_coo_sparse_mat) 
    ! Quick route from coo 
    call altmp%mv_from_coo(aa, info)
  class default
    call altmp%mv_from_fmt(aa, info)
  end select
  
  if (info /= 0) then
    info = 1121
    call psb_errpush(info,name)
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
