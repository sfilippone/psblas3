subroutine psb_d_mat_renum(alg,mat,info)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_mat_renum
  implicit none 
  integer, intent(in) :: alg
  type(psb_dspmat_type), intent(inout) :: mat
  integer, intent(out) :: info
  
  integer :: err_act
  character(len=20)           :: name

  info = psb_success_
  name = 'mat_renum'
  call psb_erractionsave(err_act)

  info = psb_success_

  select case (alg)
  case(psb_renum_gps_) 

    call psb_mat_renum_gps(mat,info)

  case default
    info = psb_err_input_value_invalid_i_
    call psb_errpush(info,name,i_err=(/1,alg,0,0,0/))
    goto 9999 
  end select
  
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_non_
    call psb_errpush(info,name)
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
  return

contains

  subroutine psb_mat_renum_gps(a,info)
    use psb_base_mod
    use psb_gps_mod
    implicit none 
    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(out) :: info

    ! 
    class(psb_d_base_sparse_mat), allocatable :: aa
    type(psb_d_csr_sparse_mat)  :: acsr
    type(psb_d_coo_sparse_mat)  :: acoo
    
    integer :: err_act
    character(len=20)           :: name
    integer, allocatable :: ndstk(:,:), iold(:), ndeg(:), perm(:) 
    integer :: i, j, k, ideg, nr, ibw, ipf, idpth

    info = psb_success_
    name = 'mat_renum'
    call psb_erractionsave(err_act)

    info = psb_success_

    call a%mold(aa)
    call a%mv_to(aa)
    call aa%mv_to_fmt(acsr,info)
    ! Insert call to gps_reduce
    nr   = acsr%get_nrows()
    ideg = 0
    do i=1, nr
      ideg = max(ideg,acsr%irp(i+1)-acsr%irp(i))
    end do
    allocate(ndstk(nr,ideg), iold(nr), perm(nr+1), ndeg(nr),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, name) 
      goto 9999
    end if
    do i=1, nr
      iold(i) = i 
    end do
    perm = 0

    call psb_gps_reduce(ndstk,nr,ideg,iold,perm,ndeg,ibw,ipf,idpth)
    
    if (.not.psb_isaperm(nr,perm)) then 
      write(0,*) 'Something wrong: bad perm from gps_reduce'
      info = psb_err_from_subroutine_
      call psb_errpush(info,name) 
      goto 9999
    end if
    ! Move to coordinate to apply renumbering
    call acsr%mv_to_coo(acoo,info)
    do i=1, acoo%get_nzeros()
      acoo%ia(i) = perm(acoo%ia(i))
      acoo%ja(i) = perm(acoo%ja(i))
    end do
    call acoo%fix(info) 

    ! Get back to where we started from
    call aa%mv_from_coo(acoo,info)

    call a%mv_from(aa)
    
    deallocate(aa)
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_mat_renum_gps

end subroutine psb_d_mat_renum
