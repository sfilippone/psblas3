subroutine psb_d_mat_renum(alg,mat,info,perm)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_mat_renum
  implicit none 
  integer, intent(in) :: alg
  type(psb_dspmat_type), intent(inout) :: mat
  integer, intent(out) :: info
  integer, allocatable, optional, intent(out) :: perm(:)
  
  integer            :: err_act, nr, nc
  character(len=20)  :: name

  info = psb_success_
  name = 'mat_renum'
  call psb_erractionsave(err_act)

  info = psb_success_
  
  nr = mat%get_nrows()
  nc = mat%get_ncols()
  if (nr /= nc) then 
    info = psb_err_rectangular_mat_unsupported_
    call psb_errpush(info,name,i_err=(/nr,nc,0,0,0/))
    goto 9999
  end if

  select case (alg)
  case(psb_mat_renum_gps_) 

    call psb_mat_renum_gps(mat,info,perm)

  case(psb_mat_renum_amd_) 

    call psb_mat_renum_amd(mat,info,perm)

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

  subroutine psb_mat_renum_gps(a,info,operm)
    use psb_base_mod
    use psb_gps_mod
    implicit none 
    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(out) :: info
    integer, allocatable, optional, intent(out) :: operm(:)

    ! 
    class(psb_d_base_sparse_mat), allocatable :: aa
    type(psb_d_csr_sparse_mat)  :: acsr
    type(psb_d_coo_sparse_mat)  :: acoo
    
    integer :: err_act
    character(len=20)           :: name
    integer, allocatable :: ndstk(:,:), iold(:), ndeg(:), perm(:) 
    integer :: i, j, k, ideg, nr, ibw, ipf, idpth

    info = psb_success_
    name = 'mat_renum_gps'
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
      ndstk(i,:) = 0
      k  = 0
      do j=acsr%irp(i),acsr%irp(i+1)-1
        k = k + 1
        ndstk(i,k) = acsr%ja(j)
      end do
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
    if (present(operm)) then 
      call psb_realloc(nr,operm,info)
      if (info /= psb_success_) then 
        info = psb_err_alloc_dealloc_
        call psb_errpush(info,name) 
        goto 9999
      end if
      operm(1:nr) = perm(1:nr)
    end if

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


  subroutine psb_mat_renum_amd(a,info,operm)
    use psb_base_mod    
    implicit none 
    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(out) :: info
    integer, allocatable, optional, intent(out) :: operm(:)

    ! 
    class(psb_d_base_sparse_mat), allocatable :: aa
    type(psb_d_csr_sparse_mat)  :: acsr
    type(psb_d_coo_sparse_mat)  :: acoo
    
    integer :: err_act
    character(len=20)           :: name
    integer, allocatable :: ndstk(:,:), iold(:), ndeg(:), perm(:) 
    integer :: i, j, k, ideg, nr, ibw, ipf, idpth

    info = psb_success_
    name = 'mat_renum_amd'
    call psb_erractionsave(err_act)

#if defined(HAVE_AMD) && defined(HAVE_ISO_C_BINDING)

    info = psb_success_
    nr = a%get_nrows()
    allocate(operm(nr))
    do i=1, nr
      operm(i) = i
    end do
!!$    call a%mold(aa)
!!$    call a%mv_to(aa)
!!$    call aa%mv_to_fmt(acsr,info)
!!$    ! Insert call to gps_reduce
!!$    nr   = acsr%get_nrows()
!!$    ideg = 0
!!$    do i=1, nr
!!$      ideg = max(ideg,acsr%irp(i+1)-acsr%irp(i))
!!$    end do
!!$    allocate(ndstk(nr,ideg), iold(nr), perm(nr+1), ndeg(nr),stat=info)
!!$    if (info /= 0) then
!!$      info = psb_err_alloc_dealloc_
!!$      call psb_errpush(info, name) 
!!$      goto 9999
!!$    end if
!!$    do i=1, nr
!!$      iold(i) = i 
!!$      ndstk(i,:) = 0
!!$      k  = 0
!!$      do j=acsr%irp(i),acsr%irp(i+1)-1
!!$        k = k + 1
!!$        ndstk(i,k) = acsr%ja(j)
!!$      end do
!!$    end do
!!$    perm = 0
!!$
!!$    call psb_gps_reduce(ndstk,nr,ideg,iold,perm,ndeg,ibw,ipf,idpth)
!!$
!!$    if (.not.psb_isaperm(nr,perm)) then 
!!$      write(0,*) 'Something wrong: bad perm from gps_reduce'
!!$      info = psb_err_from_subroutine_
!!$      call psb_errpush(info,name) 
!!$      goto 9999
!!$    end if
!!$    ! Move to coordinate to apply renumbering
!!$    call acsr%mv_to_coo(acoo,info)
!!$    do i=1, acoo%get_nzeros()
!!$      acoo%ia(i) = perm(acoo%ia(i))
!!$      acoo%ja(i) = perm(acoo%ja(i))
!!$    end do
!!$    call acoo%fix(info) 
!!$
!!$    ! Get back to where we started from
!!$    call aa%mv_from_coo(acoo,info)
!!$    call a%mv_from(aa)
!!$    if (present(operm)) then 
!!$      call psb_realloc(nr,operm,info)
!!$      if (info /= psb_success_) then 
!!$        info = psb_err_alloc_dealloc_
!!$        call psb_errpush(info,name) 
!!$        goto 9999
!!$      end if
!!$      operm(1:nr) = perm(1:nr)
!!$    end if
!!$
!!$    deallocate(aa)
#else 
    info = psb_err_missing_aux_lib_
    call psb_errpush(info, name) 
    goto 9999
#endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_mat_renum_amd

end subroutine psb_d_mat_renum


subroutine psb_d_cmp_bwpf(mat,bwl,bwu,prf,info)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_cmp_bwpf
  implicit none 
  type(psb_dspmat_type), intent(in) :: mat
  integer, intent(out) :: bwl, bwu
  integer, intent(out) :: prf
  integer, intent(out) :: info
  !
  integer, allocatable :: irow(:), icol(:)
  real(psb_dpk_), allocatable :: val(:)
  integer :: nz
  integer :: i, j, lrbu, lrbl
  
  info = psb_success_
  bwl = 0
  bwu = 0
  prf = 0
  select type (aa=>mat%a)
  class is (psb_d_csr_sparse_mat)
    do i=1, aa%get_nrows()
      lrbl = 0
      lrbu = 0
      do j = aa%irp(i), aa%irp(i+1) - 1
        lrbl = max(lrbl,i-aa%ja(j))
        lrbu = max(lrbu,aa%ja(j)-i)
      end do
      prf = prf + lrbl+lrbu
      bwu  = max(bwu,lrbu)
      bwl  = max(bwl,lrbu)
    end do
    
  class default
    do i=1, aa%get_nrows()
      lrbl = 0
      lrbu = 0
      call aa%csget(i,i,nz,irow,icol,val,info)
      if (info /= psb_success_) return
      do j=1, nz
        lrbl = max(lrbl,i-icol(j))
        lrbu = max(lrbu,icol(j)-i)
      end do
      prf = prf + lrbl+lrbu
      bwu  = max(bwu,lrbu)
      bwl  = max(bwl,lrbu)
    end do
  end select
  
end subroutine psb_d_cmp_bwpf
