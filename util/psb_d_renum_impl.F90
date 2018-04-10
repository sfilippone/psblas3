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
subroutine psb_d_mat_renums(alg,mat,info,perm)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_mat_renums
  implicit none 
  character(len=*), intent(in) :: alg
  type(psb_dspmat_type), intent(inout) :: mat
  integer(psb_ipk_), intent(out) :: info
  integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
  
  integer(psb_ipk_) :: err_act, nr, nc, ialg
  character(len=20)  :: name

  info = psb_success_
  name = 'mat_renum'
  call psb_erractionsave(err_act)

  info = psb_success_
  select case (psb_toupper(alg))
  case ('GPS')
    ialg = psb_mat_renum_gps_
  case ('AMD')
    ialg = psb_mat_renum_amd_
  case ('NONE', 'ID') 
    ialg = psb_mat_renum_identity_
  case default
    write(0,*) 'Unknown algorithm "',psb_toupper(alg),'"'
    ialg = -1
  end select

  call psb_mat_renum(ialg,mat,info,perm)

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999 
  end if
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_d_mat_renums
  
subroutine psb_d_mat_renum(alg,mat,info,perm)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_mat_renum
  implicit none 
  integer(psb_ipk_), intent(in) :: alg
  type(psb_dspmat_type), intent(inout) :: mat
  integer(psb_ipk_), intent(out) :: info
  integer(psb_ipk_), allocatable, optional, intent(out) :: perm(:)
  
  integer(psb_ipk_) :: err_act, nr, nc, i, ierr(5)
  character(len=20)  :: name

  info = psb_success_
  name = 'mat_renum'
  call psb_erractionsave(err_act)

  info = psb_success_
  
  nr = mat%get_nrows()
  nc = mat%get_ncols()
  if (nr /= nc) then 
    info = psb_err_rectangular_mat_unsupported_
    ierr(1) = nr; ierr(2) = nc;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  select case (alg)
  case(psb_mat_renum_gps_) 

    call psb_mat_renum_gps(mat,info,perm)

  case(psb_mat_renum_amd_) 

    call psb_mat_renum_amd(mat,info,perm)

  case(psb_mat_renum_identity_) 
    nr = mat%get_nrows()
    allocate(perm(nr),stat=info) 
    if (info == 0) then 
      do i=1,nr
        perm(i) = i
      end do
    else
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif
  case default
    info = psb_err_input_value_invalid_i_
    ierr(1) = 1; ierr(2) = alg;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999 
  end select
  
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999 
  end if
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_mat_renum_gps(a,info,operm)
    use psb_base_mod
    use psb_gps_mod
    implicit none 
    type(psb_dspmat_type), intent(inout) :: a
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), allocatable, optional, intent(out) :: operm(:)

    ! 
    class(psb_d_base_sparse_mat), allocatable :: aa
    type(psb_d_csr_sparse_mat)  :: acsr
    type(psb_d_coo_sparse_mat)  :: acoo
    
    integer(psb_ipk_) :: err_act
    character(len=20)           :: name
    integer(psb_ipk_), allocatable :: ndstk(:,:), iold(:), ndeg(:), perm(:) 
    integer(psb_ipk_) :: i, j, k, ideg, nr, ibw, ipf, idpth

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

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_mat_renum_gps


  subroutine psb_mat_renum_amd(a,info,operm)
#if defined(HAVE_AMD) 
    use iso_c_binding
#endif
    use psb_base_mod    
    implicit none 
    type(psb_dspmat_type), intent(inout) :: a
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), allocatable, optional, intent(out) :: operm(:)

    ! 
#if defined(HAVE_AMD) 
    interface 
      function psb_amd_order(n,ap,ai,p)&
           & result(res) bind(c,name='psb_amd_order')
        use iso_c_binding
        integer(c_int) :: res
        integer(c_int), value :: n
        integer(c_int) :: ap(*), ai(*), p(*)
      end function psb_amd_order
    end interface
#endif

    type(psb_d_csc_sparse_mat)  :: acsc
    class(psb_d_base_sparse_mat), allocatable :: aa
    type(psb_d_coo_sparse_mat)  :: acoo

    integer(psb_ipk_), allocatable :: perm(:) 
    integer(psb_ipk_) :: err_act
    character(len=20)           :: name
    integer(psb_ipk_) :: i, j, k, ideg, nr, ibw, ipf, idpth, nz

    info = psb_success_
    name = 'mat_renum_amd'
    call psb_erractionsave(err_act)

#if defined(HAVE_AMD) 

    info = psb_success_
    nr   = a%get_nrows()
    nz   = a%get_nzeros()
    allocate(perm(nr),stat=info)
    if (info /= 0) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name) 
      goto 9999
    end if


    allocate(aa, mold=a%a)
    call a%mv_to(acsc)

    acsc%ia(:)  = acsc%ia(:) - 1
    acsc%icp(:) = acsc%icp(:) - 1

    info = psb_amd_order(nr,acsc%icp,acsc%ia,perm)
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_amd_order') 
      goto 9999
    end if

    perm(:)     = perm(:) + 1
    acsc%ia(:)  = acsc%ia(:) + 1
    acsc%icp(:) = acsc%icp(:) + 1

    call acsc%mv_to_coo(acoo,info)
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

    deallocate(aa,perm)

#else 

    info = psb_err_missing_aux_lib_
    call psb_errpush(info, name) 
    goto 9999
#endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return
  end subroutine psb_mat_renum_amd

end subroutine psb_d_mat_renum


subroutine psb_d_cmp_bwpf(mat,bwl,bwu,prf,info)
  use psb_base_mod
  use psb_renum_mod, psb_protect_name => psb_d_cmp_bwpf
  implicit none 
  type(psb_dspmat_type), intent(in) :: mat
  integer(psb_ipk_), intent(out) :: bwl, bwu
  integer(psb_ipk_), intent(out) :: prf
  integer(psb_ipk_), intent(out) :: info
  !
  integer(psb_ipk_), allocatable :: irow(:), icol(:)
  real(psb_dpk_), allocatable :: val(:)
  integer(psb_ipk_) :: nz
  integer(psb_ipk_) :: i, j, lrbu, lrbl
  
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
