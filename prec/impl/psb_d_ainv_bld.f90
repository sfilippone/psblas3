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
!    Moved here from AMG-AINV, original copyright below.
!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
!
subroutine psb_d_ainv_bld(a,alg,fillin,thresh,wmat,d,zmat,desc,info,blck,iscale)

  use psb_base_mod
  use psb_prec_const_mod
  use psb_d_biconjg_mod

  implicit none

  ! Arguments
  type(psb_dspmat_type), intent(in), target   :: a
  integer(psb_ipk_), intent(in)                 :: fillin,alg
  real(psb_dpk_), intent(in)                  :: thresh
  type(psb_dspmat_type), intent(inout)        :: wmat, zmat
  real(psb_dpk_), allocatable                  :: d(:)
  Type(psb_desc_type), Intent(in)               :: desc
  integer(psb_ipk_), intent(out)                :: info
  type(psb_dspmat_type), intent(in), optional :: blck
  integer(psb_ipk_), intent(in), optional       :: iscale
  !
  integer(psb_ipk_)             :: i, nztota, err_act, n_row, nrow_a
  type(psb_d_coo_sparse_mat)  :: acoo
  type(psb_d_csr_sparse_mat)  :: acsr
  type(psb_dspmat_type)       :: atmp
  real(psb_dpk_), allocatable :: arws(:), acls(:)
  real(psb_dpk_), allocatable  :: pq(:), ad(:)
  integer(psb_ipk_)             :: debug_level, debug_unit
  type(psb_ctxt_type)           :: ctxt
  integer(psb_ipk_)             :: np,me
  integer(psb_ipk_)             :: nzrmax, iscale_
  real(psb_dpk_)              :: sp_thresh
  real(psb_dpk_)               :: weight
  character(len=20)             :: name, ch_err


  info = psb_success_
  name = 'psb_dainv_bld'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt        = psb_cd_get_context(desc)
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  iscale_ = psb_ilu_scale_none_
  if (present(iscale)) iscale_ = iscale
  weight = done
  !
  ! Check the memory available to hold the W and Z factors
  ! and allocate it if needed
  !
  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if (present(blck)) then
    nztota = nztota + blck%get_nzeros()
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ': out get_nnzeros',nrow_a,nztota,&
       & a%get_nrows(),a%get_ncols(),a%get_nzeros()


  n_row  = desc%get_local_rows()
  allocate(pq(n_row),stat=info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  end if

  nzrmax    = fillin
  sp_thresh = thresh

  !
  ! Ok, let's start first with Z (i.e. Upper)
  !
  call a%csclip(acoo,info,imax=n_row,jmax=n_row)
  call acsr%mv_from_coo(acoo,info)
  select case(iscale_)
  case(psb_ilu_scale_none_)
    ! Ok, do nothing.

  case(psb_ilu_scale_maxval_)
    weight = acsr%maxval()
    weight = done/weight
    !call acsr%scal(weight,info)
    call psb_d_csr_scals(weight,a,info)

  case(psb_ilu_scale_arcsum_)
    allocate(arws(n_row),acls(n_row),ad(n_row),stat=info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    end if
    !call acsr%arwsum(arws)
    !call acsr%aclsum(acls)
    call psb_d_csr_arwsum(arws,acsr)
    call psb_d_csr_aclsum(acls,acsr)
    ad(1:n_row) = sqrt(sqrt(arws(1:n_row)*acls(1:n_row)))
    ad(1:n_row) = done/ad(1:n_row)
    !call acsr%scal(ad,info,side='L')
    !call acsr%scal(ad,info,side='R')
    call psb_d_csr_scal(ad,acsr,info,'L')
    call psb_d_csr_scal(ad,acsr,info,'R')
  case default
    call psb_errpush(psb_err_from_subroutine_,name,a_err='wrong iscale')
    goto 9999
  end select

  !
  ! Here for the actual workhorses.
  ! Only biconjg is surviving for now....
  !
    call psb_sparse_biconjg(alg,n_row,acsr,pq,&
         &   zmat,wmat,nzrmax,sp_thresh,info)
  ! Done. Hopefully....

  if (info /= psb_success_) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='sparse_biconjg')
    goto 9999
  end if
  call atmp%mv_from(acsr)

  !
  ! Is this right???
  !
  do i=1, n_row
    if (abs(pq(i)) < d_epstol) then
      pq(i) = done
    else
      pq(i) = done/pq(i)
    end if
  end do

  select case(iscale_)
  case(psb_ilu_scale_none_)
    ! Ok, do nothing.
  case(psb_ilu_scale_maxval_)
    pq(:) = pq(:)*weight

  case(psb_ilu_scale_arcsum_)
    call zmat%scal(ad,info,side='L')
    call wmat%scal(ad,info,side='R')

  case default
    call psb_errpush(psb_err_from_subroutine_,name,a_err='wrong iscale')
    goto 9999
  end select

  call psb_move_alloc(pq,d,info)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_d_ainv_bld
