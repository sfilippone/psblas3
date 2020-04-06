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
! File: psb_cgetmatinfo.f90
!
! This function containts the implementation for obtaining information on the
! paralle sparse matrix
!
function  psb_cget_nnz(a,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_cget_nnz
  use psi_mod
  use mpi

  implicit none

  integer(psb_lpk_)                     :: res
  type(psb_cspmat_type), intent(in)   :: a
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iia, jja
  integer(psb_lpk_) :: m,n,ia,ja,localnnz
  character(len=20) :: name, ch_err
  !
  name='psb_cget_nnz'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()

  ! Check for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  localnnz = a%get_nzeros()

  call MPI_ALLREDUCE(localnnz, res, 1, MPI_LONG, MPI_SUM, ictxt, info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end function

function  psb_c_is_matupd(a,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_c_is_matupd
  use psi_mod

  implicit none

  logical                               :: res
  type(psb_cspmat_type), intent(in)   :: a
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iia, jja
  integer(psb_lpk_) :: m,n,ia,ja,localnnz
  character(len=20) :: name, ch_err
  !
  name='psb_cis_matupd'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()

  ! Check for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  res = a%is_upd()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end function

function  psb_c_is_matasb(a,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_c_is_matasb
  use psi_mod

  implicit none

  logical                               :: res
  type(psb_cspmat_type), intent(in)   :: a
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iia, jja
  integer(psb_lpk_) :: m,n,ia,ja,localnnz
  character(len=20) :: name, ch_err
  !
  name='psb_cis_matasb'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()

  ! Check for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  res = a%is_asb()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end function

function  psb_c_is_matbld(a,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_c_is_matbld
  use psi_mod

  implicit none

  logical                               :: res
  type(psb_cspmat_type), intent(in)   :: a
  type(psb_desc_type), intent (in)      :: desc_a
  integer(psb_ipk_), intent(out)        :: info

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iia, jja
  integer(psb_lpk_) :: m,n,ia,ja,localnnz
  character(len=20) :: name, ch_err
  !
  name='psb_cis_matbld'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ictxt=desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()

  ! Check for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  res = a%is_bld()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end function
