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
!         software without specific prior written permission.
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
! File: psb_sglobtril.f90
!
! Subroutine: psb_sglobtril
!
Subroutine psb_sglobtril(a,desc_a,b,info,diag,imin,imax,jmin,jmax)

  use psb_base_mod, psb_protect_name => psb_sglobtril
  use psi_mod

  Implicit None

  !     .. Array Arguments ..
  Type(psb_sspmat_type), Intent(inout)       ::  a
  Type(psb_desc_type), Intent(inout), target :: desc_a
  Type(psb_sspmat_type), Intent(out)       ::  b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax

  !     .. Local Scalars ..
  integer(psb_ipk_) ::  i, j, err_act,m,&
       &  nz
  integer(psb_lpk_) :: gidx, lnz, gnr, gnc,lnr,lnc
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: me, np
  integer(psb_mpk_) :: icomm, minfo

  type(psb_s_coo_sparse_mat)  :: dtcoo
  type(psb_ls_coo_sparse_mat) :: ldtcoo, ldcootril
  type(psb_lsspmat_type)     :: ldglob,ldtril
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err
  character(len=50) :: fname
  integer :: iout

  name='psb_sglobtril'
  info  = psb_success_
  if (psb_errstatus_fatal()) return
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if
  ctxt = desc_a%get_context()
  icomm = ctxt%get_mpic()
  Call psb_info(ctxt, me, np)

  If (debug_level >= psb_debug_outer_) &
       & Write(debug_unit,*) me,' ',trim(name),&
       & ': start',diag

  gnr = desc_a%get_global_rows()
  gnc = desc_a%get_global_cols()
  call a%a%cp_to_lcoo(ldtcoo,info)
  lnz = ldtcoo%get_nzeros()
  call desc_a%l2gip(ldtcoo%ia(1:lnz),info,owned=.false.)
  call desc_a%l2gip(ldtcoo%ja(1:lnz),info,owned=.false.)
  call ldtcoo%tril(ldcootril,info,&
       & diag=diag,imax=gnr,jmax=gnc)
  lnz = ldcootril%get_nzeros()

  lnr = desc_a%get_local_rows()
  lnc = desc_a%get_local_cols()
  call desc_a%g2lip(ldcootril%ia(1:lnz),info,owned=.false.)
  call desc_a%g2lip(ldcootril%ja(1:lnz),info,owned=.false.)
  call ldcootril%set_nrows(lnr)
  call ldcootril%set_ncols(lnc)
  call ldcootril%fix(info)
  call b%mv_from_lb(ldcootril)
  call b%cscnv(info,mold=a%a)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': end'
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

End Subroutine psb_sglobtril

