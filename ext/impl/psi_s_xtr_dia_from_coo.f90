!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  
subroutine psi_s_xtr_dia_from_coo(nr,nc,nz,ia,ja,val,d,nrd,ncd,data,info,&
     & initdata, rdisp)    
  use psb_base_mod, only : psb_ipk_, psb_success_, psb_spk_, szero

  implicit none 
  integer(psb_ipk_), intent(in)  :: nr, nc, nz, nrd,ncd,ia(:), ja(:), d(:)
  real(psb_spk_),    intent(in)  :: val(:)
  real(psb_spk_),    intent(out) :: data(nrd,ncd)
  integer(psb_ipk_), intent(out) :: info
  logical, intent(in), optional  :: initdata
  integer(psb_ipk_), intent(in), optional :: rdisp

  !locals
  logical                        :: initdata_
  integer(psb_ipk_) :: rdisp_
  integer(psb_ipk_) :: i,ir,ic,k
  logical, parameter :: debug=.false.

  info = psb_success_
  initdata_ = .true.
  if (present(initdata)) initdata_ = initdata
  rdisp_ = 0
  if (present(rdisp)) rdisp_ = rdisp

  if (debug) write(0,*) 'Start xtr_dia_from_coo',nr,nc,nz,nrd,ncd,initdata_, rdisp_

  if (initdata_) data(1:nrd,1:ncd) = szero

  do i=1,nz
    ir = ia(i) 
    k  = ja(i) - ir
    ic = d(nr+k)
    if (debug) write(0,*) 'loop xtr_dia_from_coo :',ia(i),ja(i),k,ir-rdisp_,ic
    data(ir-rdisp_,ic) = val(i)
  enddo


end subroutine psi_s_xtr_dia_from_coo
