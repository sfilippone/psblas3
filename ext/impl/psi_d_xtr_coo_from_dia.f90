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
  
subroutine psi_d_xtr_coo_from_dia(nr,nc,ia,ja,val,nz,nrd,ncd,data,offsets,info,rdisp)
  use psb_base_mod, only : psb_ipk_, psb_success_, psb_dpk_, dzero

  implicit none 

  integer(psb_ipk_), intent(in)    :: nr,nc, nrd,ncd, offsets(:) 
  integer(psb_ipk_), intent(inout) :: ia(:), ja(:),nz
  real(psb_dpk_),    intent(inout) :: val(:)
  real(psb_dpk_),    intent(in)    :: data(nrd,ncd)
  integer(psb_ipk_), intent(out)   :: info
  integer(psb_ipk_), intent(in), optional :: rdisp

  !locals
  integer(psb_ipk_) :: rdisp_, nrcmdisp, rdisp1
  integer(psb_ipk_) :: i,j,ir1, ir2, ir, ic,k
  logical, parameter :: debug=.false.

  info = psb_success_
  rdisp_ = 0
  if (present(rdisp)) rdisp_ = rdisp

  if (debug) write(0,*) 'Start xtr_coo_from_dia',nr,nc,nrd,ncd, rdisp_
  nrcmdisp = min(nr-rdisp_,nc-rdisp_) 
  rdisp1   = 1-rdisp_
  nz = 0 
  do j=1, ncd
    if (offsets(j)>=0) then 
      ir1 = 1
      ! ir2 = min(nrd,nr - offsets(j) - rdisp_,nc-offsets(j)-rdisp_)
      ir2 = min(nrd, nrcmdisp - offsets(j))
    else
      ! ir1 = max(1,1-offsets(j)-rdisp_) 
      ir1 = max(1, rdisp1 - offsets(j))
      ir2 = min(nrd, nrcmdisp) 
    end if
    if (debug) write(0,*) ' Loop  J',j,ir1,ir2, offsets(j)      
    do i=ir1,ir2 
      ir = i + rdisp_
      ic = i + rdisp_ + offsets(j)
      if (debug) write(0,*) ' Loop  I',i,ir,ic
      nz = nz + 1
      ia(nz) = ir
      ja(nz) = ic 
      val(nz) = data(i,j)
    enddo
  end do

end subroutine psi_d_xtr_coo_from_dia

