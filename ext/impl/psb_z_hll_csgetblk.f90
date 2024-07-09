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
  

subroutine psb_z_hll_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  use psb_base_mod
  use psb_z_hll_mat_mod, psb_protect_name => psb_z_hll_csgetblk
  implicit none

  class(psb_z_hll_sparse_mat), intent(in)    :: a
  class(psb_z_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                :: imin,imax
  integer(psb_ipk_),intent(out)                :: info
  logical, intent(in), optional                :: append
  integer(psb_ipk_), intent(in), optional      :: iren(:)
  integer(psb_ipk_), intent(in), optional      :: jmin,jmax
  logical, intent(in), optional                :: rscale,cscale
  Integer(Psb_ipk_)  :: err_act, nzin, nzout
  character(len=20)  :: name='hll_getblk'
  logical            :: append_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(append)) then 
    append_ = append
  else
    append_ = .false.
  endif
  if (append_) then 
    nzin = a%get_nzeros()
  else
    nzin = 0
  endif

  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
       & nzin=nzin, rscale=rscale, cscale=cscale)

  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%set_host()
  call b%fix(info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_z_hll_csgetblk
