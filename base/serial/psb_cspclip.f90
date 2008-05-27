!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$   
! File:  psb_cspclip.f90 
! Subroutine: psb_cspclip
!    Creates a "clipped" copy of input matrix A. Output is always in COO. 
! Arguments:

!*****************************************************************************
!*                                                                           *
!*****************************************************************************
subroutine psb_cspclip(a,b,info,imin,imax,jmin,jmax,rscale,cscale)
  use psb_spmat_type
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_cspclip
  implicit none 
  type(psb_cspmat_type), intent(in)  :: a
  type(psb_cspmat_type), intent(out) :: b
  integer, intent(out)               :: info
  integer, intent(in), optional      :: imin,imax,jmin,jmax
  logical, intent(in), optional      :: rscale,cscale

  integer               :: err_act
  character(len=20)     :: name
  integer      :: imin_,imax_,jmin_,jmax_
  logical      :: rscale_,cscale_
  integer      :: sizeb, nzb, mb, kb, ifst, ilst, nrt, nzt, i, j
  integer, parameter :: irbk=40

  name='psb_csp_clip'
  info  = 0
  call psb_erractionsave(err_act)
  call psb_set_erraction(0)  

  call psb_nullify_sp(b)


  if (present(imin)) then 
    imin_ = imin
  else
    imin_ = 1
  endif
  if (present(imax)) then 
    imax_ = imax
  else
    imax_ = a%m
  endif
  if (present(jmin)) then 
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then 
    jmax_ = jmax
  else
    jmax_ = a%k
  endif
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .true.
  endif
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .true.
  endif
  
  if (rscale_) then 
    mb = imax_ - imin_ +1
  else 
    mb = a%m  ! Should this be imax_ ?? 
  endif
  if (cscale_) then 
    kb = jmax_ - jmin_ +1
  else 
    kb = a%k  ! Should this be jmax_ ?? 
  endif


  sizeb = psb_sp_get_nnzeros(a)
  call psb_sp_all(mb,kb,b,sizeb,info)
  b%fida   = 'COO'
  b%descra = a%descra
  nzb = 0 
  do i=imin_, imax_, irbk
    nrt = min(irbk,imax_-i+1)
    ifst = i
    ilst = ifst + nrt - 1
    call psb_sp_getrow(ifst,a,nzt,b%ia1,b%ia2,b%aspk,info,&
         &  lrw=ilst, append=.true.,nzin=nzb)
    do j=nzb+1,nzb+nzt 
      if ((jmin_ <= b%ia2(j)).and.(b%ia2(j) <= jmax_)) then 
        nzb         = nzb + 1
        b%aspk(nzb) = b%aspk(j) 
        b%ia1(nzb)  = b%ia1(j) 
        b%ia2(nzb)  = b%ia2(j) 
      end if
    end do
    
  end do
  b%infoa(psb_nnz_) = nzb 

  if (rscale_) then 
    do i=1, nzb
      b%ia1(i) = b%ia1(i) - imin_ + 1
    end do
  end if
  if (cscale_) then 
    do i=1, nzb
      b%ia2(i) = b%ia2(i) - jmin_ + 1
    end do
  end if
  call psb_fixcoo(b,info) 
  call psb_sp_trim(b,info)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_cspclip

