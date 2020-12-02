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
subroutine psb_c_sp_drop(idiag,nzrmax,sp_thresh,nz,iz,valz,info)
  use psb_base_mod
  implicit none
  real(psb_spk_), intent(in)       :: sp_thresh
  integer(psb_ipk_), intent(in)    :: idiag, nzrmax
  integer(psb_ipk_), intent(inout) :: nz
  integer(psb_ipk_), intent(inout) :: iz(:)
  complex(psb_spk_), intent(inout)    :: valz(:)
  integer(psb_ipk_), intent(out)   :: info
  !
  integer(psb_ipk_)              :: i, j, idf, nw
  complex(psb_spk_)                 :: witem
  integer(psb_ipk_)              :: widx
  complex(psb_spk_), allocatable    :: xw(:)
  integer(psb_ipk_), allocatable :: xwid(:), indx(:)


  info = psb_success_

  if (nz > min(size(iz),size(valz))) then
    write(0,*) 'Serious size problem ',nz,size(iz),size(valz)
    info = -2
    return
  end if
  allocate(xw(nz),xwid(nz),indx(nz),stat=info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) ' Memory allocation failure in sp_drop',nz,info
    return
  endif

  ! Always keep the diagonal element
  idf = -1
  do i=1, nz
    if (iz(i) == idiag) then
      idf     = i
      witem   = valz(i)
      widx    = iz(i)
      valz(i) = valz(1)
      iz(i)   = iz(1)
      valz(1) = witem
      iz(1)   = widx
      exit
    end if
  end do

  if (idf == -1) then

    xw(1:nz) = valz(1:nz)
    call psb_qsort(xw(1:nz),indx(1:nz),dir=psb_asort_down_)
    do i=1, nz
      xwid(i) = iz(indx(i))
    end do
    nw = min(nz,nzrmax)
    do
      if (nw <= 1) exit
      if (abs(xw(nw)) < sp_thresh) then
        nw = nw - 1
      else
        exit
      end if
    end do
    nw = max(nw, 1)

  else

    nw = nz-1

    xw(1:nw) = valz(2:nz)

    call psb_qsort(xw(1:nw),indx(1:nw),dir=psb_asort_down_)
    nw = min(nw,nzrmax-1)
    do
      if (nw <= 1) exit
      if (abs(xw(nw)) < sp_thresh) then
        nw = nw - 1
      else
        exit
      end if
    end do

    do i=1, nw
      xwid(i) = iz(1+indx(i))
    end do
    nw       = nw + 1
    xw(nw)   = valz(1)
    xwid(nw) = iz(1)
  end if

  call psb_msort(xwid(1:nw),indx(1:nw),dir=psb_sort_up_)

  do i=1, nw
    valz(i) = xw(indx(i))
    iz(i)   = xwid(i)
  end do
  nz = nw
  if (nz>nzrmax) write(0,*) 'in sp_drop: ',nw,nzrmax,nz
  deallocate(xw,xwid,indx,stat=info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) ' Memory deallocation failure in sp_drop',info
    return
  endif
  return
end subroutine psb_c_sp_drop
