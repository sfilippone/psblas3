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
subroutine amg_c_sparsify(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,info,istart,iheap,ikr)
  use psb_base_mod
  implicit none

  real(psb_spk_), intent(in)              :: sp_thresh
  integer(psb_ipk_), intent(in)           :: idiag, n, nzrmax
  complex(psb_spk_), intent(inout)           :: zw(:)
  integer(psb_ipk_), intent(out)          :: nz
  integer(psb_ipk_), intent(out)          :: iz(:)
  complex(psb_spk_), intent(out)             :: valz(:)
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: istart
  type(psb_i_heap), optional              :: iheap
  integer(psb_ipk_), optional             :: ikr(:)
  !
  integer(psb_ipk_)              :: i, istart_, last_i, iret,k
  complex(psb_spk_)                :: witem
  integer(psb_ipk_)              :: widx
  complex(psb_spk_), allocatable   :: xw(:)
  integer(psb_ipk_), allocatable :: xwid(:), indx(:)
  type(psb_c_idx_heap)         :: heap


  info = psb_success_
  istart_ = 1
  if (present(istart)) istart_ = max(1,istart)
  if (.false.) then
    nz = 0
    do i=istart_, n
      if ((i == idiag).or.(abs(zw(i)) >= sp_thresh)) then
        nz       = nz + 1
        iz(nz)   = i
        valz(nz) = zw(i)
      end if
    end do

  else

    allocate(xw(nzrmax),xwid(nzrmax),indx(nzrmax),stat=info)
    if (info /= psb_success_) then
      return
    end if

    call heap%init(info,dir=psb_asort_down_)

    ! Keep at least the diagonal
    nz = 0

    if (present(iheap)) then
      if (.not.(present(ikr))) then
        write(psb_err_unit,*) 'Error: if IHEAP then also IKR'
        info = -1
        return
      end if
      last_i = -1
      do
        call iheap%get_first(i,iret)
        if (iret < 0) exit
        ! An index may have been put on the heap more than once.
        if (i == last_i) cycle
        last_i = i
        if (i == idiag) then
          xw(1)   = zw(i)
          xwid(1) = i
        else if (abs(zw(i)) >= sp_thresh) then
          call heap%insert(zw(i),i,info)
        end if
        zw(i)  = dzero
        ikr(i) = 0
      end do

    else

      do i=istart_, n
        if (i == idiag) then
          xw(1)   = zw(i)
          xwid(1) = i
        else if (abs(zw(i)) >= sp_thresh) then
          call heap%insert(zw(i),i,info)
        end if
        zw(i) = dzero
      end do
    end if

    k = 1
    do
      if (k == nzrmax) exit
      call heap%get_first(witem,widx,info)
      if (info == -1) then
        info = psb_success_
        exit
      endif
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        return
      end if
      k = k + 1
      xw(k)   = witem
      xwid(k) = widx
    end do
    call heap%free(info)
    nz = k
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)
!!$    write(0,*) 'sparsify output for idiag ',idiag,' :',nz,sp_thresh
    do i=1, nz
      valz(i) = xw(indx(i))
      iz(i)   = xwid(i)
!!$      write(0,*) '         ',iz(i),valz(i)
    end do

  end if

  return

end subroutine amg_c_sparsify


subroutine amg_c_sparsify_list(idiag,nzrmax,sp_thresh,n,zw,nz,iz,valz,lhead,listv,ikr,info)
  use psb_base_mod
  implicit none

  real(psb_spk_), intent(in)      :: sp_thresh
  integer(psb_ipk_), intent(in)     :: idiag, n, nzrmax
  complex(psb_spk_), intent(inout)    :: zw(:)
  integer(psb_ipk_), intent(out)    :: nz
  integer(psb_ipk_), intent(out)    :: iz(:)
  complex(psb_spk_), intent(out)       :: valz(:)
  integer(psb_ipk_), intent(out)    :: info
  integer(psb_ipk_), intent(inout)  :: lhead, listv(:)
  integer(psb_ipk_)                 :: ikr(:)
  !
  integer(psb_ipk_)              :: i, istart_, last_i, iret,k,current, next
  complex(psb_spk_)                :: witem
  integer(psb_ipk_)              :: widx
  complex(psb_spk_), allocatable    :: xw(:)
  integer(psb_ipk_), allocatable :: xwid(:), indx(:)


  info = psb_success_
  istart_ = 1
  allocate(xw(n),xwid(n),indx(n),stat=info)

  current = lhead
  lhead = -1
  i = 0
  do while (current >0)
    i = i + 1
    xw(i)   = zw(current)
    xwid(i) = current

    if (current == idiag) then
      ! Bring the diagona into first position
      witem = xw(1)
      widx  = xwid(1)
      xw(1)   = xw(i)
      xwid(1) = xwid(i)
      xw(i)   = witem
      xwid(i) = widx
    end if

    zw(current)  = dzero
    ikr(current) = 0

    next = listv(current)
    listv(current) = -1
    current = next
  end do
  nz = i
  if (nz > 2) call psb_hsort(xw(2:nz),ix=xwid(2:nz),&
       & dir=psb_asort_down_,flag=psb_sort_keep_idx_)
!!$  write(0,*) 'Done first msort '
!!$  write(0,*) '   after first msort for idiag ',idiag,' :',nz,sp_thresh
!!$  do i=1, nz
!!$    write(0,*) '         ',xwid(i),xw(i)
!!$  end do

  i = 2
  do while (i<=nz)
    if (abs(xw(i)) < sp_thresh) exit
    i = i + 1
  end do
!!$  write(0,*) 'NZ ',nz, i, nzrmax
  nz = max(1,min(i-1,nzrmax))
  call psb_msort(xwid(1:nz),ix=indx(1:nz),dir=psb_sort_up_)
!!$  write(0,*) 'Done second msort '

!!$  write(0,*) 'sparsify output for idiag ',idiag,' :',nz,i,sp_thresh
  do i=1, nz
    valz(i) = xw(indx(i))
    iz(i)   = xwid(i)
!!$      write(0,*) '         ',iz(i),valz(i)
  end do

  return

end subroutine amg_c_sparsify_list
