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
subroutine psb_dsparse_biconjg_mlk(n,a,p,z,w,nzrmax,sp_thresh,info)
  use psb_base_mod
  use psb_ainv_tools_mod
  use psb_d_biconjg_mod, psb_protect_name => psb_dsparse_biconjg_mlk
  !
  ! Left-looking variant
  !
  !
  implicit none
  integer(psb_ipk_), intent(in)             :: n
  type(psb_d_csr_sparse_mat), intent(in)    :: a
  type(psb_d_csc_sparse_mat), intent(inout) :: z,w
  integer(psb_ipk_), intent(in)             :: nzrmax
  real(psb_dpk_), intent(in)                :: sp_thresh
  real(psb_dpk_), intent(out)               :: p(:)
  integer(psb_ipk_), intent(out)            :: info

  ! Locals
  integer(psb_ipk_), allocatable :: ia(:), ja(:), izkr(:), izcr(:), hlist(:), bfr(:), rwlist(:)
  real(psb_dpk_), allocatable    :: zval(:),val(:), q(:)
  integer(psb_ipk_) :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
       & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn, ipz1, ipz2,&
       &  ipj, lastj, nextj, nzw, hlhead, li, mj, kkc, ifrst, ilst, rwhead
  type(psb_i_heap)   :: heap, rheap
  type(psb_d_csc_sparse_mat) :: ac
  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='psb_biconjg_mlk'
  logical, parameter :: debug=.false., test_merge=.true.

  allocate(zval(n),ia(n),val(n),izkr(n),izcr(n),q(n), &
       & hlist(n),rwlist(n),bfr(n),stat=info)
  if (info == psb_success_) call ac%cp_from_fmt(a,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    return
  end if
  !
  ! izkr(i): flag nonzeros in ZVAL. To minimize traffic into heap.
  ! izcr(i): flag rows to be used for the dot products. Used to minimize
  !               traffic in rheap.
  !
  do i=1,n
    izkr(i)  = 0
    izcr(i)  = 0
    zval(i)  = dzero
    hlist(i) = -1
    rwlist(i) = -1
  end do

  ! Init z_1=e_1 and p_1=a_11
  p(1) = dzero
  i   = 1
  nz  = a%irp(i+1) - a%irp(i)
  do j=1,nz
    if (a%ja(j) == 1) then
      p(1) = a%val(j)
      exit
    end if
  end do
  if (abs(p(1)) < d_epstol) &
       & p(1) = 1.d-3

  q(1) = p(1)
  !
  !
  call z%allocate(n,n,n*nzrmax)

  z%icp(1)  = 1
  z%icp(2)  = 2
  z%ia(1)  = 1
  z%val(1) = done
  nzz       = 1

  call w%allocate(n,n,n*nzrmax)
  w%icp(1)  = 1
  w%icp(2)  = 2
  w%ia(1)  = 1
  w%val(1) = done
  nzw       = 1


  do i = 2, n
    if (debug) write(0,*) 'Main loop iteration ',i,n

    !
    ! Update loop on Z.
    ! Must be separated from update loop of W because of
    ! the conflict on J that would result.
    !

    ! ZVAL = e_i
    ! !$        do j=1, i-1
    ! !$          zval(j) = dzero
    ! !$        end do
    zval(i)  = done
    izkr(i) = 1
    rwhead  = i

    hlhead = -1

    kkc = 0
    ilst  = ac%icp(i)-1
    ifrst = ac%icp(i)
    do j = ac%icp(i+1)-1, ac%icp(i), -1
      if (ac%ia(j) < i) then
        ilst = j
        exit
      end if
    end do
    kkc = ilst-ifrst+1

    if (.true..or.debug) then
!!$      write(0,*) 'Outer Before insertion  : ',hlhead
      call printlist(hlhead,hlist)
    end if
    if (kkc > 0) then
!!$      write(0,*) i,' Outer Inserting : ',kkc,':',ac%ia(ifrst:ilst)

      !call hlmerge(hlhead,hlist,bfr(1:kkc))
      call hlmerge(hlhead,hlist,ac%ia(ifrst:ilst))
    end if
    if (debug) then
      write(0,*) 'Outer After insertion: ',hlhead
      call printlist(hlhead,hlist)
    end if

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='init_lists')
      return
    end if

    ! Update loop
    ! The idea is to keep track of the indices of the nonzeros in zval,
    ! so as to only do the dot products on the rows which have nonzeros
    ! in their positions; to do this we keep an extra
    ! copy of A in CSC, and the row indices to be considered are in rheap.
    lastj = -1
    outer: do
      mj      = hlhead
      if (mj > 0)  then
        hlhead = hlist(mj)
        hlist(mj) = -1
      end if
      j = mj
      if (j < 0)  exit outer

      izcr(j) = 0
      if (j>=i) cycle outer

      if (debug) write(0,*) 'update loop, using row: ',j,i,mj
      ip1 = a%irp(j)
      ip2 = a%irp(j+1) - 1
      do
        if (ip2 < ip1 ) exit
        if (a%ja(ip2) <= n) exit
        ip2 = ip2 -1
      end do
      nzra = max(0,ip2 - ip1 + 1)
      p(i) = psb_spge_dot(nzra,a%ja(ip1:ip2),a%val(ip1:ip2),zval)
      ! !$          write(psb_err_unit,*) j,i,p(i)

      alpha = (-p(i)/p(j))
!!$      write(0,*) 'At step ',i,j,' p(i) ',p(i),alpha
!!$      write(0,*) '   Current list is : ',hlhead
!!$      call printlist(hlhead,hlist)
!!$


      if (.false..or.(abs(alpha) > sp_thresh)) then
        ifrst=z%icp(j)
        ilst=z%icp(j+1)-1
        call hlmerge(rwhead,rwlist,z%ia(ifrst:ilst))
!!$        write(0,*) 'At step ',i,j,' range ',z%icp(j), z%icp(j+1)-1, &
!!$             & ' vals ',z%ia(z%icp(j):z%icp(j+1)-1)

        do k=z%icp(j), z%icp(j+1)-1
          kr     = z%ia(k)
          zval(kr) = zval(kr) + alpha*z%val(k)

          if (izkr(kr) == 0) then
!!$            write(0,*) ' main inner Inserting ',kr
!!$            call hlmerge(rwhead,rwlist,(/kr/))
            izkr(kr) = 1
            ! We have just added a new nonzero in KR. Thus, we will
            ! need to explicitly compute the dot products on all
            ! rows j<k<i with nonzeros in column kr; we keep  them in
            ! a heap.
            !
            ilst  = ac%icp(kr)-1
            ifrst = ac%icp(kr+1)
            kkc = 0
            do kc = ac%icp(kr), ac%icp(kr+1)-1
              if ((ac%ia(kc) < i).and.(ac%ia(kc) >j)) then
                ifrst = min(ifrst,kc )
                ilst  = max(ilst,kc)
              end if
            end do
            kkc = ilst-ifrst+1
            if (debug) then
              write(0,*) 'Inner Before insertion: '
              call printlist(hlhead,hlist)
              write(0,*) 'Inner Inserting : ',kkc,':',ac%ia(ifrst:ilst)
            end if
            if (ilst >= ifrst) then
!!$              write(0,*) j,i,' Inner inserting ',ac%ia(ifrst:ilst)
              call hlmerge(hlhead,hlist,ac%ia(ifrst:ilst))
            end if

            if (debug) then
              write(0,*) 'Inner After insertion: ',hlhead
              call printlist(hlhead,hlist)
            end if

          end if
          if (info /= psb_success_) exit
        end do
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_insert_heap')
          return
        end if
      end if
    end do outer
    call a%csget(i,i,nzra,ia,ja,val,info)
    call rwclip(nzra,ia,ja,val,ione,n,ione,n)
    p(i) = psb_spge_dot(nzra,ja,val,zval)
    if (abs(p(i)) < d_epstol) &
         & p(i) = 1.d-3

!!$      write(0,*) 'Dropping from a column with: ',i,psb_howmany_heap(heap),sp_thresh

    !
    ! Sparsify current ZVAL and put into ZMAT
    !
    call sparsify(i,nzrmax,sp_thresh,n,zval,nzrz,ia,val,rwhead,rwlist,izkr,info)

    if (info /= psb_success_) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='sparsify')
      return
    end if
    call psb_ensure_size(nzz+nzrz, z%ia,  info)
    call psb_ensure_size(nzz+nzrz, z%val, info)
    ipz1 = z%icp(i)
    do j=1, nzrz
      z%ia(ipz1  + j -1) = ia(j)
      z%val(ipz1 + j -1) = val(j)
    end do
    z%icp(i+1) = ipz1 + nzrz
    nzz        = nzz + nzrz


    ! WVAL = e_i
    ! !$        do j=1, i-1
    ! !$          zval(j) = dzero
    ! !$        end do
    zval(i)  = done
    izkr(i) = 1
    rwhead  = i
    hlhead = -1

    kkc = 0
    ilst  = a%irp(i)-1
    ifrst = a%irp(i)
    do j = a%irp(i+1)-1, a%irp(i), -1
      if (a%ja(j) < i) then
        ilst = j
        exit
      end if
    end do
    kkc = ilst-ifrst+1

    if (debug) then
      write(0,*) 'Outer Before insertion: '
      call printlist(hlhead,hlist)
      write(0,*) 'Outer Inserting : ',kkc,':',a%ja(ifrst:ilst)
    end if
    if (kkc > 0 ) then
      !call hlmerge(hlhead,hlist,bfr(1:kkc))
      call hlmerge(hlhead,hlist,a%ja(ifrst:ilst))
    end if
    if (debug) then
      write(0,*) 'Outer After insertion: ',hlhead
      call printlist(hlhead,hlist)
    end if

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='init_lists')
      return
    end if

    ! Update loop
    ! The idea is to keep track of the indices of the nonzeros in zval,
    ! so as to only do the dot products on the rows which have nonzeros
    ! in their positions; to do this we keep an extra
    ! copy of A in CSC, and the row indices to be considered are in rheap.
    lastj = -1
    outerw: do
      mj      = hlhead
      if (hlhead > 0)  then
        hlhead = hlist(mj)
        hlist(mj) = -1
      end if
      j = mj
      if (j < 0)  exit outerw

      izcr(j) = 0
      if (j>=i) cycle outerw
      if (debug) write(0,*) 'update loop, using row: ',j
      ip1 = ac%icp(j)
      ip2 = ac%icp(j+1) - 1
      do
        if (ip2 < ip1 ) exit
        if (ac%ia(ip2) <= n) exit
        ip2 = ip2 -1
      end do
      nzra = max(0,ip2 - ip1 + 1)
      q(i) = psb_spge_dot(nzra,ac%ia(ip1:ip2),ac%val(ip1:ip2),zval)
      ! !$          write(psb_err_unit,*) j,i,p(i)

      alpha = (-q(i)/q(j))
      if (.false..or.(abs(alpha) > sp_thresh)) then
        ifrst=w%icp(j)
        ilst=w%icp(j+1)-1
        call hlmerge(rwhead,rwlist,w%ia(ifrst:ilst))

        do k=w%icp(j), w%icp(j+1)-1
          kr     = w%ia(k)
          zval(kr) = zval(kr) + alpha*w%val(k)
          if (izkr(kr) == 0) then
            izkr(kr) = 1
            ! We have just added a new nonzero in KR. Thus, we will
            ! need to explicitly compute the dot products on all
            ! rows j<k<i with nonzeros in column kr; we keep  them in
            ! a heap.
            !
            ilst  = a%irp(kr)-1
            ifrst = a%irp(kr+1)
            kkc = 0
            do kc = a%irp(kr), a%irp(kr+1)-1
              if ((a%ja(kc) < i).and.(a%ja(kc) >j)) then
                ifrst = min(ifrst,kc )
                ilst  = max(ilst,kc)
              end if
            end do
            kkc = ilst-ifrst+1
            if (debug) then
              write(0,*) 'Inner Before insertion: '
              call printlist(hlhead,hlist)
              write(0,*) 'Inner Inserting : ',kkc,':',a%ja(ifrst:ilst)
            end if

            call hlmerge(hlhead,hlist,a%ja(ifrst:ilst))

            if (debug) then
              write(0,*) 'Inner After insertion: ',hlhead
              call printlist(hlhead,hlist)
            end if
            if (debug) write(0,*) 'update loop, adding indices: ',&
                 &  a%ja(a%irp(kr):a%irp(kr+1)-1)

          end if
          if (info /= psb_success_) exit
        end do
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_insert_heap')
          return
        end if
      end if
    end do outerw
    ip1 = ac%icp(i)
    ip2 = ac%icp(i+1) - 1
    do
      if (ip2 < ip1 ) exit
      if (ac%ia(ip2) <= n) exit
      ip2 = ip2 -1
    end do
    nzra = max(0,ip2 - ip1 + 1)
    q(i) = psb_spge_dot(nzra,ac%ia(ip1:ip2),ac%val(ip1:ip2),zval)
    if (abs(q(i)) < d_epstol) &
         & q(i) = 1.d-3

!!$      write(0,*) 'Dropping from a column with: ',i,psb_howmany_heap(heap),sp_thresh
    !
    ! Sparsify current ZVAL and put into ZMAT
    !
    call sparsify(i,nzrmax,sp_thresh,n,zval,nzrz,ia,val,rwhead,rwlist,izkr,info)
    if (info /= psb_success_) then
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='sparsify')
      return
    end if
    call psb_ensure_size(nzw+nzrz, w%ia,  info)
    call psb_ensure_size(nzw+nzrz, w%val, info)
    ipz1 = w%icp(i)
    do j=1, nzrz
      w%ia(ipz1  + j -1) = ia(j)
      w%val(ipz1 + j -1) = val(j)
    end do
    w%icp(i+1) = ipz1 + nzrz
    nzw        = nzw + nzrz

  end do

contains

  subroutine hlmerge(head,listv,vals)
    integer(psb_ipk_), intent(inout) :: head, listv(:)
    integer(psb_ipk_), intent(in)    :: vals(:)
    integer(psb_ipk_) :: i,j,k, lh, lv, nv, vv, flh, ph

    nv  = size(vals)
    lh  = head
    flh = -1
    lv  = 1
    if ((head < 0).and.(nv > 0)) then
      ! Adjust head if empty
      head = vals(lv)
      lv = lv + 1
    else if ((head > 0) .and. (nv >0)) then
      ! Adjust head if first item less than it
      if (head > vals(lv)) then
        listv(vals(lv)) = head
        head = vals(lv)
        lv   = lv + 1
      end if
    end if

    lh = head
    ph = lh
    do while ((lh > 0) .and. (lv <= nv))
      if (lh == vals(lv)) then
        lv = lv + 1
      else if (lh > vals(lv)) then
        listv(vals(lv)) = lh
        listv(ph) = vals(lv)
        lh = vals(lv)
        lv = lv + 1
      else
        ph = lh
        lh = listv(lh)
      end if
    end do
    lh = ph
    do while (lv <= nv)
      listv(lh) = vals(lv)
      lh = listv(lh)
      lv = lv + 1
    end do
  end subroutine hlmerge


  subroutine printlist(head,listv)
    integer(psb_ipk_), intent(in) :: head, listv(:)
    integer(psb_ipk_) :: li

    li = head
    do while (li > 0)
      write(0,*) 'Item: ', li
      li = listv(li)
    end do
  end subroutine printlist

end subroutine psb_dsparse_biconjg_mlk
