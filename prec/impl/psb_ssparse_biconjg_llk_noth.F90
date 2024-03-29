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
subroutine psb_ssparse_biconjg_llk_noth(n,a,p,z,w,nzrmax,sp_thresh,info)
  use psb_base_mod
  use psb_s_ainv_tools_mod
  use psb_s_biconjg_mod, psb_protect_name => psb_ssparse_biconjg_llk_noth

  !
  ! Left-looking variant, with NO drop rule on p(i)/p(j)
  !
  !
  implicit none
  integer(psb_ipk_), intent(in)               :: n
  type(psb_s_csr_sparse_mat), intent(in)    :: a
  type(psb_s_csc_sparse_mat), intent(inout) :: z,w
  integer(psb_ipk_), intent(in)               :: nzrmax
  real(psb_spk_), intent(in)                :: sp_thresh
  real(psb_spk_), intent(out)                :: p(:)
  integer(psb_ipk_), intent(out)              :: info

  ! Locals
  integer(psb_ipk_), allocatable        :: ia(:), ja(:), izkr(:), izcr(:)
  real(psb_spk_), allocatable :: zval(:),val(:), q(:)
  integer(psb_ipk_)  :: i,j,k, kc, kr, err_act, nz, nzra, nzrz, ipzi,ipzj,&
       & nzzi,nzzj, nzz, ip1, ip2, ipza,ipzz, ipzn, nzzn, ipz1, ipz2,&
       &  ipj, lastj, nextj, nzw
  type(psb_i_heap)   :: heap, rheap
  type(psb_s_csc_sparse_mat) :: ac
  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='psb_orth_llk_noth'
  logical, parameter :: debug=.false.

  allocate(zval(n),ia(n),val(n),izkr(n),izcr(n),q(n),stat=info)
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
    izkr(i) = 0
    izcr(i) = 0
    zval(i)  = szero
  end do

  ! Init z_1=e_1 and p_1=a_11
  p(1) = szero
  i   = 1
  nz  = a%irp(i+1) - a%irp(i)
  do j=1,nz
    if (a%ja(j) == 1) then
      p(1) = a%val(j)
      exit
    end if
  end do
  if (abs(p(1)) < s_epstol) &
       & p(1) = 1.d-3

  q(1) = p(1)
  !
  !
  call z%allocate_mnnz(n,n,n*nzrmax)

  z%icp(1)  = 1
  z%icp(2)  = 2
  z%ia(1)  = 1
  z%val(1) = done
  nzz       = 1

  call w%allocate_mnnz(n,n,n*nzrmax)
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
    ! !$          zval(j) = szero
    ! !$        end do
    zval(i)  = done
    izkr(i) = 1
    call heap%init(info)
    if (info == psb_success_) call heap%insert(i,info)

    if (info == psb_success_) call rheap%init(info)
    do j = ac%icp(i), ac%icp(i+1)-1
      if (ac%ia(j) < i) then
        if (info == psb_success_) call rheap%insert(ac%ia(j),info)
        izcr(ac%ia(j)) = 1
      end if
    end do
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_init_heap')
      return
    end if

    ! Update loop
    ! The idea is to keep track of the indices of the nonzeros in zval,
    ! so as to only do the dot products on the rows which have nonzeros
    ! in their positions; to do this we keep an extra
    ! copy of A in CSC, and the row indices to be considered are in rheap.
    lastj = -1
    outer: do
      inner: do
        call rheap%get_first(j,info)
        if (debug) write(0,*) 'from get_first: ',j,info
        if (info == -1) exit outer ! Empty heap
        if (j > lastj) then
          lastj = j
          exit inner
        end if
      end do inner

      izcr(j) = 0
      if (j>=i) cycle outer
      if (debug) write(0,*) 'update loop, using row: ',j,i
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

      if (.true.) then
        do k=z%icp(j), z%icp(j+1)-1
          kr     = z%ia(k)
          zval(kr) = zval(kr) + alpha*z%val(k)
          if (izkr(kr) == 0) then

            call heap%insert(kr,info)
            if (info /= psb_success_) exit
            izkr(kr) = 1
            ! We have just added a new nonzero in KR. Thus, we will
            ! need to explicitly compute the dot products on all
            ! rows j<k<i with nonzeros in column kr; we keep  them in
            ! a heap.
            !
            do kc = ac%icp(kr), ac%icp(kr+1)-1
              nextj=ac%ia(kc)
              if ((info == psb_success_).and.(izcr(nextj)==0)&
                   & .and.(nextj>j).and.(nextj<i)) then
                call rheap%insert(nextj,info)
                izcr(nextj) = 1
              end if
            end do
            if (debug) write(0,*) 'update loop, adding indices: ',&
                 &  ac%ia(ac%icp(kr):ac%icp(kr+1)-1)

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
    !call a%csget(i,i,nzra,ia,ja,val,info)
    call psb_s_csr_csgetrow(i,i,a,nzra,ia,ja,val,info)
    call rwclip(nzra,ia,ja,val,ione,n,ione,n)
    p(i) = psb_spge_dot(nzra,ja,val,zval)
    if (abs(p(i)) < s_epstol) &
         & p(i) = 1.d-3

!!$      write(0,*) 'Dropping from a column with: ',i,psb_howmany_heap(heap),sp_thresh

    !
    ! Sparsify current ZVAL and put into ZMAT
    !
    call sparsify(i,nzrmax,sp_thresh,n,zval,nzrz,ia,val,info,iheap=heap,ikr=izkr)
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
    ! !$          zval(j) = szero
    ! !$        end do
    zval(i)  = done
    izkr(i) = 1
    call heap%init(info)
    if (info == psb_success_) call heap%insert(i,info)

    if (info == psb_success_) call rheap%init(info)
    do j = a%irp(i), a%irp(i+1)-1
      if (a%ja(j) < i) then
        if (info == psb_success_) call rheap%insert(a%ja(j),info)
        izcr(a%ja(j)) = 1
      end if
    end do
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_init_heap')
      return
    end if

    ! Update loop
    ! The idea is to keep track of the indices of the nonzeros in zval,
    ! so as to only do the dot products on the rows which have nonzeros
    ! in their positions; to do this we keep an extra
    ! copy of A in CSC, and the row indices to be considered are in rheap.
    lastj = -1
    outerw: do
      innerw: do
        call rheap%get_first(j,info)
        if (debug) write(0,*) 'from get_first: ',j,info
        if (info == -1) exit outerw ! Empty heap
        if (j > lastj) then
          lastj = j
          exit innerw
        end if
      end do innerw
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
      if (.true.) then

        do k=w%icp(j), w%icp(j+1)-1
          kr     = w%ia(k)
          zval(kr) = zval(kr) + alpha*w%val(k)
          if (izkr(kr) == 0) then
            call heap%insert(kr,info)
            if (info /= psb_success_) exit
            izkr(kr) = 1
            ! We have just added a new nonzero in KR. Thus, we will
            ! need to explicitly compute the dot products on all
            ! rows j<k<i with nonzeros in column kr; we keep  them in
            ! a heap.
            !
            do kc = a%irp(kr), a%irp(kr+1)-1
              nextj=a%ja(kc)
              if ((info == psb_success_).and.(izcr(nextj)==0)&
                   & .and.(nextj>j).and.(nextj<i)) then
                call rheap%insert(nextj,info)
                izcr(nextj) = 1
              end if
            end do
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
    if (abs(q(i)) < s_epstol) &
         & q(i) = 1.d-3

!!$      write(0,*) 'Dropping from a column with: ',i,psb_howmany_heap(heap),sp_thresh
    !
    ! Sparsify current ZVAL and put into ZMAT
    !
    call sparsify(i,nzrmax,sp_thresh,n,zval,nzrz,ia,val,info,iheap=heap,ikr=izkr)
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

end subroutine psb_ssparse_biconjg_llk_noth
