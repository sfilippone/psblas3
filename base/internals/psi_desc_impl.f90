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

subroutine psi_renum_index(iperm,idx,info)
  use psi_mod, psi_protect_name =>  psi_renum_index
  use psb_serial_mod
  implicit none

  integer(psb_ipk_), intent(out)   :: info
  integer(psb_ipk_), intent(in)    :: iperm(:)
  integer(psb_ipk_), intent(inout) :: idx(:)

  integer(psb_ipk_) :: i,j,k,nh

  i=1
  k=idx(i)
  do while (k /= -1)
    i = i+1
    nh = idx(i)
    do j = i+1, i+nh
      idx(j) = iperm(idx(j))
    enddo
    i  = i + nh + 1
    nh = idx(i)
    do j = i+1, i+nh
      idx(j) = iperm(idx(j))
    enddo
    i = i + nh + 1
    k = idx(i)
  enddo

end subroutine psi_renum_index

subroutine psi_i_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info, mold)

  use psi_mod, psi_protect_name =>  psi_i_cnv_dsc
  use psb_realloc_mod
  use mpi
  implicit none

  !     ....scalars parameters....
  integer(psb_ipk_), intent(in)      :: halo_in(:), ovrlap_in(:),ext_in(:)
  type(psb_desc_type), intent(inout) :: cdesc
  integer(psb_ipk_), intent(out)     :: info
  class(psb_i_base_vect_type), optional, intent(in) :: mold

  !     ....local scalars....
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: ictxt, err_act,nxch,nsnd,nrcv,j,k
  !     ...local array...
  integer(psb_ipk_), allocatable  :: idx_out(:), tmp_mst_idx(:)

  !     ...parameters
  integer(psb_ipk_) :: debug_level, debug_unit
  logical, parameter :: debug=.false.
  character(len=20)  :: name
  !     ...Artless
  integer(psb_ipk_) :: proc_to_comm, pnti, nerv, nesd, i, rcv_pt, snd_pt
  integer(psb_ipk_), allocatable, dimension(:) :: prcid, brvidx, rvsz, bsdidx
  integer(psb_ipk_), allocatable, dimension(:) :: iaux, sdsz
  integer(psb_ipk_), allocatable, dimension(:) :: max_degree, src_and_dest
  ! to remove: for mpix_init
  logical, parameter :: reorder=.FALSE., persistent_mpi=.TRUE.
  integer(psb_ipk_), allocatable, dimension(:) :: rcv_buf
  integer :: dist_graph_comm, degree, buf_size, req, ierr
  !     ...end Artless


  name='psi_cnv_desc'
  call psb_get_erraction(err_act)
  debug_level = psb_get_debug_level()
  debug_unit  = psb_get_debug_unit()

  info = psb_success_
  ictxt = cdesc%get_context()

  call psb_info(ictxt,me,np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif


  ! first the halo index
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on halo',&
       & size(halo_in)
  call psi_crea_index(cdesc,halo_in, idx_out,nxch,nsnd,nrcv,info)

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%halo_index,info)

  if (debug_level>0) write(debug_unit,*) me,'Done crea_index on halo'
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'


  ! then ext index
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'
  call psi_crea_index(cdesc,ext_in, idx_out,nxch,nsnd,nrcv,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ext_index,info)

  if (debug_level>0) write(debug_unit,*) me,'Done crea_index on ext'
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ovrlap'

  ! then the overlap index
  call psi_crea_index(cdesc,ovrlap_in, idx_out,nxch,nsnd,nrcv,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_index')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ovrlap_index,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_move_alloc')
    goto 9999
  end if


  ! next  ovrlap_elem
  if (debug_level>0) write(debug_unit,*) me,'Calling crea_ovr_elem'
  call psi_crea_ovr_elem(me,cdesc%ovrlap_index,cdesc%ovrlap_elem,info)
  if (debug_level>0) write(debug_unit,*) me,'Done crea_ovr_elem'
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_ovr_elem')
    goto 9999
  end if
  ! Extract ovr_mst_idx from ovrlap_elem
  if (debug_level>0) write(debug_unit,*) me,'Calling bld_ovr_mst'
  call psi_bld_ovr_mst(me,cdesc%ovrlap_elem,tmp_mst_idx,info)
  if (info == psb_success_) call psi_crea_index(cdesc,&
       & tmp_mst_idx,idx_out,nxch,nsnd,nrcv,info)
  if (debug_level>0 )write(debug_unit,*) me,'Done crea_indx'
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_bld_ovr_mst')
    goto 9999
  end if
  call psb_move_alloc(idx_out,cdesc%ovr_mst_idx,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_move_alloc')
    goto 9999
  end if

  ! finally bnd_elem
  call psi_crea_bnd_elem(idx_out,cdesc,info)
  if (info == psb_success_) call psb_move_alloc(idx_out,cdesc%bnd_elem,info)

  call cdesc%v_halo_index%bld(cdesc%halo_index,mold=mold)
  call cdesc%v_ext_index%bld(cdesc%ext_index,mold=mold)
  call cdesc%v_ovrlap_index%bld(cdesc%ovrlap_index,mold=mold)
  call cdesc%v_ovr_mst_idx%bld(cdesc%ovr_mst_idx,mold=mold)


  ! ARTLESS: start of additions
  ictxt = cdesc%get_ctxt() ! get context to then get
  call psb_info(ictxt,me,np) ! rank and number of processors
  ! print *, "---artless:psi_i_cnv_dsc, me=",me," np=",np
  iaux = cdesc%v_halo_index%get_vect() ! iaux is allocatable int array
  ! print *, "---artless:psi_i_cnv_dsc:",me," iaux = ", iaux
  ! psb_get_mpicomm  with integer(psb_mpk_) becaust standard says 4 byte
  ! copying from base/internals/psi_dswapdata_a.F90
  if (persistent_mpi) then ! artless: make a proper flag
    allocate(prcid(0:np-1), brvidx(0:np-1), rvsz(0:np-1), bsdidx(0:np-1))
    allocate(sdsz(0:np-1))
    allocate(max_degree(np-1)) ! max degree of a node is num ranks - 1
    degree = 0
    pnti   = 1
    snd_pt = 1
    rcv_pt = 1
    do while (iaux(pnti+psb_proc_id_) .ne. -1)
      degree = degree + 1
      proc_to_comm = iaux(pnti+psb_proc_id_)
      max_degree(degree) = proc_to_comm

      nerv = iaux(pnti+psb_n_elem_recv_)
      nesd = iaux(pnti+nerv+psb_n_elem_send_)
      ! print *, me, ": nerv", nerv, "nesd", nesd

      call psb_get_rank(prcid(proc_to_comm),ictxt,proc_to_comm)

      brvidx(proc_to_comm) = rcv_pt
      rvsz(proc_to_comm)   = nerv

      bsdidx(proc_to_comm) = snd_pt
      sdsz(proc_to_comm)   = nesd

      rcv_pt = rcv_pt + nerv
      snd_pt = snd_pt + nesd
      pnti   = pnti + nerv + nesd + 3
    end do

    ! source and destination nodes are the same
    allocate(src_and_dest(degree))
    src_and_dest = max_degree(1:degree)
    ! create graph comm of MPI rank's communication
    call MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, degree, src_and_dest, &
        MPI_UNWEIGHTED, degree, src_and_dest, MPI_UNWEIGHTED, MPI_INFO_NULL, &
        reorder, dist_graph_comm, ierr)
    deallocate(src_and_dest)

    allocate(cdesc%dist_graph_comm, stat=info)
    cdesc%dist_graph_comm = dist_graph_comm
    print *, "DIST_GRAPH_COMM CREATED"
    ! TESTING MPIX_INIT HERE, WILL BE MOVED ELSEWHERE
    ! THE SEND BUF, can't I just send iaux?
    ! FOR RECV BUF, best way to figure out how big?

    ! buf_size = size(iaux)
    ! allocate(rcv_buf(buf_size))
    ! rcv_buf = 0
    ! buf_size =  ! sum of
    ! buf_size = size(iaux) / 2
    ! call MPIX_Neighbor_alltoallv_init(iaux, buf_sizes, displs, MPI_INTEGER, &
    !     rcv_buf, buf_sizes, displs, MPI_INTEGER, graph_comm, MPI_INFO_NULL, &
    !     req, ierr)

    ! print *, "-----------------PRE START---------------"
    ! print *, me, ": sends", size(iaux), "to", src

    ! print *, me, ": buf_size", buf_size, "size(iaux)", size(iaux), "snd_buf", iaux
    ! call MPI_Start(req, ierr)
    ! call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)
    ! print *, me, ": buf_size", buf_size, "rcv_buf", rcv_buf

    ! print *, rank, ": iaux"
    ! print *, "====END WAIT===="

  end if
  ! ARTLESS: end of additions

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_crea_bnd_elem')
    goto 9999
  end if
  if (debug_level>0) write(debug_unit,*) me,'Done crea_bnd_elem'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psi_i_cnv_dsc


subroutine psi_i_inner_cnvs(x,hashmask,hashv,glb_lc)
  use psi_mod, psi_protect_name => psi_i_inner_cnvs

  integer(psb_ipk_), intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
  integer(psb_ipk_), intent(inout) :: x

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices
  ! are kept in a (hashed) list of ordered lists.
  ! Thus we first hash the index, then we do a binary search on the
  ! ordered sublist. The hashing is based on the low-order bits
  ! for a width of psb_hash_bits
  !

  key = x
  ih  = iand(key,hashmask)
  idx = hashv(ih)
  nh  = hashv(ih+1) - hashv(ih)
  if (nh > 0) then
    tmp = -1
    lb = idx
    ub = idx+nh-1
    do
      if (lb>ub) exit
      lm = (lb+ub)/2
      if (key == glb_lc(lm,1)) then
        tmp = lm
        exit
      else if (key<glb_lc(lm,1)) then
        ub = lm - 1
      else
        lb = lm + 1
      end if
    end do
  else
    tmp = -1
  end if
  if (tmp > 0) then
    x = glb_lc(tmp,2)
  else
    x = tmp
  end if
end subroutine psi_i_inner_cnvs

subroutine psi_i_inner_cnvs2(x,y,hashmask,hashv,glb_lc)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnvs2
  integer(psb_ipk_), intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
  integer(psb_ipk_), intent(in)  :: x
  integer(psb_ipk_), intent(out) :: y

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices
  ! are kept in a (hashed) list of ordered lists.
  ! Thus we first hash the index, then we do a binary search on the
  ! ordered sublist. The hashing is based on the low-order bits
  ! for a width of psb_hash_bits
  !

  key = x
  ih  = iand(key,hashmask)
  idx = hashv(ih)
  nh  = hashv(ih+1) - hashv(ih)
  if (nh > 0) then
    tmp = -1
    lb = idx
    ub = idx+nh-1
    do
      if (lb>ub) exit
      lm = (lb+ub)/2
      if (key == glb_lc(lm,1)) then
        tmp = lm
        exit
      else if (key<glb_lc(lm,1)) then
        ub = lm - 1
      else
        lb = lm + 1
      end if
    end do
  else
    tmp = -1
  end if
  if (tmp > 0) then
    y = glb_lc(tmp,2)
  else
    y = tmp
  end if
end subroutine psi_i_inner_cnvs2


subroutine psi_i_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnv1
  integer(psb_ipk_), intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
  logical, intent(in), optional    :: mask(:)
  integer(psb_ipk_), intent(inout) :: x(:)

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices
  ! are kept in a (hashed) list of ordered lists.
  ! Thus we first hash the index, then we do a binary search on the
  ! ordered sublist. The hashing is based on the low-order bits
  ! for a width of psb_hash_bits
  !
  if (present(mask)) then
    do i=1, n
      if (mask(i)) then
        key = x(i)
        ih  = iand(key,hashmask)
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih)
        if (nh > 0) then
          tmp = -1
          lb = idx
          ub = idx+nh-1
          do
            if (lb>ub) exit
            lm = (lb+ub)/2
            if (key == glb_lc(lm,1)) then
              tmp = lm
              exit
            else if (key<glb_lc(lm,1)) then
              ub = lm - 1
            else
              lb = lm + 1
            end if
          end do
        else
          tmp = -1
        end if
        if (tmp > 0) then
          x(i) = glb_lc(tmp,2)
        else
          x(i) = tmp
        end if
      end if
    end do
  else
    do i=1, n
      key = x(i)
      ih  = iand(key,hashmask)
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih)
      if (nh > 0) then
        tmp = -1
        lb = idx
        ub = idx+nh-1
        do
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm,1)) then
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else
        tmp = -1
      end if
      if (tmp > 0) then
        x(i) = glb_lc(tmp,2)
      else
        x(i) = tmp
      end if
    end do
  end if
end subroutine psi_i_inner_cnv1

subroutine psi_i_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnv2
  integer(psb_ipk_), intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
  logical, intent(in),optional  :: mask(:)
  integer(psb_ipk_), intent(in)  :: x(:)
  integer(psb_ipk_), intent(out) :: y(:)

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices
  ! are kept in a (hashed) list of ordered lists.
  ! Thus we first hash the index, then we do a binary search on the
  ! ordered sublist. The hashing is based on the low-order bits
  ! for a width of psb_hash_bits
  !
  if (present(mask)) then
    do i=1, n
      if (mask(i)) then
        key = x(i)
        ih  = iand(key,hashmask)
        if (ih > ubound(hashv,1) ) then
          write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
        end if
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih)
        if (nh > 0) then
          tmp = -1
          lb = idx
          ub = idx+nh-1
          do
            if (lb>ub) exit
            lm = (lb+ub)/2
            if (key == glb_lc(lm,1)) then
              tmp = lm
              exit
            else if (key<glb_lc(lm,1)) then
              ub = lm - 1
            else
              lb = lm + 1
            end if
          end do
        else
          tmp = -1
        end if
        if (tmp > 0) then
          y(i) = glb_lc(tmp,2)
        else
          y(i) = tmp
        end if
      end if
    end do
  else
    do i=1, n
      key = x(i)
      ih  = iand(key,hashmask)
      if (ih > ubound(hashv,1) ) then
        write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
      end if
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih)
      if (nh > 0) then
        tmp = -1
        lb = idx
        ub = idx+nh-1
        do
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm,1)) then
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else
        tmp = -1
      end if
      if (tmp > 0) then
        y(i) = glb_lc(tmp,2)
      else
        y(i) = tmp
      end if
    end do
  end if
end subroutine psi_i_inner_cnv2

subroutine psi_i_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
  use psi_mod, psi_protect_name =>  psi_i_bld_ovr_mst

  use psb_realloc_mod
  implicit none

  !     ....scalars parameters....
  integer(psb_ipk_), intent(in)               :: me, ovrlap_elem(:,:)
  integer(psb_ipk_), allocatable, intent(out) :: mst_idx(:)
  integer(psb_ipk_), intent(out)              :: info

  integer(psb_ipk_) :: i, j, proc, nov,isz, ip, err_act, idx
  character(len=20)  :: name

  name='psi_bld_ovr_mst'
  call psb_get_erraction(err_act)

  nov = size(ovrlap_elem,1)
  isz = 3*nov+1
  call psb_realloc(isz,mst_idx,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='reallocate')
    goto 9999
  end if
  mst_idx = -1
  j = 1
  do i=1, nov
    proc = ovrlap_elem(i,3)
    if (me /= proc) then
      idx = ovrlap_elem(i,1)
      mst_idx(j+0) = proc
      mst_idx(j+1) = 1
      mst_idx(j+2) = idx
      j = j + 3
    end if
  end do
  mst_idx(j) = -1

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psi_i_bld_ovr_mst
