!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_zcdbldext.f90
!
! Subroutine: psb_zcdbldext
!    This routine takes a matrix A with its descriptor, and builds the 
!    auxiliary descriptor corresponding to the number of overlap levels
!    specified on input. 
! 
! Arguments: 
!    a        - type(psb_zspmat_type).       The input sparse matrix.
!    desc_a   - type(psb_desc_type).         The input communication descriptor.
!    novr     - integer.                       The number of overlap levels.
!    desc_ov  - type(psb_desc_type).         The auxiliary output communication 
!                                              descriptor.
!    info     - integer.                       Return code.
!    extype   - integer.                       Choice of type of overlap:
!                               psb_ovt_xhal_: build a descriptor with an extended 
!                                              stencil, i.e. enlarge the existing 
!                                              halo by novr additional layers.
!                               psb_ovt_asov_: build a descriptor suitable 
!                                              for Additive Schwarz preconditioner.
!                                              This last choice implies that: 
!                                              a. The novr halo layers are added to
!                                                 the overlap;
!                                              b. The novr layers are also copied to 
!                                                 the ext_ structure to provide 
!                                                 the mapping between the base 
!                                                 descriptor and the overlapped one.
!                                              c. The (novr+1)-th layer becomes the
!                                                 new halo.
!
Subroutine psb_zcdbldext(a,desc_a,novr,desc_ov,info, extype)

  use psb_base_mod, psb_protect_name => psb_zcdbldext
  use psi_mod

#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif

  !     .. Array Arguments ..
  integer(psb_ipk_), intent(in)                      :: novr
  Type(psb_zspmat_type), Intent(in)       ::  a
  Type(psb_desc_type), Intent(inout), target :: desc_a
  Type(psb_desc_type), Intent(out)        :: desc_ov
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_), intent(in),optional            :: extype

  !     .. Local Scalars ..
  integer(psb_ipk_) ::  i, j, err_act,m,&
       &  lovr, lworks,lworkr, n_row,n_col, n_col_prev, &
       &  index_dim,elem_dim, l_tmp_ovr_idx,l_tmp_halo, nztot,nhalo
  integer(psb_ipk_) :: counter,counter_h, counter_o, counter_e,idx,gidx,proc,n_elem_recv,&
       & n_elem_send,tot_recv,tot_elem,cntov_o,&
       & counter_t,n_elem,i_ovr,jj,proc_id,isz, &
       & idxr, idxs, iszr, iszs, nxch, nsnd, nrcv,lidx, extype_
  integer(psb_mpik_) :: icomm, ictxt, me, np, minfo

  integer(psb_ipk_), allocatable :: irow(:), icol(:)
  integer(psb_ipk_), allocatable :: tmp_halo(:),tmp_ovr_idx(:), orig_ovr(:)
  integer(psb_ipk_), allocatable :: halo(:),ovrlap(:),works(:),workr(:),&
       & t_halo_in(:), t_halo_out(:),temp(:),maskr(:)
  integer(psb_mpik_),allocatable :: brvindx(:),rvsz(:), bsdindx(:),sdsz(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name, ch_err

  name='psb_zcdbldext'
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
  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  Call psb_info(ictxt, me, np)

  If (debug_level >= psb_debug_outer_) &
       & Write(debug_unit,*) me,' ',trim(name),&
       & ': start',novr

  if (present(extype)) then
    extype_ = extype
  else
    extype_ = psb_ovt_xhal_  
  endif
  m      = desc_a%get_local_rows()
  n_row  = desc_a%get_local_rows()
  n_col  = desc_a%get_local_cols()
  nhalo  = n_col-m

  if (novr<0) then
    info=psb_err_iarg_neg_
    ierr(1)=1; ierr(2)=novr
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif

  select case(extype_) 
  case(psb_ovt_xhal_,psb_ovt_asov_)
  case default
    ierr(1)=5; ierr(2)=extype_
    call psb_errpush(psb_err_input_value_invalid_i_,&
         & name,i_err=ierr)
    goto 9999
  end select

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ': Calling desccpy'

  call psb_cdcpy(desc_a,desc_ov,info)

  if (psb_errstatus_fatal()) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_cdcpy')
    goto 9999
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ': From desccpy'

  if ((novr == 0).or.(np==1)) then 
    !
    ! Just copy the input.  
    ! Should we return also when is_repl() ? 
    !
    return
  endif


  if (extype_ == psb_ovt_asov_) then 
    ! Need to switch to a format that can support overlap,
    ! so far: LIST or HASH. This will also reinitialize properly
    ! the inex map contents. Encapsulate choice
    ! in a separate method.
    call psb_cd_switch_ovl_indxmap(desc_ov,info) 
  end if
  if (info == 0) call psb_cd_set_ovl_bld(desc_ov,info)
  if (info /= 0) goto 9999

  If (debug_level >= psb_debug_outer_)then 
    Write(debug_unit,*) me,' ',trim(name),&
         & ': BEGIN ',nhalo, desc_ov%indxmap%get_state()
    call psb_barrier(ictxt)
  endif
  !
  ! Ok, since we are only estimating, do it as follows: 
  ! LOVR= (NNZ/NROW)*N_HALO*NOVR  This assumes that the local average 
  ! nonzeros per row is the same as the global. 
  !
  ! Allow for empty matrices.
  nztot = max(ione,a%get_nzeros())
  if (nztot>0) then 
    lovr   = ((nztot+m-1)/m)*nhalo*novr
    lworks = ((nztot+m-1)/m)*nhalo
    lworkr = ((nztot+m-1)/m)*nhalo
  endif
  If (debug_level >= psb_debug_outer_)&
       & Write(debug_unit,*) me,' ',trim(name),':ovr_est done',novr,lovr

  index_dim = max(desc_a%v_halo_index%get_nrows(),1_psb_ipk_)
  elem_dim  = index_dim

  l_tmp_ovr_idx = novr*(3*Max(2*index_dim,1)+1)
  l_tmp_halo    = novr*(3*index_dim)

  desc_ov%base_desc => desc_a

  If (debug_level >= psb_debug_outer_) then
    Write(debug_unit,*) me,' ',trim(name),':Start',&
         & lworks,lworkr, desc_ov%indxmap%get_state()
    call psb_barrier(ictxt)
  endif


  Allocate(brvindx(np+1),rvsz(np),sdsz(np),bsdindx(np+1),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  Allocate(works(lworks),workr(lworkr),t_halo_in(l_tmp_halo),&
       & t_halo_out(l_tmp_halo), temp(lworkr),stat=info)
  if (info == psb_success_) allocate(orig_ovr(l_tmp_ovr_idx),&
       & tmp_ovr_idx(l_tmp_ovr_idx), &
       & tmp_halo(l_tmp_halo),stat=info)

  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999      
  end if

  halo           = desc_a%v_halo_index%get_vect()
  if (.not.allocated(halo)) halo = (/ -ione /)
  ovrlap         = desc_a%v_ovrlap_index%get_vect()
  if (.not.allocated(ovrlap)) ovrlap = (/ -ione /)    

  tmp_ovr_idx(:) = -1
  orig_ovr(:)    = -1
  tmp_halo(:)    = -1
  counter_e      =  1
  tot_recv       =  0
  counter_t      =  1
  counter_h      =  1
  counter_o      =  1
  cntov_o        =  1
  ! Init overlap with desc_a%ovrlap (if any) 
  counter = 1
  Do While (ovrlap(counter) /= -1)
    proc        = ovrlap(counter+psb_proc_id_)
    n_elem_recv = ovrlap(counter+psb_n_elem_recv_)
    n_elem_send = ovrlap(counter+n_elem_recv+psb_n_elem_send_)

    Do j=0,n_elem_recv-1

      idx = ovrlap(counter+psb_elem_recv_+j)
      call desc_ov%indxmap%l2g(idx,gidx,info) 
      If (gidx < 0) then 
        info=-3
        call psb_errpush(info,name)
        goto 9999
      endif
      call psb_ensure_size((cntov_o+3),orig_ovr,info,pad=-ione)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_ensure_size')
        goto 9999
      end if
      orig_ovr(cntov_o)=proc
      orig_ovr(cntov_o+1)=1
      orig_ovr(cntov_o+2)=gidx
      orig_ovr(cntov_o+3)=-1
      cntov_o=cntov_o+3
    end Do
    counter=counter+n_elem_recv+n_elem_send+3
  end Do




  !
  ! A picture is in order to understand what goes on here.
  ! I is the internal part; H is halo, R row, C column. The final
  ! matrix with N levels of overlap looks like this
  !
  !             I   | Hc1 |  0  |  0  |
  !          -------|-----|-----|-----|
  !            Hr1  | Hd1 | Hc2 |  0  |
  !          -------|-----|-----|-----|
  !             0   | Hr2 | Hd2 | Hc2 |
  !          -------|-----|-----|-----|
  !             0   |  0  | Hr3 | Hd3 | Hc3
  !
  ! At the start we already have I and Hc1, so we know the row
  ! indices that will make up Hr1, and also who owns them. As we
  ! actually get those rows, we receive the column indices in Hc2;
  ! these define the row indices for Hr2, and so on. When we have
  ! reached the desired level HrN. 
  !
  !
  Do i_ovr = 1, novr

    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ':Running on overlap level ',i_ovr,' of ',novr

    !
    ! At this point, halo contains a valid halo corresponding to the
    ! matrix enlarged with the elements in the frontier for I_OVR-1.
    ! At the start, this is just the halo for A; the rows for indices in
    ! the first halo will contain column indices defining the second halo
    ! level and so on.
    !
    bsdindx(:) = 0
    sdsz(:)    = 0
    brvindx(:) = 0
    rvsz(:)    = 0
    idxr       = 0
    idxs       = 0
    counter    = 1
    counter_t  = 1

    n_col_prev = desc_ov%get_local_cols() 

    Do While (halo(counter) /= -1)
      tot_elem=0
      proc=halo(counter+psb_proc_id_)
      n_elem_recv=halo(counter+psb_n_elem_recv_)
      n_elem_send=halo(counter+n_elem_recv+psb_n_elem_send_)
      If ((counter+n_elem_recv+n_elem_send) > Size(halo)) then
        info = -1
        call psb_errpush(info,name)
        goto 9999
      end If
      tot_recv=tot_recv+n_elem_recv
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           &  ': tot_recv:',proc,n_elem_recv,tot_recv
      !
      !
      ! The format of the halo vector exists in two forms: 1. Temporary 
      ! 2. Assembled. In this loop we are using the (assembled) halo_in and 
      ! copying it into (temporary) halo_out; this is because tmp_halo will
      ! be enlarged with the new column indices received, and will reassemble
      ! everything for the next iteration.
      ! 
      !
      ! add recv elements in halo_index into ovrlap_index
      !
      Do j=0,n_elem_recv-1
        If ((counter+psb_elem_recv_+j)>Size(halo)) then 
          info=-2
          call psb_errpush(info,name)
          goto 9999
        end If

        idx = halo(counter+psb_elem_recv_+j)
        call desc_ov%l2g(idx,gidx,info) 
        If (gidx < 0) then 
          info=-3
          call psb_errpush(info,name)
          goto 9999
        endif
        call psb_ensure_size((counter_o+3),tmp_ovr_idx,info,pad=-ione)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if

        tmp_ovr_idx(counter_o)   = proc
        tmp_ovr_idx(counter_o+1) = 1
        tmp_ovr_idx(counter_o+2) = gidx
        tmp_ovr_idx(counter_o+3) = -1
        counter_o=counter_o+3
        call psb_ensure_size((counter_h+3),tmp_halo,info,pad=-ione)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if

        tmp_halo(counter_h)   = proc
        tmp_halo(counter_h+1) = 1
        tmp_halo(counter_h+2) = idx
        tmp_halo(counter_h+3) = -1

        counter_h=counter_h+3

      Enddo
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ':Checktmp_o_i Loop Mid1',tmp_ovr_idx(1:10)
      counter   = counter+n_elem_recv

      !
      ! add send elements in halo_index into ovrlap_index
      !
      Do j=0,n_elem_send-1

        idx = halo(counter+psb_elem_send_+j)
        call desc_ov%l2g(idx,gidx,info) 
        If (gidx < 0) then 
          info=-3
          call psb_errpush(info,name)
          goto 9999
        endif
        call psb_ensure_size((counter_o+3),tmp_ovr_idx,info,pad=-ione)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if

        tmp_ovr_idx(counter_o)   = proc
        tmp_ovr_idx(counter_o+1) = 1
        tmp_ovr_idx(counter_o+2) = gidx
        tmp_ovr_idx(counter_o+3) = -1
        counter_o=counter_o+3

        !
        ! Prepare to exchange the halo rows with the other proc.
        !
        If (i_ovr <= (novr)) Then
          call a%csget(idx,idx,n_elem,irow,icol,info)
          if (info /= psb_success_) then
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='csget')
            goto 9999
          end if

          call psb_ensure_size((idxs+tot_elem+n_elem),works,info)
          if (info /= psb_success_) then
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='psb_ensure_size')
            goto 9999
          end if
          call desc_ov%l2g(icol(1:n_elem),&
               & works(idxs+tot_elem+1:idxs+tot_elem+n_elem),&
               & info) 

          tot_elem=tot_elem+n_elem
        End If

      Enddo


      if (i_ovr <= novr) then 
        if (tot_elem > 1) then 
          call psb_msort_unique(works(idxs+1:idxs+tot_elem),i)
          tot_elem=i
        endif

        sdsz(proc+1) = tot_elem
        idxs         = idxs + tot_elem
      end if
      counter   = counter+n_elem_send+3

    Enddo
    if (debug_level >= psb_debug_outer_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ':End phase 1', m, n_col, tot_recv

    if (i_ovr <= novr) then
      ! 
      ! Exchange data requests with everybody else: so far we have 
      ! accumulated RECV requests, we have an all-to-all to build
      ! matchings SENDs.
      !       
      call mpi_alltoall(sdsz,1,psb_mpi_def_integer,rvsz,1, &
           & psb_mpi_def_integer,icomm,minfo)
      if (minfo /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='mpi_alltoall')
        goto 9999
      end if
      idxs = 0
      idxr = 0
      counter = 1
      Do 
        proc=halo(counter)
        if (proc == -1) exit
        n_elem_recv = halo(counter+psb_n_elem_recv_)
        counter     = counter+n_elem_recv
        n_elem_send = halo(counter+psb_n_elem_send_)

        bsdindx(proc+1) = idxs
        idxs = idxs + sdsz(proc+1)
        brvindx(proc+1) = idxr
        idxr = idxr + rvsz(proc+1)
        counter   = counter+n_elem_send+3
      Enddo

      iszr=sum(rvsz)
      if (max(iszr,1) > lworkr) then 
        call psb_realloc(max(iszr,1),workr,info)

        if (psb_errstatus_fatal()) then
          info=psb_err_alloc_dealloc_
          call psb_errpush(info,name)
          goto 9999
        end if
        lworkr = max(iszr,1)
      end if

      call mpi_alltoallv(works,sdsz,bsdindx,psb_mpi_ipk_integer,&
           & workr,rvsz,brvindx,psb_mpi_ipk_integer,icomm,minfo)
      if (minfo /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='mpi_alltoallv')
        goto 9999
      end if

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),': ISZR :',iszr

      call psb_ensure_size(iszr,maskr,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_ensure_size')
        goto 9999
      end if
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ': going for first idx_cnv', desc_ov%indxmap%get_state()

      call desc_ov%indxmap%g2l(workr(1:iszr),maskr(1:iszr),info)
      iszs = count(maskr(1:iszr)<=0)
      if (iszs > size(works)) call psb_realloc(iszs,works,info)
      j = 0
      do i=1,iszr
        if (maskr(i) < 0) then 
          j=j+1
          works(j) = workr(i)
        end if
      end do
      ! Eliminate duplicates from request
      call psb_msort_unique(works(1:j),iszs)

      !
      ! fnd_owner on desc_a because we want the procs who
      ! owned the rows from the beginning!
      !
      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ': going for fnd_owner', desc_ov%indxmap%get_state()
      call desc_a%fnd_owner(works(1:iszs),temp,info)
      n_col = desc_ov%get_local_cols()

      if (debug_level >= psb_debug_outer_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ': Done fnd_owner', desc_ov%indxmap%get_state()

      do i=1,iszs
        idx = works(i)
        n_col   = desc_ov%get_local_cols()
        call desc_ov%indxmap%g2l_ins(idx,lidx,info)
        if (desc_ov%get_local_cols() >  n_col ) then
          !
          ! This is a new index. Assigning a local index as
          ! we receive them guarantees that all indices for HALO(I)
          ! will be less than those for HALO(J) whenever I<J
          !
          proc_id = temp(i) 

          call psb_ensure_size((counter_t+3),t_halo_in,info,pad=-ione)
          if (info /= psb_success_) then
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='psb_ensure_size')
            goto 9999
          end if

          t_halo_in(counter_t)   = proc_id
          t_halo_in(counter_t+1) = 1
          t_halo_in(counter_t+2) = lidx
          t_halo_in(counter_t+3) = -1
          counter_t              = counter_t+3
        endif
      end Do
      n_col   = desc_ov%get_local_cols()

    end if

    !
    ! Ok, now we have a temporary halo with all the info for the
    ! next round. If we need to keep going, convert the halo format
    ! from temporary to final, so that we can work out the next iteration.
    ! This uses one of the convert_comm internals, i.e. we are doing
    ! the equivalent of a partial call to convert_comm
    !
    t_halo_in(counter_t)=-1

    If (i_ovr < (novr)) Then


      if (debug_level >= psb_debug_outer_) then
        write(debug_unit,*) me,' ',trim(name),':Checktmp_o_i 1',tmp_ovr_idx(1:10)      
        write(debug_unit,*) me,' ',trim(name),':Calling Crea_index'
      end if

      call psi_crea_index(desc_ov,t_halo_in,t_halo_out,.false.,&
           & nxch,nsnd,nrcv,info)

      if (debug_level >= psb_debug_outer_) then 
        write(debug_unit,*) me,' ',trim(name),':Done Crea_Index'
        call psb_barrier(ictxt)
      end if
      call psb_move_alloc(t_halo_out,halo,info)
      !
      ! At this point we have built the halo necessary for I_OVR+1.
      !
    End If
  End Do

  select case(extype_) 
  case(psb_ovt_xhal_)
    !
    ! Build an extended-stencil halo, but no overlap enlargement. 
    ! Here we need: 1. orig_ovr -> ovrlap
    !               2. (tmp_halo + t_halo_in) -> halo
    !               3. (t_ovr_idx) -> /dev/null 
    !               4. n_row(ov) = n_row(a)
    !               5. n_col(ov)  current. 
    !
    call psb_move_alloc(orig_ovr,desc_ov%ovrlap_index,info)
    call psb_ensure_size((counter_h+counter_t+1),tmp_halo,info,pad=-ione)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_ensure_size')
      goto 9999
    end if
    tmp_halo(counter_h:counter_h+counter_t-1) = t_halo_in(1:counter_t)
    counter_h            = counter_h+counter_t-1
    tmp_halo(counter_h:) = -1
    call psb_move_alloc(tmp_halo,desc_ov%halo_index,info)
    deallocate(tmp_ovr_idx,stat=info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='deallocate')
      goto 9999
    end if

  case(psb_ovt_asov_)
    !
    ! Build an overlapped descriptor for Additive Schwarz 
    !  with overlap enlargement; we need the overlap, 
    !  the (new) halo and the mapping into the new index space. 
    ! Here we need: 1. (orig_ovr + t_ovr_idx) -> ovrlap
    !               2. (tmp_halo) -> ext
    !               3. (t_halo_in) -> halo 
    !               4. n_row(ov)  current.
    !               5. n_col(ov)  current. 
    ! 
    call desc_ov%indxmap%set_lr(n_col_prev)
    call psb_ensure_size((cntov_o+counter_o+1),orig_ovr,info,pad=-ione)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_ensure_size')
      goto 9999
    end if
    orig_ovr(cntov_o:cntov_o+counter_o-1) = tmp_ovr_idx(1:counter_o)
    cntov_o = cntov_o+counter_o-1
    orig_ovr(cntov_o:) = -1
    call psb_move_alloc(orig_ovr,desc_ov%ovrlap_index,info)
    deallocate(tmp_ovr_idx,stat=info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='deallocate')
      goto 9999
    end if
    tmp_halo(counter_h:) = -1
    call psb_move_alloc(tmp_halo,desc_ov%ext_index,info)
    call psb_move_alloc(t_halo_in,desc_ov%halo_index,info)
  case default
    ierr(1)=5; ierr(2)=extype_
    call psb_errpush(psb_err_input_value_invalid_i_,name,i_err=ierr)
    goto 9999
  end select

  !
  ! At this point we have gathered all the indices in the halo at
  ! N levels of overlap. Just call icdasb forcing to use 
  ! the halo_index provided. This is the same routine as gets 
  ! called inside CDASB.
  !

  if (debug_level >= psb_debug_outer_) then
    write(debug_unit,*) me,' ',trim(name),': converting indexes'
    call psb_barrier(ictxt)
  end if

  call psb_icdasb(desc_ov,info,ext_hv=.true.)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='icdasdb')
    goto 9999
  end if

  call psb_cd_set_ovl_asb(desc_ov,info)

  if (info == psb_success_)  then 
    if (allocated(irow)) deallocate(irow,stat=info)
    if ((info == psb_success_).and.allocated(icol)) &
         & deallocate(icol,stat=info)
    if (info /= psb_success_) then
      ierr(1) = info
      call psb_errpush(psb_err_from_subroutine_ai_,name, &
           & a_err='deallocate',i_err=ierr)
      goto 9999
    end if
  end if

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

End Subroutine psb_zcdbldext

