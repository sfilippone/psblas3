!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
! File: psb_cdovr.f90
!
! Subroutine: psb_cdovr
!    This routine takes a matrix A with its descriptor, and builds the 
!    auxiliary descriptor corresponding to the number of overlap levels
!    specified on input. It really is just a size estimation/allocation
!    front end for <psb_cdovrbld>.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).       The input sparse matrix.
!    desc_a   - type(<psb_desc_type>).         The input communication descriptor.
!    norv     - integer.                       The number of overlap levels.
!    desc_ov  - type(<psb_desc_type>).         The auxiliary output communication 
!                                              descriptor.
!    info     - integer.                       Eventually returns an error code.
!
Subroutine psb_dcdovr(a,desc_a,novr,desc_ov,info, extype)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_error_mod
  use psb_penv_mod
  use psb_tools_mod, only : psb_cdcpy
  use psb_realloc_mod
  use psi_mod
  use mpi
  Implicit None

  !     .. Array Arguments ..
  integer, intent(in)                :: novr
  Type(psb_dspmat_type), Intent(in)  ::  a
  Type(psb_desc_type), Intent(in)    :: desc_a
  Type(psb_desc_type), Intent(inout) :: desc_ov
  integer, intent(out)               :: info
  integer, intent(in),optional       :: extype

  interface psb_icdasb
    subroutine psb_icdasb(desc_a,info,ext_hv)
      use psb_descriptor_type
      Type(psb_desc_type), intent(inout) :: desc_a
      integer, intent(out)               :: info
      logical, intent(in),optional       :: ext_hv
    end subroutine psb_icdasb
  end interface


  integer   icomm, err_act

  !     .. Local Scalars ..
  Integer ::  i, j, k, np, me,m,nnzero,&
       &  ictxt, lovr, lworks,lworkr, n_row,n_col, int_err(5),&
       &  index_dim,elem_dim, l_tmp_ovr_idx,l_tmp_halo, nztot,nhalo
  Integer :: counter,counter_h, counter_o, counter_e,idx,gidx,proc,n_elem_recv,&
       & n_elem_send,tot_recv,tot_elem,cntov_o,&
       & counter_t,n_elem,i_ovr,jj,proc_id,isz, mglob, glx, &
       & idxr, idxs, lx, iszr, iszs, nxch, nsnd, nrcv,lidx,irsv, extype_

  type(psb_dspmat_type) :: blk
  Integer, allocatable  :: tmp_halo(:),tmp_ovr_idx(:), orig_ovr(:)
  Integer,allocatable   :: halo(:),works(:),workr(:),t_halo_in(:),&
       & t_halo_out(:),temp(:),maskr(:)
  Integer,allocatable  :: brvindx(:),rvsz(:), bsdindx(:),sdsz(:)

  Logical,Parameter :: debug=.false.
  character(len=20) :: name, ch_err

  name='psb_cdovr'
  info  = 0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  Call psb_info(ictxt, me, np)

  If(debug) Write(0,*)'in psb_cdovr',novr

  if (present(extype)) then
    extype_ = extype
  else
    extype_ = psb_ovt_xhal_  
  endif
  m      = psb_cd_get_local_rows(desc_a)
  nnzero = Size(a%aspk)
  n_row  = psb_cd_get_local_rows(desc_a)
  n_col  = psb_cd_get_local_cols(desc_a)
  nhalo  = n_col-m

  If(debug) Write(0,*)'IN CDOVR1',novr ,m,nnzero,n_col
  if (novr<0) then
    info=10
    int_err(1)=1
    int_err(2)=novr
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  if (debug) write(0,*) 'Calling desccpy'
  call psb_cdcpy(desc_a,desc_ov,info)
  if (info /= 0) then
    info=4010
    ch_err='psb_cdcpy'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debug) write(0,*) 'From desccpy'
  if (novr==0) then 
    !
    ! Just copy the input.  
    !
    return
  endif

  call psb_get_mpicomm(ictxt,icomm )

  If(debug)then 
    Write(0,*)'BEGIN cdovr',me,nhalo
    call psb_barrier(ictxt)
  endif

  !
  ! Ok, since we are only estimating, do it as follows: 
  ! LOVR= (NNZ/NROW)*N_HALO*NOVR  This assumes that the local average 
  ! nonzeros per row is the same as the global. 
  !
  nztot = psb_sp_get_nnzeros(a)
  if (nztot>0) then 
    lovr   = ((nztot+m-1)/m)*nhalo*novr
    lworks = ((nztot+m-1)/m)*nhalo
    lworkr = ((nztot+m-1)/m)*nhalo
  else
    info=-1
    call psb_errpush(info,name)
    goto 9999
  endif
  If(debug)Write(0,*)'ovr_est done',me,novr,lovr

  index_dim = size(desc_a%halo_index)
  elem_dim  = size(desc_a%halo_index)

  l_tmp_ovr_idx = novr*(3*Max(2*index_dim,1)+1)
  l_tmp_halo    = novr*(3*Size(desc_a%halo_index))

!!$  write(0,*) 'Size of desc_ov ',    desc_ov%matrix_data(psb_desc_size_), &
!!$       & psb_desc_normal_,psb_desc_large_
  call psb_cd_set_bld(desc_ov,info)

  If(debug) then
    Write(0,*)'Start cdovrbld',me,lworks,lworkr
    call psb_barrier(ictxt)
  endif


  Allocate(brvindx(np+1),rvsz(np),sdsz(np),bsdindx(np+1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  Allocate(works(lworks),workr(lworkr),t_halo_in(l_tmp_halo),&
       & t_halo_out(l_tmp_halo), temp(lworkr),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if


  call psb_sp_all(blk,max(lworks,lworkr),info)
  if (info /= 0) then
    info=4010
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  blk%fida='COO'

  Allocate(orig_ovr(l_tmp_ovr_idx),tmp_ovr_idx(l_tmp_ovr_idx),&
       & tmp_halo(l_tmp_halo), halo(size(desc_a%halo_index)),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  halo(:)                = desc_a%halo_index(:)
  desc_ov%ovrlap_elem(:) = -1
  tmp_ovr_idx(:)         = -1
  orig_ovr(:)            = -1
  tmp_halo(:)            = -1
  counter_e              = 1
  tot_recv               = 0
  counter_h              = 1
  counter_o              = 1
  cntov_o                = 1
  ! Init overlap with desc_a%ovrlap (if any) 
  counter = 1
  Do While (desc_a%ovrlap_index(counter) /= -1)
    proc        = desc_a%ovrlap_index(counter+psb_proc_id_)
    n_elem_recv = desc_a%ovrlap_index(counter+psb_n_elem_recv_)
    n_elem_send = desc_a%ovrlap_index(counter+n_elem_recv+psb_n_elem_send_)

    Do j=0,n_elem_recv-1

      idx = desc_a%ovrlap_index(counter+psb_elem_recv_+j)
      If(idx > Size(desc_ov%loc_to_glob)) then 
        info=-3
        call psb_errpush(info,name)
        goto 9999
      endif

      gidx = desc_ov%loc_to_glob(idx)

      call psb_check_size((cntov_o+3),orig_ovr,info,pad=-1)
      if (info /= 0) then
        info=4010
        call psb_errpush(info,name,a_err='psb_check_size')
        goto 9999
      end if

      orig_ovr(cntov_o)=proc
      orig_ovr(cntov_o+1)=1
      orig_ovr(cntov_o+2)=gidx
      orig_ovr(cntov_o+3)=-1
      cntov_o=cntov_o+3
    end Do
    counter=counter+n_elem_recv+n_elem_send+2
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

    if (debug) write(0,*) me,'Running on overlap level ',i_ovr,' of ',novr

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

    desc_ov%matrix_data(psb_n_row_)=desc_ov%matrix_data(psb_n_col_)

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
      if (debug) write(0,*) me,' CDOVRBLD tot_recv:',proc,n_elem_recv,tot_recv
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
        If((counter+psb_elem_recv_+j)>Size(halo)) then 
          info=-2
          call psb_errpush(info,name)
          goto 9999
        end If

        idx = halo(counter+psb_elem_recv_+j)
        If(idx > Size(desc_ov%loc_to_glob)) then 
          info=-3
          call psb_errpush(info,name)
          goto 9999
        endif

        gidx = desc_ov%loc_to_glob(idx)

        call psb_check_size((counter_o+3),tmp_ovr_idx,info,pad=-1)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_check_size')
          goto 9999
        end if

        tmp_ovr_idx(counter_o)=proc
        tmp_ovr_idx(counter_o+1)=1
        tmp_ovr_idx(counter_o+2)=gidx
        tmp_ovr_idx(counter_o+3)=-1
        counter_o=counter_o+3
        call psb_check_size((counter_h+3),tmp_halo,info,pad=-1)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_check_size')
          goto 9999
        end if

        tmp_halo(counter_h)=proc
        tmp_halo(counter_h+1)=1
        tmp_halo(counter_h+2)=idx
        tmp_halo(counter_h+3)=-1

        counter_h=counter_h+3

      Enddo
      if (debug) write(0,*) me,'Checktmp_o_i Loop Mid1',tmp_ovr_idx(1:10)
      counter   = counter+n_elem_recv

      !
      ! add send elements in halo_index into ovrlap_index
      !
      Do j=0,n_elem_send-1

        idx = halo(counter+psb_elem_send_+j)
        gidx = desc_ov%loc_to_glob(idx)
        if (idx > psb_cd_get_local_rows(Desc_a)) &
             & write(0,*) me,i_ovr,'Out of local rows ',&
             & idx,psb_cd_get_local_rows(Desc_a)

        call psb_check_size((counter_o+3),tmp_ovr_idx,info,pad=-1)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_check_size')
          goto 9999
        end if

        tmp_ovr_idx(counter_o)=proc
        tmp_ovr_idx(counter_o+1)=1
        tmp_ovr_idx(counter_o+2)=gidx
        tmp_ovr_idx(counter_o+3)=-1
        counter_o=counter_o+3

        !
        ! Prepare to exchange the halo rows with the other proc.
        !
        If (i_ovr <= (novr)) Then
          n_elem = psb_sp_get_nnz_row(idx,a)

          call psb_check_size((idxs+tot_elem+n_elem),works,info)
          if (info /= 0) then
            info=4010
            call psb_errpush(info,name,a_err='psb_check_size')
            goto 9999
          end if

          If((n_elem) > size(blk%ia2)) Then
            isz = max((3*size(blk%ia2))/2,(n_elem))
            if (debug) write(0,*) me,'Realloc blk',isz
            call psb_sp_reall(blk,isz,info)
            if (info /= 0) then
              info=4010
              ch_err='psb_sp_reall'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          End If

          call psb_sp_getblk(idx,a,blk,info)
          if (info /= 0) then
            info=4010
            ch_err='psb_sp_getblk'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
!!$          write(0,*) me,'Iteration: ',j,i_ovr
          Do jj=1,n_elem
            works(idxs+tot_elem+jj)=desc_ov%loc_to_glob(blk%ia2(jj))
          End Do
          tot_elem=tot_elem+n_elem
        End If

      Enddo


      if (i_ovr <= novr) then 
        if (tot_elem > 1) then 
          call psb_msort(works(idxs+1:idxs+tot_elem))
          lx = works(idxs+1)
          i = 1 
          j = 1
          do 
            j = j + 1 
            if (j > tot_elem) exit
            if (works(idxs+j) /= lx) then 
              i = i + 1
              works(idxs+i) = works(idxs+j) 
              lx = works(idxs+i) 
            end if
          end do
          tot_elem = i
        endif
        if (debug) write(0,*) me,'Checktmp_o_i Loop Mid2',tmp_ovr_idx(1:10)
        sdsz(proc+1) = tot_elem
        idxs         = idxs + tot_elem
      end if
      counter   = counter+n_elem_send+3
      if (debug) write(0,*) me,'Checktmp_o_i Loop End',tmp_ovr_idx(1:10)
    Enddo
    if (debug) write(0,*)me,'End phase 1 CDOVRBLD', m, n_col, tot_recv

    if (i_ovr <= novr) then
      ! 
      ! Exchange data requests with everybody else: so far we have 
      ! accumulated RECV requests, we have an all-to-all to build
      ! matchings SENDs.
      !       
      call mpi_alltoall(sdsz,1,mpi_integer,rvsz,1,mpi_integer,icomm,info)
      if (info /= 0) then
        info=4010
        ch_err='mpi_alltoall'
        call psb_errpush(info,name,a_err=ch_err)
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
        if (info /= 0) then
          info=4010
          ch_err='psb_realloc'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        lworkr=max(iszr,1)
      end if

      call mpi_alltoallv(works,sdsz,bsdindx,mpi_integer,&
           & workr,rvsz,brvindx,mpi_integer,icomm,info)
      if (info /= 0) then
        info=4010
        ch_err='mpi_alltoallv'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      if (debug) write(0,*) 'ISZR :',iszr

      if (psb_is_large_desc(desc_ov)) then 
        call psb_check_size(iszr,maskr,info)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_check_size')
          goto 9999
        end if
        call psi_idx_cnv(iszr,workr,maskr,desc_ov,info)
        iszs = count(maskr(1:iszr)<=0)
        if (iszs > size(works)) call psb_realloc(iszs,works,info)
        j = 0
        do i=1,iszr
          if (maskr(i) < 0) then 
            j=j+1
            works(j) = workr(i)
          end if
        end do
        !
        ! fnd_owner on desc_a because we want the procs who
        ! owned the rows from the beginning!
        !
        call psi_fnd_owner(iszs,works,temp,desc_a,info)
        n_col=psb_cd_get_local_cols(desc_ov)

        do i=1,iszs
          idx = works(i)
          n_col   = psb_cd_get_local_cols(desc_ov)
          call psi_idx_ins_cnv(idx,lidx,desc_ov,info)
          if (psb_cd_get_local_cols(desc_ov) >  n_col ) then
            !
            ! This is a new index. Assigning a local index as
            ! we receive them guarantees that all indices for HALO(I)
            ! will be less than those for HALO(J) whenever I<J
            !
            proc_id = temp(i) 

            call psb_check_size((counter_t+3),t_halo_in,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_check_size')
              goto 9999
            end if

            t_halo_in(counter_t)=proc_id
            t_halo_in(counter_t+1)=1
            t_halo_in(counter_t+2)=lidx
            t_halo_in(counter_t+3)=-1
            counter_t=counter_t+3
            if (.false.) write(0,*) me,' CDOVRBLD: Added t_halo_in ',&
                 &proc_id,lidx,idx
          endif
        end Do
        n_col   = psb_cd_get_local_cols(desc_ov)

      else

        Do i=1,iszr
          idx=workr(i)
          if (idx <1) then 
            write(0,*) me,'Error in CDOVRBLD level',i_ovr,' : ',idx,i,iszr
          else  If (desc_ov%glob_to_loc(idx) < -np) Then
            !
            ! This is a new index. Assigning a local index as
            ! we receive them guarantees that all indices for HALO(I)
            ! will be less than those for HALO(J) whenever I<J
            !
            n_col=n_col+1
            proc_id=-desc_ov%glob_to_loc(idx)-np-1
            call psb_check_size(n_col,desc_ov%loc_to_glob,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_check_size')
              goto 9999
            end if

            desc_ov%glob_to_loc(idx)=n_col
            desc_ov%loc_to_glob(n_col)=idx

            call psb_check_size((counter_t+3),t_halo_in,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_check_size')
              goto 9999
            end if

            t_halo_in(counter_t)=proc_id
            t_halo_in(counter_t+1)=1
            t_halo_in(counter_t+2)=n_col
            t_halo_in(counter_t+3)=-1
            counter_t=counter_t+3
            if (debug) write(0,*) me,' CDOVRBLD: Added into t_halo_in from recv',&
                 &proc_id,n_col,idx
          else if (desc_ov%glob_to_loc(idx) < 0) Then
            if (debug) write(0,*) me,'Wrong input to cdovrbld??',&
                 &idx,desc_ov%glob_to_loc(idx)
          End If
        End Do
        desc_ov%matrix_data(psb_n_col_)=n_col

      end if

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


      if (debug) write(0,*) me,'Checktmp_o_i 1',tmp_ovr_idx(1:10)
      if (debug) write(0,*) me,'Calling Crea_Halo'

      call psi_crea_index(desc_ov,t_halo_in,t_halo_out,.false.,&
           & nxch,nsnd,nrcv,info)

      if (debug) then 
        write(0,*) me,'Done Crea_Index'
        call psb_barrier(ictxt)
      end if
      if (debug) write(0,*) me,'Checktmp_o_i 2',tmp_ovr_idx(1:10)
      if (debug) write(0,*) me,'Done Crea_Halo'
      call psb_transfer(t_halo_out,halo,info)
      !
      ! At this point we have built the halo necessary for I_OVR+1.
      !
    End If
    if (debug) write(0,*) me,'Checktmp_o_i ',tmp_ovr_idx(1:10)

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
    desc_ov%matrix_data(psb_n_row_) = desc_a%matrix_data(psb_n_row_)
    call psb_transfer(orig_ovr,desc_ov%ovrlap_index,info)
    call psb_check_size((counter_h+counter_t+1),tmp_halo,info,pad=-1)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_check_size')
      goto 9999
    end if
    tmp_halo(counter_h:counter_h+counter_t-1) = t_halo_in(1:counter_t)
    counter_h            = counter_h+counter_t-1
    tmp_halo(counter_h:) = -1
    call psb_transfer(tmp_halo,desc_ov%halo_index,info)
    deallocate(tmp_ovr_idx,stat=info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='deallocate')
      goto 9999
    end if
    
  case(psb_ovt_asov_)
    !
    ! Build an overlapped descriptor for Additive Schwarz 
    !  with overlap enlargement; we need the overlap, 
    !  the (new) halo and the mapping into the new index space. 
    ! Here we need: 1. (orig_ovr + t_ovr_idx -> ovrlap
    !               2. (tmp_halo) -> ext
    !               3. (t_halo_in) -> halo 
    !               4. n_row(ov)  current.
    !               5. n_col(ov)  current. 
    ! 
    call psb_check_size((cntov_o+counter_o+1),orig_ovr,info,pad=-1)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_check_size')
      goto 9999
    end if
    orig_ovr(cntov_o:cntov_o+counter_o-1) = tmp_ovr_idx(1:counter_o)
    cntov_o = cntov_o+counter_o-1
    orig_ovr(cntov_o:) = -1
    call psb_transfer(orig_ovr,desc_ov%ovrlap_index,info)
    deallocate(tmp_ovr_idx,stat=info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='deallocate')
      goto 9999
    end if
    tmp_halo(counter_h:) = -1
    call psb_transfer(tmp_halo,desc_ov%ext_index,info)
    call psb_transfer(t_halo_in,desc_ov%halo_index,info)

  case default
    call psb_errpush(30,name,i_err=(/5,extype_,0,0,0/))
    goto 9999
  end select

  !
  ! At this point we have gathered all the indices in the halo at
  ! N levels of overlap. Just call icdasb forcing to use 
  ! the halo_index provided. This is the same routine as gets 
  ! called inside CDASB.
  !

  if (debug) then
    write(0,*) 'psb_cdovrbld: converting indexes'
    call psb_barrier(ictxt)
  end if
  call psb_icdasb(desc_ov,info,ext_hv=.true.)

  if (debug) then 
    write(0,*) me,'Done CDASB'
    call psb_barrier(ictxt)
  end if

  if (info == 0) call psb_sp_free(blk,info)
  if (info /= 0) then
    ch_err='sp_free'
    call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  Return

End Subroutine psb_dcdovr

