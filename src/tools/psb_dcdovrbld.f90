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
! File: psb_dcdovrbld.f90
!
! Subroutine: psb_dcdovrbld
!     This routine takes a matrix A with its descriptor, and builds the  
!     auxiliary descriptor corresponding to the number of overlap levels 
!     specified on input. This is the actual worker horse.....           
!     Note that n_ovr > 0 thanks to the caller routine.                  
! 
! Parameters: 
!    n_ovr         - integer.                The number of overlap levels               
!    desc_p        - type(<psb_desc_type>).  The communication descriptor 
!                                            for the preconditioner.
!    desc_a        - type(<psb_desc_type>).  The communication descriptor.
!    a             - type(<psb_dspmat_type). The matrix upon which the preconditioner
!                                            will be built.
!    l_tmp_halo    - integer.                Input estimate for allocation sizes.
!    l_tmp_ovr_idx - integer.                Input estimate for allocation sizes.
!    lworkr        - integer.                Input estimate for allocation sizes.
!    lworks        - integer.                Input estimate for allocation sizes.
!    info          - integer.                Possibly return an error code
Subroutine psb_dcdovrbld(n_ovr,desc_p,desc_a,a,&
     &       l_tmp_halo,l_tmp_ovr_idx,lworks,lworkr,info)
  use psb_descriptor_type
  use psb_serial_mod
  Use psi_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_const_mod
  use psb_penv_mod
  use psb_tools_mod
  use mpi
  Implicit None

  type(psb_dspmat_type),intent(in)  :: a
  type(psb_desc_type),intent(in)    :: desc_a
  type(psb_desc_type),intent(inout) :: desc_p
  integer,intent(in)                :: n_ovr
  ! Input estimates for allocation sizes. Could we do without these two???
  Integer, Intent(in)     :: l_tmp_halo,l_tmp_ovr_idx
  Integer, Intent(inout)  :: lworks, lworkr
  integer, intent(out)    :: info

  type(psb_dspmat_type)  :: blk
  Integer, allocatable   :: tmp_halo(:),tmp_ovr_idx(:)
  Integer :: counter,counter_h, counter_o, counter_e,j,idx,gidx,proc,n_elem_recv,&
       & n_elem_send,tot_recv,tot_elem,n_col,m,ictxt,np,me,dl_lda,lwork,&
       & counter_t,n_elem,i_ovr,jj,i,proc_id,isz, mglob, glx,n_row, &
       & idxr, idxs, lx, iszr, iszs, err_act, icomm, nxch, nsnd, nrcv,lidx
  Integer,allocatable   :: halo(:),works(:),workr(:),t_halo_in(:),&
       & t_halo_out(:),temp(:),maskr(:)
  Integer,allocatable  :: brvindx(:),rvsz(:), bsdindx(:),sdsz(:)

  Logical,Parameter :: debug=.false.
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,t7, tl, tch
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_cdovrbld'
  call psb_erractionsave(err_act)

  If(debug) Write(0,*)'cdovrbld begin'
  ictxt = psb_cd_get_context(desc_a)

  Call psb_info(ictxt,me,np)

  Allocate(brvindx(np+1),rvsz(np),sdsz(np),bsdindx(np+1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  tl = 0.0
  tch = 0.0
  t4 = 0.0
  call psb_get_mpicomm(ictxt,icomm )

  mglob = psb_cd_get_global_rows(desc_a)
  m     = psb_cd_get_local_rows(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)
  if (debug) write(0,*) me,' On entry to CDOVRBLD n_col:',n_col

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

  Allocate(tmp_ovr_idx(l_tmp_ovr_idx),tmp_halo(l_tmp_halo),&
       & halo(size(desc_a%halo_index)),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  halo(:)               = desc_a%halo_index(:)
  desc_p%ovrlap_elem(:) = -1
  tmp_ovr_idx(:)        = -1
  tmp_halo(:)           = -1
  counter_e             = 1
  tot_recv              = 0
  counter_h             = 1
  counter_o             = 1

  ! Init overlap with desc_a%ovrlap (if any) 
  counter = 1
  Do While (desc_a%ovrlap_index(counter) /= -1)
    proc        = desc_a%ovrlap_index(counter+psb_proc_id_)
    n_elem_recv = desc_a%ovrlap_index(counter+psb_n_elem_recv_)
    n_elem_send = desc_a%ovrlap_index(counter+n_elem_recv+psb_n_elem_send_)

    Do j=0,n_elem_recv-1

      idx = desc_a%ovrlap_index(counter+psb_elem_recv_+j)
      If(idx > Size(desc_p%loc_to_glob)) then 
        info=-3
        call psb_errpush(info,name)
        goto 9999
      endif

      gidx = desc_p%loc_to_glob(idx)

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
  ! reached the desired level HrN, we may ignore HcN.
  !
  !
  Do i_ovr=1,n_ovr

    if (debug) write(0,*) me,'Running on overlap level ',i_ovr,' of ',n_ovr

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
    counter      = 1
    counter_t    = 1

    t1 = mpi_wtime()
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

      t3 = mpi_wtime()
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
        If(idx > Size(desc_p%loc_to_glob)) then 
          info=-3
          call psb_errpush(info,name)
          goto 9999
        endif

        gidx = desc_p%loc_to_glob(idx)

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
        if (.not.psb_is_large_desc(desc_p)) then 
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
        end if

      Enddo
      if (debug) write(0,*) me,'Checktmp_o_i Loop Mid1',tmp_ovr_idx(1:10)
      counter   = counter+n_elem_recv

      !
      ! add send elements in halo_index into ovrlap_index
      !
      Do j=0,n_elem_send-1

        idx = halo(counter+psb_elem_send_+j)
        gidx = desc_p%loc_to_glob(idx)
        if (idx > psb_cd_get_local_rows(Desc_a)) &
             & write(0,*) me,i_ovr,'Out of local rows ',idx,psb_cd_get_local_rows(Desc_a)
        
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
        If (i_ovr < (n_ovr)) Then
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

          Do jj=1,n_elem
            works(idxs+tot_elem+jj)=desc_p%loc_to_glob(blk%ia2(jj))
          End Do
          tot_elem=tot_elem+n_elem
        End If

      Enddo

      t4 = t4 + mpi_wtime() -t3

      if (i_ovr < n_ovr) then 
        if (tot_elem > 1) then 
          call imsr(tot_elem,works(idxs+1))
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

    if (i_ovr < n_ovr) then
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
      if (.false.) then 
        open(70+me,position='append')
        write(70+me,*) ' Workr on iteration ',i_ovr
        write(70+me,'(8(i7,1x))') workr(1:iszr)
        call flush(70+me)
        close(70+me)
      end if

      if (debug) write(0,*) 'ISZR :',iszr

      if (psb_is_large_desc(desc_a)) then 
        call psb_check_size(iszr,maskr,info)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_check_size')
          goto 9999
        end if
        call psi_idx_cnv(iszr,workr,maskr,desc_p,info)
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
        n_col=psb_cd_get_local_cols(desc_p)

        do i=1,iszs
          idx = works(i)
          n_col   = psb_cd_get_local_cols(desc_p)
          call psi_idx_ins_cnv(idx,lidx,desc_p,info)
          if (psb_cd_get_local_cols(desc_p) >  n_col ) then
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
            counter_t=counter_t+3
            if (.false.) write(0,*) me,' CDOVRBLD: Added t_halo_in ',&
                 &proc_id,lidx,idx
          endif
        end Do

      else

        Do i=1,iszr
          idx=workr(i)
          if (idx <1) then 
            write(0,*) me,'Error in CDOVRBLD level',i_ovr,' : ',idx,i,iszr
          else  If (desc_p%glob_to_loc(idx) < -np) Then
            !
            ! This is a new index. Assigning a local index as
            ! we receive them guarantees that all indices for HALO(I)
            ! will be less than those for HALO(J) whenever I<J
            !
            n_col=n_col+1
            proc_id=-desc_p%glob_to_loc(idx)-np-1
            call psb_check_size(n_col,desc_p%loc_to_glob,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_check_size')
              goto 9999
            end if

            desc_p%glob_to_loc(idx)=n_col
            desc_p%loc_to_glob(n_col)=idx

            call psb_check_size((counter_t+3),t_halo_in,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_check_size')
              goto 9999
            end if
            
            t_halo_in(counter_t)=proc_id
            t_halo_in(counter_t+1)=1
            t_halo_in(counter_t+2)=n_col
            counter_t=counter_t+3
            if (debug) write(0,*) me,' CDOVRBLD: Added into t_halo_in from recv',&
                 &proc_id,n_col,idx
          else if (desc_p%glob_to_loc(idx) < 0) Then
            if (debug) write(0,*) me,'Wrong input to cdovrbld??',&
                 &idx,desc_p%glob_to_loc(idx)
          End If
        End Do
        desc_p%matrix_data(psb_n_col_)=n_col

      end if

    end if
    t2 = mpi_wtime()
!!$    desc_p%matrix_data(psb_n_row_)=desc_p%matrix_data(psb_n_col_)
    !
    ! Ok, now we have a temporary halo with all the info for the
    ! next round. If we need to keep going, convert the halo format
    ! from temporary to final, so that we can work out the next iteration.
    ! This uses one of the convert_comm internals, i.e. we are doing
    ! the equivalent of a partial call to convert_comm
    !

    If (i_ovr < (n_ovr)) Then

      t_halo_in(counter_t)=-1

      if (debug) write(0,*) me,'Checktmp_o_i 1',tmp_ovr_idx(1:10)
      if (debug) write(0,*) me,'Calling Crea_Halo'

      call psi_crea_index(desc_p,t_halo_in,t_halo_out,.false.,nxch,nsnd,nrcv,info)

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
    t3 = mpi_wtime()
    tl = tl +(t2-t1)
    tch = tch +(t3-t2)

  End Do
  t1 = mpi_wtime()

  desc_p%matrix_data(psb_m_)=psb_cd_get_global_rows(desc_a)
  desc_p%matrix_data(psb_n_)=psb_cd_get_global_cols(desc_a)

  tmp_halo(counter_h:)=-1
  tmp_ovr_idx(counter_o:)=-1


  !
  ! At this point we have gathered all the indices in the halo at
  ! N levels of overlap. Just call cnv_dsc. This is
  ! the same routine as gets called inside CDASB.
  !

  if (debug) then
    write(0,*) 'psb_cdovrbld: converting indexes'
    call psb_barrier(ictxt)
  end if
  !.... convert comunication stuctures....
  ! Note that we have to keep local_rows until the very end, 
  ! because otherwise the halo build mechanism of cdasb
  ! will not work. 
  call psb_transfer(tmp_ovr_idx,desc_p%ovrlap_index,info)
  call psb_transfer(tmp_halo,desc_p%halo_index,info)
  call psb_cdasb(desc_p,info)
  desc_p%matrix_data(psb_n_row_)=desc_p%matrix_data(psb_n_col_)

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
  t2 = mpi_wtime()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == act_abort) then
    call psb_error(ictxt)
    return
  end if
  return

End Subroutine psb_dcdovrbld
