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
!    desc_p        - type(<psb_desc_type>).  The communication descriptor for the preconditioner.
!    desc_a        - type(<psb_desc_type>).  The communication descriptor.
!    a             - type(<psb_dspmat_type). The matrix upon which the preconditioner will be built.
!    l_tmp_halo    - integer.                Input estimate for allocation sizes.
!    l_tmp_ovr_idx - integer.                Input estimate for allocation sizes.
!    lworkr        - integer.                Input estimate for allocation sizes.
!    lworks        - integer.                Input estimate for allocation sizes.
!    info          - integer.                Eventually returns an error code
Subroutine psb_dcdovrbld(n_ovr,desc_p,desc_a,a,&
     &       l_tmp_halo,l_tmp_ovr_idx,lworks,lworkr,info)
  use psb_descriptor_type
  use psb_serial_mod
  Use psi_mod
  use psb_realloc_mod
  use psb_tools_mod, only : psb_cdprt
  use psb_error_mod
  use psb_const_mod
  use psb_penv_mod
  Implicit None

  include 'mpif.h'
  type(psb_dspmat_type),intent(in)  :: a
  type(psb_desc_type),intent(in)    :: desc_a
  type(psb_desc_type),intent(inout) :: desc_p
  integer,intent(in)                   :: n_ovr
  ! Input estimates for allocation sizes. Could we do without these two???
  Integer, Intent(in) :: l_tmp_halo,l_tmp_ovr_idx
  Integer, Intent(inout)  :: lworks, lworkr
  integer, intent(out)    :: info

  type(psb_dspmat_type)            :: blk


  Integer, Pointer :: tmp_halo(:),tmp_ovr_idx(:)

  Integer :: counter,counter_h, counter_o, counter_e,j,idx,gidx,proc,n_elem_recv,&
       & n_elem_send,tot_recv,tot_elem,n_col,m,ictxt,np,me,dl_lda,lwork,&
       & counter_t,n_elem,i_ovr,jj,i,proc_id,isz, mglob, glx,n_row, &
       & idxr, idxs, lx, iszr, err_act, icomm

  Integer,Pointer  :: halo(:),length_dl(:),works(:),workr(:),t_halo_in(:),&
       & t_halo_out(:),work(:),dep_list(:),temp(:)
  Integer,Pointer :: brvindx(:),rvsz(:), bsdindx(:),sdsz(:)
  integer :: pairtree(2)

  Logical,Parameter :: debug=.false.
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,t7, tl, tch
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  name='psb_cdovrbld'
  call psb_erractionsave(err_act)

  If(debug) Write(0,*)'cdovrbld begin'
  ictxt = desc_a%matrix_data(psb_ctxt_)

  Call psb_info(ictxt,me,np)

  call psb_nullify_sp(blk)

  Allocate(brvindx(np+1),rvsz(np),sdsz(np),bsdindx(np+1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  tl = 0.0
  tch = 0.0
  t4 = 0.0
  call psb_get_mpicomm(ictxt,icomm )

  mglob = desc_a%matrix_data(psb_m_)
  m     = desc_a%matrix_data(psb_n_row_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)
  if (debug) write(0,*) me,' On entry to CDOVRBLD n_col:',n_col

  dl_lda=np*5
  lwork=5*(5*np+2)*np+10
  Allocate(works(lworks),workr(lworkr),t_halo_in(3*Size(desc_p%halo_index)),&
       & t_halo_out(Size(desc_p%halo_index)), work(lwork),&
       & length_dl(np+1),dep_list(dl_lda*np),temp(lworkr),stat=info)
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
  halo => desc_a%halo_index

  Allocate(tmp_ovr_idx(l_tmp_ovr_idx),tmp_halo(l_tmp_halo),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  desc_p%ovrlap_elem(:) = -1
  tmp_ovr_idx(:)        = -1
  tmp_halo(:)           = -1
  counter_e             = 1
  tot_recv              = 0
  counter_h             = 1
  counter_o             = 1

  ! See comment in main loop below.
  call InitPairSearchTree(pairtree,info)
  if (info /= 0) then
    info=4010
    ch_err='InitPairSearhTree'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (debug) write(0,*) me,'Done InitPairSearchTree',info

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
!!$    t_halo_in(:) = -1

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
        info=-1
        call psb_errpush(info,name)
        goto 9999
      end If
      tot_recv=tot_recv+n_elem_recv
      if (debug) write(0,*) me,' CDOVRBLD tot_recv:',proc,n_elem_recv,tot_recv
      !
      ! While running through the column indices exchanged with other procs
      ! we have to keep track of which elements actually are overlap and halo 
      ! ones to record them in overlap_elem.  We do this by maintaining  
      ! an AVL balanced search tree: at each point counter_e is the next
      ! free index element. The search routine for gidx will return
      ! glx if gidx was already assigned a local index (glx<counter_e)
      ! but if gidx was a new index for this process, then it creates
      ! a new pair (gidx,counter_e), and glx==counter_e. In this case we
      ! need to record this for the overlap exchange. In the first iteration
      ! this gets filled with the first halo indices.
      !
      ! The format of the halo vector exists in two forms: 1. Temporary 
      ! 2. Assembled. In this loop we are using the (assembled) halo_in and 
      ! copying it into (temporary) tmp_halo; this is because tmp_halo will
      ! be enlarged with the new column indices received, and will reassemble
      ! everything for the next iteration.
      ! 

!!$      if (me==0) Write(0,*)'Loop ',size(halo), counter, n_elem_recv,n_elem_send 
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
        idx = halo(counter+psb_elem_recv_+j)
        If(idx > Size(desc_p%loc_to_glob)) then 
          info=-3
          call psb_errpush(info,name)
          goto 9999
        endif

        gidx = desc_p%loc_to_glob(idx)

        If((counter_o+2) > Size(tmp_ovr_idx)) Then
          isz = max((3*Size(tmp_ovr_idx))/2,(counter_o+3))
          if (debug) write(0,*) me,'Realloc tmp_ovr',isz
          call psb_realloc(isz,tmp_ovr_idx,info,pad=-1)
          if (info /= 0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        End If

        tmp_ovr_idx(counter_o)=proc
        tmp_ovr_idx(counter_o+1)=1
        tmp_ovr_idx(counter_o+2)=gidx
        tmp_ovr_idx(counter_o+3)=-1
        counter_o=counter_o+3

        If((counter_h+2) > Size(tmp_halo)) Then
          isz = max((3*Size(tmp_halo))/2,(counter_h+3))
          if (debug) write(0,*) me,'Realloc tmp_halo',isz
          call psb_realloc(isz,tmp_halo,info)
          if (info /= 0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        End If

        tmp_halo(counter_h)=proc
        tmp_halo(counter_h+1)=1
        tmp_halo(counter_h+2)=idx
        tmp_halo(counter_h+3)=-1

        counter_h=counter_h+3

        call SearchInsKeyVal(pairtree,gidx,counter_e,glx,info)
!!$        if (debug) write(0,*) 'From searchInsKey ',gidx,glx,counter_e,info
        if (info>=0) then
          If (glx < counter_e)  Then
            desc_p%ovrlap_elem(glx+psb_n_dom_ovr_)= &
                 &                 desc_p%ovrlap_elem(glx+psb_n_dom_ovr_)+1
          Else
            If((counter_e+2) > Size(desc_p%ovrlap_elem)) Then
              isz = max((3*Size(desc_p%ovrlap_elem))/2,(counter_e+3))
              if (debug) write(0,*) me,'Realloc ovr_El',isz
              call psb_realloc(isz,desc_p%ovrlap_elem,info)
              if (info /= 0) then
                info=4010
                ch_err='psrealloc'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            End If
!!$            if (debug) write(0,*) 'Adding into ovrlap ',gidx,glx,counter_e,info

            desc_p%ovrlap_elem(counter_e)=gidx
            desc_p%ovrlap_elem(counter_e+psb_n_dom_ovr_)=2
            desc_p%ovrlap_elem(counter_e+2)=-1
            counter_e  = counter_e + 2
          End If
        else
          write(0,*) me, 'Cdovrbld From SearchInsKeyVal: ',info
        endif
      Enddo
      if (debug) write(0,*) me,'Checktmp_o_i Loop Mid1',tmp_ovr_idx(1:10)
      counter   = counter+n_elem_recv

      !
      ! add send elements in halo_index into ovrlap_index
      !
      Do j=0,n_elem_send-1

        idx = halo(counter+psb_elem_send_+j)
        gidx = desc_p%loc_to_glob(idx)

        If((counter_o+2) > Size(tmp_ovr_idx)) Then
          isz = max((3*Size(tmp_ovr_idx))/2,(counter_o+3))
          if (debug) write(0,*) me,'Realloc tmp_ovr',isz
          call psb_realloc(isz,tmp_ovr_idx,info)
          if (info /= 0) then
            info=4010
            ch_err='psrealloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        End If

        tmp_ovr_idx(counter_o)=proc
        tmp_ovr_idx(counter_o+1)=1
        tmp_ovr_idx(counter_o+2)=gidx
        tmp_ovr_idx(counter_o+3)=-1
        counter_o=counter_o+3

        call SearchInsKeyVal(pairtree,gidx,counter_e,glx,info)
!!$          if (debug) write(0,*) 'From searchInsKey ',gidx,glx,counter_e,info
        if (info>=0) then
          If (glx < counter_e)  Then
            desc_p%ovrlap_elem(glx+psb_n_dom_ovr_)= &
                 &                 desc_p%ovrlap_elem(glx+psb_n_dom_ovr_)+1
          Else
            If((counter_e+2) > Size(desc_p%ovrlap_elem)) Then
              isz = max((3*Size(desc_p%ovrlap_elem))/2,(counter_e+3))
              if (debug) write(0,*) me,'Realloc ovr_el',isz
              call psb_realloc(isz,desc_p%ovrlap_elem,info)
              if (info /= 0) then
                info=4010
                ch_err='psb_realloc'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            End If
!!$            if (debug) write(0,*) 'Adding into ovrlap ',gidx,glx,counter_e,info
            desc_p%ovrlap_elem(counter_e)=gidx
            desc_p%ovrlap_elem(counter_e+psb_n_dom_ovr_)=2
            desc_p%ovrlap_elem(counter_e+2)=-1
            counter_e  = counter_e + 2
          End If
        else
          write(0,*) me,'Cdovrbld From SearchInsKeyVal: ',info
        endif

        !
        ! Prepare to exchange the halo rows with the other proc.
        !
        If (i_ovr < (n_ovr)) Then
          call psb_spinfo(psb_nzrowreq_,a,n_elem,info,iaux=idx)
          if (info /= 0) then
            info=4010
            ch_err='psb_spinfo'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if

          If((idxs+tot_elem+n_elem) > lworks) Then
            isz = max((3*lworks)/2,(idxs+tot_elem+n_elem))
            if (debug) write(0,*) me,'Realloc works',isz
            call psb_realloc(isz,works,info)
            if (info /= 0) then
              info=4010
              ch_err='psb_realloc'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            lworks = isz
          End If
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

          call psb_spgtblk(idx,a,blk,info)
          if (info /= 0) then
            info=4010
            ch_err='psb_spgtblk'
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
!!$          write(0,*) me,'Realloc temp',tot_elem+2
          if (tot_elem+2 > size(temp)) then 
!!$            write(0,*) me,'Realloc temp',tot_elem+2
            deallocate(temp)
            allocate(temp(tot_elem+2),stat=info)   
            if (info /= 0) then 
              call psb_errpush(4010,name,a_err='Allocate')
              goto 9999      
            end if
          endif
          Call mrgsrt(tot_elem,works(idxs+1),temp,info)
          If (info == 0) Call ireordv1(tot_elem,works(idxs+1),temp)
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

      Do i=1,iszr
        idx=workr(i)
        if (idx <1) then 
          write(0,*) me,'Error in CDOVRBLD ',idx,i,iszr
!!$          write(0,*) me, ' WORKR :',workr(1:iszr)
        else  If (desc_p%glob_to_loc(idx) < -np) Then
          !
          ! This is a new index. Assigning a local index as
          ! we receive the guarantees that all indices for HALO(I)
          ! will be less than those for HALO(J) whenever I<J
          !
          n_col=n_col+1
          proc_id=-desc_p%glob_to_loc(idx)-np-1
          If(n_col > Size(desc_p%loc_to_glob)) Then
            isz = 3*n_col/2
            if (debug) write(0,*) me,'Realloc loc_to_glob'
            call psb_realloc(isz,desc_p%loc_to_glob,info)
            if (info /= 0) then
              info=4010
              ch_err='psrealloc'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          End If
          desc_p%glob_to_loc(idx)=n_col
          desc_p%loc_to_glob(n_col)=idx
          If((counter_t+3) > Size(t_halo_in))Write(0,*)'bingo'
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
    end if
    t2 = mpi_wtime()
    n_row=m+tot_recv
    desc_p%matrix_data(psb_n_row_)=n_row
    desc_p%matrix_data(psb_n_col_)=n_col

    !
    ! Ok, now we have a temporary halo with all the info for the
    ! next round. If we need to keep going, convert the halo format
    ! from temporary to final, so that we can work out the next iteration.
    ! This uses one of the convert_comm internals, i.e. we are doing
    ! the equivalent of a partial call to convert_comm
    !

    If (i_ovr < (n_ovr)) Then

      If(lwork < (counter_t/3+np*3)) Then
        isz = max((3*lwork)/2,(counter_t/3+np*3))
        if (debug) write(0,*) me,'Realloc work',isz
        deallocate(work)
        allocate(work(isz),stat=info)
        if (info /= 0) then 
          call psb_errpush(4010,name,a_err='Allocate')
          goto 9999      
        end if

        lwork=size(work)
      Endif
      t_halo_in(counter_t)=-1

      if (debug) write(0,*) me,'Checktmp_o_i 1',tmp_ovr_idx(1:10)
      if (debug) write(0,*) me,'Calling Crea_Halo'

      call psi_crea_index(desc_p,t_halo_in,t_halo_out,.false.,info)

      if (debug) then 
        write(0,*) me,'Done Crea_Index'
        call psb_barrier(ictxt)
      end if
      if (debug) write(0,*) me,'Checktmp_o_i 2',tmp_ovr_idx(1:10)
      if (debug) write(0,*) me,'Done Crea_Halo'

      halo => t_halo_out
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
  call FreePairSearchTree(pairtree)

  desc_p%matrix_data(psb_m_)=desc_a%matrix_data(psb_m_)
  desc_p%matrix_data(psb_n_)=desc_a%matrix_data(psb_n_)

  tmp_halo(counter_h)=-1
  tmp_ovr_idx(counter_o)=-1


  !
  ! At this point we have gathered all the indices in the halo at
  ! N levels of overlap. Just call convert_comm. This is
  ! the same routine as gets called inside SPASB.
  !

  if (debug) then
    write(0,*) 'psb_cdasb: converting indexes'
    call psb_barrier(ictxt)
  end if
  !.... convert comunication stuctures....
  ! first the halo index
  call psi_crea_index(desc_p,tmp_halo,&
       & desc_p%halo_index,.false.,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='psi_crea_index')
    goto 9999
  end if

  ! then the overlap index
  call psi_crea_index(desc_p,tmp_ovr_idx,&
       & desc_p%ovrlap_index,.true.,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='psi_crea_index')
    goto 9999
  end if

  ! next is the ovrlap_elem index
  call psi_crea_ovr_elem(desc_p%ovrlap_index,desc_p%ovrlap_elem)

  ! finally bnd_elem
  call psi_crea_bnd_elem(desc_p,info)
  if(info /= 0) then
    call psb_errpush(4010,name,a_err='psi_crea_bnd_elem')
    goto 9999
  end if

  ! Ok, register into MATRIX_DATA &  free temporary work areas
  desc_p%matrix_data(psb_dec_type_) = psb_desc_asb_

  allocate(desc_p%lprm(1),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  desc_p%lprm(1) = 0

  if (debug) then 
    write(0,*) me,'Done Convert_comm'
    call psb_barrier(ictxt)
  end if

  if (.false.) then
    call psb_cdprt(70+me,desc_p,.false.)
  end if

  if (debug) write(0,*) me,'Done ConvertComm'

  Deallocate(works,workr,t_halo_in,t_halo_out,work,&
       & length_dl,dep_list,tmp_ovr_idx,tmp_halo,&
       & brvindx,rvsz,sdsz,bsdindx,temp,stat=info)
  if (info == 0) call psb_sp_free(blk,info)
  if (info /= 0) then
    info=4010
    ch_err='sp_free'
    call psb_errpush(info,name,a_err=ch_err)
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
