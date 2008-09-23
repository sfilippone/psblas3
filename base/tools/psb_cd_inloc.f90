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
! File: psb_cd_inloc.f90
!
! Subroutine: psb_cd_inloc
!    Allocate descriptor with a local vector V containing the list 
!    of indices that are assigned to the current process. The global size 
!    is equal to the largest index found on any process. 
! 
! Arguments: 
!    v       - integer, dimension(:).         The array containg the partitioning scheme.
!    ictxt - integer.                         The communication context.
!    desc  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
subroutine psb_cd_inloc(v, ictxt, desc, info, globalcheck)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psi_mod
  use psb_penv_mod
  implicit None
  !....Parameters...
  Integer, intent(in)               :: ictxt, v(:)
  integer, intent(out)              :: info
  type(psb_desc_type), intent(out)  :: desc
  logical, intent(in), optional     :: globalcheck

  !locals
  Integer             :: i,j,np,me,loc_row,err,&
       & loc_col,nprocs,n, k,glx,nlu,&
       & idx, flag_, err_act,m, novrl, norphan,&
       & npr_ov, itmpov, i_pnt, nrt
  integer              :: int_err(5),exch(3)
  Integer, allocatable :: temp_ovrlap(:), tmpgidx(:,:), vl(:),&
       & nov(:), ov_idx(:,:)
  integer              :: debug_level, debug_unit
  logical              :: check_, islarge
  character(len=20)    :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  err=0
  name = 'psb_cd_inloc'
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': start',np


  loc_row = size(v)
  if (.false.) then 
    m = loc_row
    call psb_sum(ictxt,m)
  else
    m = maxval(v)
    nrt = loc_row
    call psb_sum(ictxt,nrt)
    call psb_max(ictxt,m)
  end if
  if (present(globalcheck)) then 
    check_ = globalcheck
  else
    check_ = .true.
  end if

  n = m

  !... check m and n parameters....
  if (m < 1) then
    info = 10
    int_err(1) = 1
    int_err(2) = m
  else if (n < 1) then
    info = 10
    int_err(1) = 2
    int_err(2) = n
  endif

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (me == psb_root_) then
    exch(1)=m
    exch(2)=n
    exch(3)=psb_cd_get_large_threshold()
    call psb_bcast(ictxt,exch(1:3),root=psb_root_)
  else
    call psb_bcast(ictxt,exch(1:3),root=psb_root_)
    if (exch(1) /= m) then
      err=550
      int_err(1)=1
      call psb_errpush(err,name,int_err)
      goto 9999
    else if (exch(2) /= n) then
      err=550
      int_err(1)=2
      call psb_errpush(err,name,int_err)
      goto 9999
    endif
    call psb_cd_set_large_threshold(exch(3))
  endif

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  doing global checks'  

  islarge = psb_cd_choose_large_state(ictxt,m)

  allocate(vl(loc_row),stat=info) 
  if (info /=0) then 
    info=4000
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  
  !
  ! Checks for valid input:
  !  1. legal range
  !  2. no orphans
  !  3. any overlap?
  ! Checks 2 and 3 are controlled by globalcheck
  !  

  if (check_.or.(.not.islarge)) then 
    allocate(tmpgidx(m,2),stat=info) 
    if (info /=0) then 
      info=4000
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    tmpgidx = 0
    flag_   = 1
    do i=1,loc_row
      if ((v(i)<1).or.(v(i)>m)) then 
        info = 551
        int_err(1) = i
        int_err(2) = v(i)
        int_err(3) = loc_row
        int_err(4) = m
      else
        tmpgidx(v(i),1) = me+flag_
        tmpgidx(v(i),2) = 1
      endif
      vl(i) = v(i) 
    end do

    if (info ==0) then 
      call psb_amx(ictxt,tmpgidx(:,1))
      call psb_sum(ictxt,tmpgidx(:,2))
      novrl   = 0
      npr_ov  = 0
      norphan = 0
      do i=1, m
        if (tmpgidx(i,2) < 1) then 
          norphan = norphan + 1 
        else if (tmpgidx(i,2) > 1) then 
          novrl  = novrl + 1 
          npr_ov = npr_ov + tmpgidx(i,2)
        end if
      end do
      if (norphan > 0) then 
        int_err(1) = norphan
        int_err(2) = m
        info = 552
      end if
    end if
  else
    novrl   = 0
    norphan = 0
    npr_ov  = 0
    do i=1,loc_row
      if ((v(i)<1).or.(v(i)>m)) then 
        info = 551
        int_err(1) = i
        int_err(2) = v(i)
        int_err(3) = loc_row
        int_err(4) = m
        exit
      endif
      vl(i) = v(i) 
    end do
    
    if ((m /= nrt).and.(me==psb_root_))  then 
      write(0,*) trim(name),' Warning: globalcheck=.false., but there is a mismatch'
      write(0,*) trim(name),'        : in the global sizes!',m,nrt
    end if
  end if

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  !
  ! Now sort the input items, and eliminate duplicates
  ! (unlikely, but possible)
  !
  call psb_msort_unique(vl,nlu)
  loc_row = nlu

  call psb_nullify_desc(desc)

  !
  ! Figure out overlap in the input
  ! 
  if (novrl > 0) then 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': code for NOVRL>0',novrl,npr_ov

    allocate(nov(0:np),ov_idx(npr_ov,2),stat=info)
    if (info /= 0) then 
      info=4025
      int_err(1)=np + 2*npr_ov
      call psb_errpush(info,name,i_err=int_err,a_err='integer')
      goto 9999
    endif
    nov=0
    do i=1, nlu
      k = vl(i) 
      if (tmpgidx(k,2) > 1) then 
        nov(me) = nov(me) + 1 
      end if
    end do
    call psb_sum(ictxt,nov)
    nov(1:np) = nov(0:np-1)
    nov(0) = 1
    do i=1, np
      nov(i) = nov(i) + nov(i-1)
    end do
    ov_idx = 0
    j = nov(me)
    do i=1, nlu
      k = vl(i) 
      if (tmpgidx(k,2) > 1) then 
        ov_idx(j,1) = k
        ov_idx(j,2) = me
        j = j + 1 
      end if
    end do

    if (j /= nov(me+1)) then 
      info=4001
      call psb_errpush(info,name,a_err='overlap count')
      goto 9999
    end if
    call psb_max(ictxt,ov_idx)
    call psb_msort(ov_idx(:,1),ix=ov_idx(:,2),flag=psb_sort_keep_idx_)

  end if

  ! allocate work vector
  if (islarge) then 
    allocate(desc%matrix_data(psb_mdata_size_),&
         &temp_ovrlap(max(1,2*loc_row)),desc%lprm(1),&
         & stat=info)
    if (info == 0) then 
      desc%lprm(1)        = 0   
      desc%matrix_data(:) = 0
      desc%idxmap%state   = psb_desc_large_
    end if
  else
    allocate(desc%idxmap%glob_to_loc(m),desc%matrix_data(psb_mdata_size_),&
         &temp_ovrlap(max(1,2*loc_row)),desc%lprm(1),&
         & stat=info)
    if (info == 0) then 
      desc%lprm(1)        = 0   
      desc%matrix_data(:) = 0
      desc%idxmap%state   = psb_desc_normal_
    end if
  end if
  if (info /= 0) then     
    info=4025
    int_err(1)=2*m+psb_mdata_size_
    call psb_errpush(info,name,i_err=int_err,a_err='integer')
    goto 9999
  endif

  ! estimate local cols number 
  loc_col = min(2*loc_row,m)

  allocate(desc%idxmap%loc_to_glob(loc_col),stat=info)  
  if (info /= 0) then
    info=4025
    int_err(1)=loc_col
    call psb_errpush(info,name,i_err=int_err,a_err='integer')
    goto 9999
  end if
  desc%idxmap%loc_to_glob(:) = -1
  temp_ovrlap(:) = -1
  desc%matrix_data(psb_m_)        = m
  desc%matrix_data(psb_n_)        = n
  ! This has to be set BEFORE any call to SET_BLD
  desc%matrix_data(psb_ctxt_)     = ictxt
  call psb_get_mpicomm(ictxt,desc%matrix_data(psb_mpi_c_))


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  starting main loop' ,info

  !
  ! We have to decide whether we have a "large" index space.
  !
  if (islarge) then 
    !
    ! Yes, we do have a large index space. Therefore we are 
    ! keeping on the local process a map of only the global 
    ! indices ending up here; this map is stored in an AVL
    ! tree during the build stage, so as to guarantee log-time
    ! serch and insertion of new items. At assembly time it 
    ! is transferred to a series of ordered linear lists, 
    ! hashed by the low order bits of the entries.
    !

    if (info /= 0) then
      info=4025
      int_err(1)=loc_col
      call psb_errpush(info,name,i_err=int_err,a_err='integer')
      goto 9999
    end if

    ! Use sorted indices to fill in loc_to_glob
    j      = 1
    itmpov = 0
    do k=1, loc_row
      i = vl(k)
      desc%idxmap%loc_to_glob(k) = i

      if (check_) then 
        nprocs = tmpgidx(i,2) 
        if (nprocs > 1) then 
          do 
            if (j > size(ov_idx,dim=1)) then 
              info=4001
              call psb_errpush(info,name,a_err='search ov_idx')
              goto 9999
            end if
            if (ov_idx(j,1) == i) exit
            j = j + 1 
          end do
          call psb_ensure_size((itmpov+3+nprocs),temp_ovrlap,info,pad=-1)
          if (info /= 0) then
            info=4010
            call psb_errpush(info,name,a_err='psb_ensure_size')
            goto 9999
          end if
          itmpov = itmpov + 1
          temp_ovrlap(itmpov) = i
          itmpov = itmpov + 1
          temp_ovrlap(itmpov) = nprocs
          temp_ovrlap(itmpov+1:itmpov+nprocs) = ov_idx(j:j+nprocs-1,2)
          itmpov = itmpov + nprocs
        end if
      end if

    enddo

    if (info /= 0) then 
      info=4001
      call psb_errpush(info,name,a_err='insert loop')
      goto 9999
    endif

  else 

    !
    ! No, we don't have a large index space. Therefore we can  
    ! afford to keep on the local process a map of all global 
    ! indices; for those we know are here we immediately store 
    ! the corresponding local index, for the others we encode 
    ! the index of the process owning them, so that during the 
    ! insertion phase we can use the information to build the 
    ! data exchange lists "on-the-fly".
    !

    do i=1,m

      if (((tmpgidx(i,1)-flag_) > np-1).or.((tmpgidx(i,1)-flag_) < 0)) then
        info=580
        int_err(1)=3
        int_err(2)=tmpgidx(i,1) - flag_
        int_err(3)=i
        exit
      end if

      desc%idxmap%glob_to_loc(i) = -(np+(tmpgidx(i,1)-flag_)+1)
    enddo

    if (info /= 0) then 
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! Use sorted indices to fill in loc_to_glob
    j      = 1
    itmpov = 0
    do k=1, loc_row
      i = vl(k)
      desc%idxmap%loc_to_glob(k) = i
      desc%idxmap%glob_to_loc(i) = k

      nprocs = tmpgidx(i,2) 
      if (nprocs > 1) then 
        do 
          if (j > size(ov_idx,dim=1)) then 
            info=4001
            call psb_errpush(info,name,a_err='search ov_idx')
            goto 9999
          end if
          if (ov_idx(j,1) == i) exit
          j = j + 1 
        end do
        call psb_ensure_size((itmpov+3+nprocs),temp_ovrlap,info,pad=-1)
        if (info /= 0) then
          info=4010
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if
        itmpov = itmpov + 1
        temp_ovrlap(itmpov) = i
        itmpov = itmpov + 1
        temp_ovrlap(itmpov) = nprocs
        temp_ovrlap(itmpov+1:itmpov+nprocs) = ov_idx(j:j+nprocs-1,2)
        itmpov = itmpov + nprocs
      end if
    enddo

  end if

  call psi_bld_tmpovrl(temp_ovrlap,desc,info)

  if (info == 0) deallocate(temp_ovrlap,vl,stat=info)
  if ((info == 0).and.(allocated(tmpgidx)))&
       &  deallocate(tmpgidx,stat=info)
  if ((info == 0) .and.(allocated(ov_idx)))  &
       & deallocate(ov_idx,nov,stat=info)
  if (info /= 0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  endif

  ! set fields in desc%MATRIX_DATA....
  desc%matrix_data(psb_n_row_)  = loc_row
  desc%matrix_data(psb_n_col_)  = loc_row

  call psb_realloc(max(1,loc_row/2),desc%halo_index, info)
  if (info == 0) call psb_realloc(1,desc%ext_index, info)
  if (info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='psb_realloc')
    Goto 9999
  end if
  desc%matrix_data(psb_pnt_h_) = 1
  desc%halo_index(:)           = -1
  desc%ext_index(:)            = -1

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_cd_set_bld(desc,info)
  if (info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='psb_cd_set_bld')
    Goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cd_inloc
