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
! File: psb_cd_inloc.f90
!
! Subroutine: psb_cd_inloc
!    Allocate descriptor with a local vector V containing the list 
!    of indices that are assigned to the current process. The global size 
!    is equal to the largest index found on any process. 
! 
! Arguments: 
!    v       - integer(psb_ipk_), dimension(:).         The array containg the partitioning scheme.
!    ictxt - integer.                         The communication context.
!    desc  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
subroutine psb_cd_inloc(v, ictxt, desc, info, globalcheck,idx,usehash)
  use psb_base_mod
  use psi_mod
  use psb_repl_map_mod
  use psb_list_map_mod
  use psb_hash_map_mod
  implicit None
  !....Parameters...
  integer(psb_ipk_), intent(in)               :: ictxt
  integer(psb_lpk_), intent(in)               :: v(:)
  integer(psb_ipk_), intent(out)              :: info
  type(psb_desc_type), intent(out)  :: desc
  logical, intent(in), optional     :: globalcheck,usehash
  integer(psb_ipk_), intent(in), optional     :: idx(:)

  !locals
  integer(psb_ipk_) :: i,j,np,me,loc_row,err,&
       & loc_col,nprocs,k,glx,nlu,&
       & flag_, err_act, novrl, norphan,&
       & npr_ov, itmpov, i_pnt
  integer(psb_lpk_) :: m, n, nrt, il
  integer(psb_lpk_) :: l_err(5),exch(3)
  integer(psb_ipk_), allocatable :: tmpgidx(:,:), &
       & nov(:), ov_idx(:,:), temp_ovrlap(:)
  integer(psb_lpk_), allocatable :: vl(:), ix(:), l_temp_ovrlap(:)
  integer(psb_ipk_)  :: debug_level, debug_unit
  real(psb_dpk_)     :: t0, t1, t2, t3, t4, t5
  logical, parameter :: debug_size=.false.
  logical            :: do_timings=.false.
  logical            :: check_, islarge, usehash_
  character(len=20)  :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  err=0
  name = 'psb_cd_inloc'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': start',np
  if (do_timings) then 
    call psb_barrier(ictxt)
    t0 = psb_wtime()
  end if
  loc_row = size(v)
  m       = maxval(v)
  nrt     = loc_row
  call psb_sum(ictxt,nrt)
  call psb_max(ictxt,m)
  
  if (present(globalcheck)) then 
    check_ = globalcheck
  else
    check_ = .false.
  end if
  if (present(usehash)) then 
    usehash_ = usehash
  else
    usehash_ = .false.
  end if

  n = m

  !... check m and n parameters....
  if (m < 1) then
    info = psb_err_iarg_neg_
    l_err(1) = 1
    l_err(2) = m
  else if (n < 1) then
    info = psb_err_iarg_neg_
    l_err(1) = 2
    l_err(2) = n
  endif

  if (info /= psb_success_) then 
    call psb_errpush(info,name,l_err=l_err)
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
      l_err(1)=1
      call psb_errpush(err,name,l_err=l_err)
      goto 9999
    else if (exch(2) /= n) then
      err=550
      l_err(1)=2
      call psb_errpush(err,name,l_err=l_err)
      goto 9999
    endif
    call psb_cd_set_large_threshold(exch(3))
  endif
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  doing global checks'  

  islarge = psb_cd_is_large_size(m)

  allocate(vl(max(loc_row,ione)),ix(max(loc_row,ione)),stat=info) 
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name,l_err=l_err)
    goto 9999
  end if
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': sizes',loc_row,m,nrt,check_

  !
  ! Checks for valid input:
  !  1. legal range
  !  2. no orphans
  !  3. any overlap?
  ! Checks 2 and 3 are controlled by globalcheck
  !  

  if (check_.or.(.not.islarge)) then
    if (debug_size) &
         & write(debug_unit,*) me,' ',trim(name),': Going for global checks'
  
    allocate(tmpgidx(m,2),stat=info) 
    if (info /= psb_success_) then 
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name,l_err=l_err)
      goto 9999
    end if
    tmpgidx = 0
    flag_   = 1
    do i=1,loc_row
      if ((v(i)<1).or.(v(i)>m)) then 
        info = psb_err_entry_out_of_bounds_
        l_err(1) = i
        l_err(2) = v(i)
        l_err(3) = loc_row
        l_err(4) = m
      else
        tmpgidx(v(i),1) = me+flag_
        tmpgidx(v(i),2) = 1
      endif
      vl(i) = v(i) 
    end do

    if (info == psb_success_) then 
      call psb_amx(ictxt,tmpgidx(:,1))
      call psb_sum(ictxt,tmpgidx(:,2))
      novrl   = 0
      npr_ov  = 0
      norphan = 0
      do il=1, m
        if (tmpgidx(il,2) < 1) then 
          norphan = norphan + 1 
        else if (tmpgidx(il,2) > 1) then 
          novrl  = novrl + 1 
          npr_ov = npr_ov + tmpgidx(il,2)
        end if
      end do
      if (norphan > 0) then 
        l_err(1) = norphan
        l_err(2) = m
        info = psb_err_inconsistent_index_lists_
      end if
    end if
    
  else
    novrl   = 0
    norphan = 0
    npr_ov  = 0
    do i=1,loc_row
      if ((v(i)<1).or.(v(i)>m)) then 
        info = psb_err_entry_out_of_bounds_
        l_err(1) = i
        l_err(2) = v(i)
        l_err(3) = loc_row
        l_err(4) = m
        exit
      endif
      vl(i) = v(i) 
    end do

    if ((m /= nrt).and.(me == psb_root_))  then 
      write(psb_err_unit,*) trim(name),' Warning: globalcheck=.false., but there is a mismatch'
      write(psb_err_unit,*) trim(name),'        : in the global sizes!',m,nrt
    end if
  end if
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': After global checks '

  if (do_timings) then 
    call psb_barrier(ictxt)
    t1 = psb_wtime()
  end if

  if (info /= psb_success_) then 
    call psb_errpush(info,name,l_err=l_err)
    goto 9999
  end if

  ! Sort, eliminate duplicates, then
  ! scramble back into original position.
  ix(1) = -1 
  if (present(idx)) then 
    if (size(idx) >= loc_row) then 
      do i=1, loc_row
        ix(i) = idx(i) 
      end do
    end if
  end if
  if (ix(1) == -1) then 
    do i=1, loc_row
      ix(i) = i 
    end do
  end if
  call psb_msort(vl,ix,flag=psb_sort_keep_idx_)
  nlu = min(1,loc_row)
  do i=2,loc_row
    if (vl(i) /= vl(nlu)) then
      nlu = nlu + 1 
      vl(nlu) = vl(i)
      ix(nlu) = ix(i)
    end if
  end do
  call psb_msort(ix(1:nlu),vl(1:nlu),flag=psb_sort_keep_idx_)

  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': After sort ',nlu

  call psb_nullify_desc(desc)
  if (do_timings) then 
    call psb_barrier(ictxt)
    t2 = psb_wtime()
  end if

  !
  ! Figure out overlap in the input.
  ! Note: the code above guarantees that if mpgidx was not allocated,
  !       then novrl = 0, hence all accesses to tmpgidx
  !       are safe. 
  ! 
  if (novrl > 0) then
    if (debug_size) &
         & write(debug_unit,*) me,' ',trim(name),': Check overlap '
    
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': code for NOVRL>0',novrl,npr_ov

    allocate(nov(0:np),ov_idx(npr_ov,2),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      l_err(1)=np + 2*npr_ov
      call psb_errpush(info,name,l_err=l_err,a_err='integer')
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
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err='overlap count')
      goto 9999
    end if
    call psb_max(ictxt,ov_idx)
    call psb_msort(ov_idx(:,1),ix=ov_idx(:,2),flag=psb_sort_keep_idx_)

  end if
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': Done overlap '
    

  ! allocate work vector
  allocate(l_temp_ovrlap(max(1,2*loc_row)),desc%lprm(1),&
       & stat=info)
  if (info == psb_success_) then 
    desc%lprm(1)        = 0   
  end if
  if (info /= psb_success_) then     
    info=psb_err_alloc_request_
    l_err(1)=2*m+psb_mdata_size_
    call psb_errpush(info,name,l_err=l_err,a_err='integer')
    goto 9999
  endif

  l_temp_ovrlap(:) = -1

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  starting main loop' ,info

  j      = 1
  itmpov = 0
  if (check_) then 
    do k=1, loc_row
      il = v(k)
      nprocs = tmpgidx(il,2) 
      if (nprocs > 1) then 
        do 
          if (j > size(ov_idx,dim=1)) then 
            info=psb_err_internal_error_
            call psb_errpush(info,name,a_err='search ov_idx')
            goto 9999
          end if
          if (ov_idx(j,1) == i) exit
          j = j + 1 
        end do
        call psb_ensure_size((itmpov+3+nprocs),l_temp_ovrlap,info,pad=-1_psb_lpk_)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_ensure_size')
          goto 9999
        end if
        itmpov = itmpov + 1
        l_temp_ovrlap(itmpov) = il
        itmpov = itmpov + 1
        l_temp_ovrlap(itmpov) = nprocs
        l_temp_ovrlap(itmpov+1:itmpov+nprocs) = ov_idx(j:j+nprocs-1,2)
        itmpov = itmpov + nprocs
      end if

    end do
  end if
  if (do_timings) then 
    call psb_barrier(ictxt)
    t3 = psb_wtime()
  end if
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': Allocate indxmap'
    
  if (np == 1) then 
    allocate(psb_repl_map :: desc%indxmap, stat=info)
  else
    if (islarge.or.usehash_) then 
      allocate(psb_hash_map :: desc%indxmap, stat=info)
    else
      allocate(psb_list_map :: desc%indxmap, stat=info)
    end if
  end if

  select type(aa => desc%indxmap) 
  type is (psb_repl_map) 
    call aa%repl_map_init(ictxt,m,info)
  class default 
    call aa%init(ictxt,vl(1:nlu),info)
  end select

  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': Done init indxmap'

  
  if (do_timings) then 
    call psb_barrier(ictxt)
    t4 = psb_wtime()
  end if

  !
  ! Now that we have initialized indxmap we can convert the
  ! indices to local numbering.
  !
  block
    integer(psb_ipk_) :: i,nprocs
    allocate(temp_ovrlap(size(l_temp_ovrlap)),stat=info)
    if (info == psb_success_) then
      temp_ovrlap = -1
      i = 1
      do while (l_temp_ovrlap(i) /= -1) 
        call desc%indxmap%g2l(l_temp_ovrlap(i),temp_ovrlap(i),info)
        i              = i + 1
        temp_ovrlap(i) = l_temp_ovrlap(i)
        nprocs         = temp_ovrlap(i)
        temp_ovrlap(i+1:i+nprocs) = l_temp_ovrlap(i+1:i+nprocs)
        i       = i + 1
        i       = i + nprocs     
      enddo
    end if
  end block
  if (info == psb_success_) call psi_bld_tmpovrl(temp_ovrlap,desc,info)
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': Done bld_tmpovrl'

  if (info == psb_success_) deallocate(temp_ovrlap,vl,ix,stat=info)
  if ((info == psb_success_).and.(allocated(tmpgidx)))&
       &  deallocate(tmpgidx,stat=info)
  if ((info == psb_success_) .and.(allocated(ov_idx)))  &
       & deallocate(ov_idx,nov,stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (do_timings) then 
    call psb_barrier(ictxt)
    t5 = psb_wtime()

    t5 = t5 - t4
    t4 = t4 - t3
    t3 = t3 - t2
    t2 = t2 - t1
    t1 = t1 - t0
    call psb_amx(ictxt,t1)
    call psb_amx(ictxt,t2)
    call psb_amx(ictxt,t3)
    call psb_amx(ictxt,t4)
    call psb_amx(ictxt,t5)
    if (me==0) then
      write(0,*) 'CD_INLOC Timings: '
      write(0,*) '    Phase 1     : ', t1
      write(0,*) '    Phase 2     : ', t2
      write(0,*) '    Phase 3     : ', t3
      write(0,*) '    Phase 4     : ', t4
      write(0,*) '    Phase 5     : ', t5
    end if
  end if
  if (debug_size) &
       & write(debug_unit,*) me,' ',trim(name),': Done cd_inloc'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_cd_inloc
