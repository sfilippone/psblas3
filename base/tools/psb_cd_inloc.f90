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
subroutine psb_cd_inloc(v, ictxt, desc, info)
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

  !locals
  Integer             :: counter,i,j,np,me,loc_row,err,&
       & loc_col,nprocs,n,itmpov, k,glx,&
       & l_ov_ix,l_ov_el,idx, flag_, err_act,m
  integer              :: int_err(5),exch(3)
  Integer, allocatable :: temp_ovrlap(:), tmpgidx(:,:)
  integer              :: debug_level, debug_unit
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
    call psb_max(ictxt,m)
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
  allocate(tmpgidx(m,2),stat=info) 
  if (info /=0) then 
    info=4000
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if



  tmpgidx=0
  flag_=1
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
  end do
  call psb_amx(ictxt,tmpgidx)
  if (info ==0) then 
    int_err(1) = count(tmpgidx(:,2) == 0)
    int_err(2) = m
    if (int_err(1)>0) then 
      info = 552
    end if
  end if

  if (info /= 0) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


  call psb_nullify_desc(desc)


  !count local rows number
  ! allocate work vector
  if (psb_cd_choose_large_state(ictxt,m)) then 
    allocate(desc%matrix_data(psb_mdata_size_),&
         &temp_ovrlap(m),stat=info)
    desc%matrix_data(psb_desc_size_) = psb_desc_large_
  else
    allocate(desc%glob_to_loc(m),desc%matrix_data(psb_mdata_size_),&
         &temp_ovrlap(m),stat=info)
    desc%matrix_data(psb_desc_size_) = psb_desc_normal_
  end if
  if (info /= 0) then     
    info=4025
    int_err(1)=2*m+psb_mdata_size_
    call psb_errpush(info,name,i_err=int_err,a_err='integer')
    goto 9999
  endif

  desc%matrix_data(psb_m_)        = m
  desc%matrix_data(psb_n_)        = n
  ! This has to be set BEFORE any call to SET_BLD
  desc%matrix_data(psb_ctxt_)     = ictxt
  call psb_get_mpicomm(ictxt,desc%matrix_data(psb_mpi_c_))


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  starting main loop' ,info
  counter = 0
  itmpov  = 0
  temp_ovrlap(:) = -1

  !
  ! We have to decide whether we have a "large" index space.
  !
  if (psb_cd_choose_large_state(ictxt,m)) then 
    !
    ! Yes, we do have a large index space. Therefore we are 
    ! keeping on the local process a map of only the global 
    ! indices ending up here; this map is stored in an AVL
    ! tree during the build stage, so as to guarantee log-time
    ! serch and insertion of new items. At assembly time it 
    ! is transferred to a series of ordered linear lists, 
    ! hashed by the low order bits of the entries.
    !
    do i=1,m

      if (((tmpgidx(i,1)-flag_) > np-1).or.((tmpgidx(i,1)-flag_) < 0)) then
        info=580
        int_err(1)=3
        int_err(2)=tmpgidx(i,1) - flag_
        int_err(3)=i
        exit
      end if

      if ((tmpgidx(i,1)-flag_) == me) then
        ! this point belongs to me
        counter=counter+1
      end if
    enddo


    loc_row=counter
    ! check on parts function
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),':  End main loop:' ,loc_row,itmpov,info

    if (info /= 0) then 
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),':  error check:' ,err

    ! estimate local cols number 
    loc_col = min(2*loc_row,m)

    allocate(desc%loc_to_glob(loc_col), desc%lprm(1),&
         & desc%ptree(2),stat=info)  
    if (info == 0) call InitPairSearchTree(desc%ptree,info)
    if (info /= 0) then
      info=4025
      int_err(1)=loc_col
      call psb_errpush(info,name,i_err=int_err,a_err='integer')
      goto 9999
    end if

    ! set LOC_TO_GLOB array to all "-1" values
    desc%lprm(1) = 0
    desc%loc_to_glob(:) = -1
    k = 0
    do i=1,m
      if ((tmpgidx(i,1)-flag_) == me) then
        k = k + 1 
        desc%loc_to_glob(k) = i
        call SearchInsKeyVal(desc%ptree,i,k,glx,info)
      endif
    enddo
    if (k /= loc_row) then 
      write(0,*) 'Large cd init: ',k,loc_row
    end if

    if (info /= 0) then 
      info=4000
      call psb_errpush(info,name)
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

      if ((tmpgidx(i,1)-flag_) == me) then
        ! this point belongs to me
        counter=counter+1
        desc%glob_to_loc(i) = counter
      else
        desc%glob_to_loc(i) = -(np+(tmpgidx(i,1)-flag_)+1)
      end if
    enddo


    loc_row=counter
    ! check on parts function
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),':  End main loop:' ,loc_row,itmpov,info

    if (info /= 0) then 
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),':  error check:' ,err

    ! estimate local cols number 
    loc_col = min(2*loc_row,m)

    allocate(desc%loc_to_glob(loc_col),&
         &desc%lprm(1),stat=info)  
    if (info /= 0) then
      info=4025
      int_err(1)=loc_col
      call psb_errpush(info,name,i_err=int_err,a_err='integer')
      goto 9999
    end if

    ! set LOC_TO_GLOB array to all "-1" values
    desc%lprm(1) = 0
    desc%loc_to_glob(:) = -1
    do i=1,m
      k = desc%glob_to_loc(i) 
      if (k > 0) then 
        desc%loc_to_glob(k) = i
      endif
    enddo

  end if


  call psi_bld_tmpovrl(temp_ovrlap,desc,info)

  deallocate(temp_ovrlap,stat=info)
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
