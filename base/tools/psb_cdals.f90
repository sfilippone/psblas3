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
! File: psb_cdals.f90
!
! Subroutine: psb_cdals
!    Allocate descriptor
!    and checks correctness of PARTS subroutine
! 
! Arguments: 
!    m       - integer.                       The number of rows.
!    n       - integer.                       The number of columns.
!    parts   - external subroutine.           The routine that contains the 
!                                                 partitioning scheme.
!    ictxt - integer.                         The communication context.
!    desc  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       Error code (if any).
subroutine psb_cdals(m, n, parts, ictxt, desc, info)
  use psb_error_mod
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_penv_mod
  implicit None
  include 'parts.fh'
  !....Parameters...
  Integer, intent(in)                 :: M,N,ictxt
  Type(psb_desc_type), intent(out)    :: desc
  integer, intent(out)                :: info

  !locals
  Integer              :: counter,i,j,np,me,loc_row,err,loc_col,nprocs,&
       & l_ov_ix,l_ov_el,idx, err_act, itmpov, k, glx 
  integer              :: int_err(5),exch(3)
  integer, allocatable  :: prc_v(:), temp_ovrlap(:)
  integer              :: debug_level, debug_unit
  character(len=20)    :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  err=0
  name = 'psb_cdall'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  

  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': ',np
  !     ....verify blacs grid correctness..

  !... check m and n parameters....
  if (m < 1) then
    info = 10
    err=info
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(err,name,int_err)
    goto 9999
  else if (n < 1) then
    info = 10
    err=info
    int_err(1) = 2
    int_err(2) = n
    call psb_errpush(err,name,int_err)
    goto 9999
  endif


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  doing global checks'  
  !global check on m and n parameters
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

  call psb_nullify_desc(desc)

  !count local rows number
  ! allocate work vector
  if (psb_cd_choose_large_state(ictxt,m)) then 
    allocate(desc%matrix_data(psb_mdata_size_),&
         & temp_ovrlap(m),prc_v(np),stat=info)
    desc%matrix_data(psb_desc_size_) = psb_desc_large_
  else
    allocate(desc%glob_to_loc(m),desc%matrix_data(psb_mdata_size_),&
         & temp_ovrlap(m),prc_v(np),stat=info)
    desc%matrix_data(psb_desc_size_) = psb_desc_normal_
  end if
  if (info /= 0) then     
    info=4025
    err=info
    int_err(1)=2*m+psb_mdata_size_+np
    call psb_errpush(err,name,int_err,a_err='integer')
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
    loc_col = (m+np-1)/np
    loc_col = min(2*loc_col,m)
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
      if (info == 0) then
        call parts(i,m,np,prc_v,nprocs)
        if (nprocs > np) then
          info=570
          int_err(1)=3
          int_err(2)=np
          int_err(3)=nprocs
          int_err(4)=i
          err=info
          call psb_errpush(err,name,int_err)
          goto 9999
        else if (nprocs <= 0) then
          info=575
          int_err(1)=3
          int_err(2)=nprocs
          int_err(3)=i
          err=info
          call psb_errpush(err,name,int_err)
          goto 9999
        else
          do j=1,nprocs
            if ((prc_v(j) > np-1).or.(prc_v(j) < 0)) then
              info=580
              int_err(1)=3
              int_err(2)=prc_v(j)
              int_err(3)=i
              err=info
              call psb_errpush(err,name,int_err)
              goto 9999
            end if
          end do
        endif
        j=1
        do 
          if (j > nprocs) exit
          if (prc_v(j) == me) exit
          j=j+1
        enddo
        
        if (j <= nprocs) then 
          if (prc_v(j) == me) then
            ! this point belongs to me
            k = k + 1 
            call psb_ensure_size((k+1),desc%loc_to_glob,info,pad=-1)
            if (info /= 0) then
              info=4010
              call psb_errpush(info,name,a_err='psb_ensure_size')
              goto 9999
            end if
            desc%loc_to_glob(k) = i
            call SearchInsKeyVal(desc%ptree,i,k,glx,info)
            if (nprocs > 1)  then
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
              temp_ovrlap(itmpov+1:itmpov+nprocs) = prc_v(1:nprocs)
              itmpov = itmpov + nprocs
            endif
          end if
        end if        
      end if
    enddo
    if (info /= 0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    endif
    loc_row = k 

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
      if (info == 0) then
        call parts(i,m,np,prc_v,nprocs)
        if (nprocs > np) then
          info=570
          int_err(1)=3
          int_err(2)=np
          int_err(3)=nprocs
          int_err(4)=i
          err=info
          call psb_errpush(err,name,int_err)
          goto 9999
        else if (nprocs <= 0) then
          info=575
          int_err(1)=3
          int_err(2)=nprocs
          int_err(3)=i
          err=info
          call psb_errpush(err,name,int_err)
          goto 9999
        else
          do j=1,nprocs
            if ((prc_v(j) > np-1).or.(prc_v(j) < 0)) then
              info=580
              int_err(1)=3
              int_err(2)=prc_v(j)
              int_err(3)=i
              err=info
              call psb_errpush(err,name,int_err)
              goto 9999
            end if
          end do
        endif
        desc%glob_to_loc(i) = -(np+prc_v(1)+1)
        j=1
        do 
          if (j > nprocs) exit
          if (prc_v(j) == me) exit
          j=j+1
        enddo
        if (j <= nprocs) then 
          if (prc_v(j) == me) then
            ! this point belongs to me
            counter=counter+1
            desc%glob_to_loc(i) = counter
            if (nprocs > 1)  then
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
              temp_ovrlap(itmpov+1:itmpov+nprocs) = prc_v(1:nprocs)
              itmpov = itmpov + nprocs
            endif
          end if
        end if
      endif
    enddo
    ! estimate local cols number 
    loc_row=counter
    loc_col=min(2*loc_row,m)

    allocate(desc%loc_to_glob(loc_col),&
         &desc%lprm(1),stat=info)  
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
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

  ! check on parts function
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  End main loop:' ,loc_row,itmpov,info


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  error check:' ,err


  call psi_bld_tmpovrl(temp_ovrlap,desc,info)
  
  if (info == 0) deallocate(prc_v,temp_ovrlap,stat=info)
  if (info /= psb_no_err_) then 
    info=4000
    err=info
    call psb_errpush(err,name)
    Goto 9999
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

  call psb_cd_set_bld(desc,info)


  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cdals
