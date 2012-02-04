!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
! File: psb_cdalv.f90
!
! Subroutine: psb_cdalv
!    Allocate descriptor by means of a global map vector V: index I 
!    is assigned to process V(I). It is assumed that V is identical 
!    on all calling processes.
!    
! 
! Arguments: 
!    v       - integer(psb_ipk_), dimension(:).         The array containg the partitioning scheme.
!    ictxt - integer.                         The communication context.
!    desc  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       Return code
!    flag    - integer.                       Are V's contents 0- or 1-based?
subroutine psb_cdalv(v, ictxt, desc, info, flag)
  use psb_base_mod
  use psi_mod
  use psb_repl_map_mod
  use psb_glist_map_mod
  use psb_hash_map_mod
  implicit None
  !....Parameters...
  integer(psb_ipk_), intent(in)               :: ictxt, v(:)
  integer(psb_ipk_), intent(in), optional     :: flag
  integer(psb_ipk_), intent(out)              :: info
  type(psb_desc_type), intent(out)  :: desc

  !locals
  integer(psb_ipk_) :: counter,i,j,np,me,loc_row,err,&
       & loc_col,nprocs,m,n,itmpov, k,glx,&
       & l_ov_ix,l_ov_el,idx, flag_, err_act
  integer(psb_ipk_) :: int_err(5),exch(3)
  integer(psb_ipk_), allocatable  :: temp_ovrlap(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info = psb_success_
  err  = 0
  name = 'psb_cdalv'

  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': ',np,me

  m = size(v)
  n = m
  !... check m and n parameters....
  if (m < 1) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
  else if (n < 1) then
    info = psb_err_iarg_neg_
    int_err(1) = 2
    int_err(2) = n
  else if (size(v)<m) then 
    info = psb_err_iarg_neg_
    int_err(1) = 2
    int_err(2) = size(v)
  endif

  if (info /= psb_success_) then 
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

  call psb_nullify_desc(desc)

  if (present(flag)) then
    flag_=flag
  else
    flag_=0
  endif

  if ((flag_<0).or.(flag_>1)) then
    info = 6
    err=info
    call psb_errpush(info,name)
    goto 9999
  end if

  ! count local rows number
  loc_row = max(1,(m+np-1)/np) 
  ! allocate work vector
  allocate(temp_ovrlap(2),stat=info)
  if (info /= psb_success_) then     
    info=psb_err_alloc_request_
    int_err(1)=2*m+psb_mdata_size_
    call psb_errpush(info,name,i_err=int_err,a_err='integer')
    goto 9999
  endif

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  starting main loop' ,info
  counter = 0
  itmpov  = 0
  temp_ovrlap(:) = -1

  do i=1,m

    if (((v(i)-flag_) > np-1).or.((v(i)-flag_) < 0)) then
      info=psb_err_partfunc_wrong_pid_
      int_err(1)=3
      int_err(2)=v(i) - flag_
      int_err(3)=i
      exit
    end if

    if ((v(i)-flag_) == me) then
      ! this point belongs to me
      counter=counter+1
    end if
  enddo
  loc_row=counter

  !
  ! We have to decide whether we have a "large" index space.
  !
  if (np == 1) then 
    allocate(psb_repl_map :: desc%indxmap, stat=info)
  else
    if (psb_cd_choose_large_state(ictxt,m)) then 
      allocate(psb_hash_map :: desc%indxmap, stat=info)
      if (info == 0) allocate(desc%indxmap%tempvg(m),stat=info)
      if (info ==0) desc%indxmap%tempvg(1:m) = v(1:m) - flag_
    else 
      allocate(psb_glist_map :: desc%indxmap, stat=info)
    end if
  end if


  select type(aa => desc%indxmap) 
  type is (psb_repl_map) 
    call aa%repl_map_init(ictxt,m,info)
  type is (psb_hash_map) 
    call aa%hash_map_init(ictxt,v,info)
  type is (psb_glist_map) 
    call aa%glist_map_init(ictxt,v,info)
  class default 
      ! This cannot happen 
    info = psb_err_internal_error_
    call psb_errpush(info,name)
    Goto 9999
  end select


  call psi_bld_tmpovrl(temp_ovrlap,desc,info)

  deallocate(temp_ovrlap,stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

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

end subroutine psb_cdalv
