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
! File: psb_ialloc.f90
!
! Function: psb_ialloc
!    Allocates dense integer matrix for PSBLAS routines
!    The descriptor may be in either the build or assembled state.
! 
! Arguments: 
!    x      - the matrix to be allocated.
!    desc_a - the communication descriptor.
!    info   - possibly returns an error code
!    n      - optional number of columns.
!    lb     - optional lower bound on column indices
subroutine psb_ialloc(x, desc_a, info, n, lb)
  use psb_base_mod, psb_protect_name => psb_ialloc
  implicit none
  
  !....parameters...
  integer(psb_ipk_), allocatable, intent(out) :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), optional, intent(in)     :: n, lb

  !locals
  integer(psb_ipk_) :: np,me,err,nr,i,j,err_act
  integer(psb_ipk_) :: ictxt,n_
  integer(psb_ipk_) :: int_err(5), exch(3)
  character(len=20)   :: name

  name='psb_geall'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  err=0
  int_err(1)=0
  call psb_erractionsave(err_act)
  
  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check m and n parameters....
  if (.not.psb_is_ok_desc(desc_a)) then 
    info = psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(n)) then 
    n_ = n
  else
    n_ = 1
  endif
  !global check on n parameters
  if (me == psb_root_) then
    exch(1)=n_
    call psb_bcast(ictxt,exch(1),root=psb_root_)
  else
    call psb_bcast(ictxt,exch(1),root=psb_root_)
    if (exch(1) /= n_) then
      info=psb_err_parm_differs_among_procs_
      int_err(1)=1
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
  endif

  !....allocate x .....
  if (psb_is_asb_desc(desc_a).or.psb_is_upd_desc(desc_a)) then
    nr = max(1,desc_a%get_local_cols())
  else if (psb_is_bld_desc(desc_a)) then
    nr = max(1,desc_a%get_local_rows())
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,int_err,a_err='Invalid desc_a')
    goto 9999
  endif
  
  call psb_realloc(nr,n_,x,info,lb2=lb)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    int_err(1)=nr*n_
    call psb_errpush(info,name,int_err,a_err='integer')
    goto 9999
  endif
  
  x(:,:) = izero

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_ialloc



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
! Function: psb_iallocv
!    Allocates dense matrix for PSBLAS routines
!    The descriptor may be in either the build or assembled state.
! 
! Arguments: 
!    x(:)   - the matrix to be allocated.
!    desc_a - the communication descriptor.
!    info   - return code
subroutine psb_iallocv(x, desc_a, info,n)
  use psb_base_mod, psb_protect_name => psb_iallocv
  implicit none

  !....parameters...
  integer(psb_ipk_), allocatable, intent(out) :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), optional, intent(in)     :: n

  !locals
  integer(psb_ipk_) :: np,me,nr,i,err_act
  integer(psb_ipk_) :: ictxt, int_err(5)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_geall'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check m and n parameters....
  if (.not.psb_is_ok_desc(desc_a)) then 
    info = psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! As this is a rank-1 array, optional parameter N is actually ignored.

  !....allocate x .....
  if (psb_is_asb_desc(desc_a).or.psb_is_upd_desc(desc_a)) then
    nr = max(1,desc_a%get_local_cols())
  else if (psb_is_bld_desc(desc_a)) then
    nr = max(1,desc_a%get_local_rows())
  else
    info = psb_err_internal_error_
    call psb_errpush(info,name,int_err,a_err='Invalid desc_a')
    goto 9999
  endif
  
  call psb_realloc(nr,x,info)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    int_err(1)=nr
    call psb_errpush(info,name,int_err,a_err='integer')
    goto 9999
  endif
  
  x(:) = izero

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_iallocv

