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
!
! File: psi_fnd_owner.f90
!
! Subroutine: psi_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    nv       - integer                 Number of indices required on  the calling
!                                       process 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       Note: the indices should be unique!
!    iprc(:)  - integer(psb_ipk_), allocatable    Output: process identifiers for
!                                       the corresponding
!                                       indices
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                return code.
! 
subroutine psi_fnd_owner(nv,idx,iprc,desc,info)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psi_mod, psb_protect_name => psi_fnd_owner
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_ipk_), intent(in) :: nv
  integer(psb_ipk_), intent(in) :: idx(:)
  integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out) :: info


  integer(psb_ipk_), allocatable :: hsz(:),hidx(:),helem(:),hproc(:),&
       & sdsz(:),sdidx(:), rvsz(:), rvidx(:),answers(:,:),idxsrch(:,:)

  integer(psb_ipk_) :: i,n_row,n_col,err_act,ih,icomm,hsize,ip,isz,k,j,&
       & last_ih, last_j
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter  :: gettime=.false.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = desc%get_context()
  icomm   = desc%get_mpic()
  n_row   = desc%get_local_rows()
  n_col   = desc%get_local_cols()


  ! check on blacs grid 
  call psb_info(ictxt, me, np)

  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (nv < 0 ) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(desc%is_ok())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid desc')
    goto 9999      
  end if

  call desc%fnd_owner(idx(1:nv),iprc,info)
  
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='desc%fnd_owner') 
    goto 9999      
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psi_fnd_owner
