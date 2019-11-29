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
!
! File: psi_fnd_owner.f90
!
! Subroutine: psi_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
! 
subroutine psi_symm_dep_list(rvsz,adj,idxmap,info,flag)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_indx_map_mod, psb_protect_name => psi_symm_dep_list
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_mpik_), intent(inout)   :: rvsz(0:)
  integer(psb_ipk_), allocatable, intent(inout) :: adj(:)
  class(psb_indx_map), intent(in) :: idxmap
  integer(psb_ipk_), intent(out)     :: info
  integer(psb_ipk_), intent(in), optional :: flag
  
  !
  integer(psb_ipk_), allocatable :: ladj(:) 
  integer(psb_mpik_), allocatable :: sdsz(:)
  integer(psb_mpik_) :: icomm, minfo, iictxt
  integer(psb_ipk_) :: i,n_row,n_col,err_act,hsize,ip,&
       & last_ih, last_j, nidx, nrecv, nadj, flag_
  integer(psb_ipk_) :: mglob, ih
  integer(psb_ipk_) :: ictxt,np,me
  character(len=20)   :: name

  info = psb_success_
  name = 'psi_symm_dep_list'
  call psb_erractionsave(err_act)

  ictxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  mglob   = idxmap%get_gr()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()
  iictxt = ictxt 

  call psb_info(ictxt, me, np)

  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(idxmap%is_valid())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid idxmap')
    goto 9999      
  end if


  nadj = size(adj)
  if (size(rvsz)<np) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid rvsz')
    goto 9999      
  end if

  if (present(flag)) then
    flag_ = flag
  else
    flag_ = psi_symm_flag_norv_
  end if

  select case(flag_)
  case(psi_symm_flag_norv_, psi_symm_flag_inrv_)
    ! Ok, keep going
  case default
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid flag')
    goto 9999      
  end select


  if (flag_ == psi_symm_flag_norv_) then 
  
  ! write(0,*) me,name,' Going through ',nidx,nadj

    Allocate(sdsz(0:np-1), stat=info)
    !
    ! First, send sizes according to adjcncy list
    !
    sdsz = 0 
    do i=1, nadj
      sdsz(adj(i)) = 1
    end do
    !write(0,*)me,' Check on sizes into a2a:',adj(:),nadj,':',sdsz(:)
    
    call mpi_alltoall(sdsz,1,psb_mpi_def_integer,&
         & rvsz,1,psb_mpi_def_integer,icomm,minfo)
  end if

  nrecv = 0 
  do ip=0, np-1
    if (rvsz(ip) > 0) nrecv = nrecv + 1
  end do

  !
  !  Now fix adj to be symmetric
  !
  call psb_realloc(nadj+nrecv,ladj,info)
  ladj(1:nadj) = adj(1:nadj)
  do ip=0, np-1
    if (rvsz(ip)>0) then
      nadj  = nadj + 1 
      ladj(nadj+1:nadj+nrecv) = ip
    end if
  end do
  call psb_msort_unique(ladj,nadj)
  call psb_realloc(nadj,adj,info)
  adj(1:nadj) = ladj(1:nadj) 
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psi_symm_dep_list
