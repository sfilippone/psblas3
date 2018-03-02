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
! File:  psb_sspgather.f90
subroutine  psb_ssp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
  use psb_desc_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_mat_mod
  use psb_tools_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  type(psb_sspmat_type), intent(inout) :: loca
  type(psb_sspmat_type), intent(inout) :: globa
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: root, dupl
  logical, intent(in), optional   :: keepnum,keeploc

  type(psb_s_coo_sparse_mat)      :: loc_coo, glob_coo
  integer(psb_ipk_) :: err_act, dupl_, nrg, ncg, nzg
  integer(psb_ipk_) :: ip,naggrm1,naggrp1, i, j, k, nzl
  logical :: keepnum_, keeploc_
  integer(psb_mpk_) :: ictxt,np,me
  integer(psb_mpk_) :: icomm, minfo, ndx
  integer(psb_mpk_), allocatable :: nzbr(:), idisp(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_gather'
  if (psb_get_errstatus().ne.0) return 
  info=psb_success_

  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  icomm = desc_a%get_mpic()
  call psb_info(ictxt, me, np)

  if (present(keepnum)) then 
    keepnum_ = keepnum
  else
    keepnum_ = .true.
  end if
  if (present(keeploc)) then 
    keeploc_ = keeploc
  else
    keeploc_ = .true.
  end if
  call globa%free()

  if (keepnum_) then 
    nrg = desc_a%get_global_rows()
    ncg = desc_a%get_global_rows()

    allocate(nzbr(np), idisp(np),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      ierr(1) = 2*np
      call psb_errpush(info,name,i_err=ierr,a_err='integer')
      goto 9999      
    end if

    
    if (keeploc_) then
      call loca%cp_to(loc_coo)
    else
      call loca%mv_to(loc_coo)
    end if
    nzl = loc_coo%get_nzeros()
    call psb_loc_to_glob(loc_coo%ia(1:nzl),desc_a,info,iact='I')
    call psb_loc_to_glob(loc_coo%ja(1:nzl),desc_a,info,iact='I')
    nzbr(:) = 0
    nzbr(me+1) = nzl
    call psb_sum(ictxt,nzbr(1:np))
    nzg = sum(nzbr)
    if (info == psb_success_) call glob_coo%allocate(nrg,ncg,nzg)
    if (info /= psb_success_) goto 9999
    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1) 
    call mpi_allgatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
         & glob_coo%val,nzbr,idisp,&
         & psb_mpi_r_spk_,icomm,minfo)
    if (minfo == psb_success_) call &
         & mpi_allgatherv(loc_coo%ia,ndx,psb_mpi_ipk_int,&
         & glob_coo%ia,nzbr,idisp,&
         & psb_mpi_ipk_int,icomm,minfo)
    if (minfo == psb_success_) call &
         & mpi_allgatherv(loc_coo%ja,ndx,psb_mpi_ipk_int,&
         & glob_coo%ja,nzbr,idisp,&
         & psb_mpi_ipk_int,icomm,minfo)
    
    if (minfo /= psb_success_) then 
      info  = minfo
      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
      goto 9999
    end if
    
    call loc_coo%free()
    call glob_coo%set_nzeros(nzg)
    if (present(dupl)) call glob_coo%set_dupl(dupl)
    call globa%mv_from(glob_coo)

  else
    write(psb_err_unit,*) 'SP_ALLGATHER: Not implemented yet with keepnum ',keepnum_
    info = -1
    goto 9999
  end if



  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(ione*ictxt,err_act)

  return
 
end subroutine psb_ssp_allgather
