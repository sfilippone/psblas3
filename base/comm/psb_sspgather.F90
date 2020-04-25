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
!
! Gathers a sparse matrix onto a single process.
! Two variants:
! 1. Gathers to PSB_s_SPARSE_MAT   (i.e. to matrix with IPK_ indices)
! 2. Gathers to PSB_ls_SPARSE_MAT  (i.e. to matrix with LPK_ indices)
!
! Note: this function uses MPI_ALLGATHERV. At this time, the size of the
! resulting matrix must be within the range of 4 bytes because of the
! restriction on MPI displacements to be 4 bytes. 
! 
!
subroutine  psb_ssp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
#if defined(HAVE_ISO_FORTRAN_ENV)
  use iso_fortran_env
#endif
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
  integer(psb_ipk_) :: nrg, ncg, nzg, nzl
  integer(psb_ipk_) :: err_act, dupl_
  integer(psb_ipk_) :: ip,naggrm1,naggrp1, i, j, k
  logical :: keepnum_, keeploc_
  integer(psb_mpk_) :: ictxt,np,me
  integer(psb_mpk_) :: icomm, minfo, ndx, root_
  integer(psb_mpk_), allocatable :: nzbr(:), idisp(:)
  integer(psb_lpk_), allocatable :: locia(:), locja(:), glbia(:), glbja(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_gather'
  info=psb_success_

  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
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
  if (present(root)) then
    root_ = root
  else
    root_ = -1
  end if
  if ((root_ == -1).or.(root_ == me)) call globa%free()

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
    call psb_realloc(nzl,locia,info)
    call psb_realloc(nzl,locja,info)    
    call psb_loc_to_glob(loc_coo%ia(1:nzl),locia(1:nzl),desc_a,info,iact='I')
    call psb_loc_to_glob(loc_coo%ja(1:nzl),locja(1:nzl),desc_a,info,iact='I')
    nzbr(:) = 0
    nzbr(me+1) = nzl
    call psb_sum(ictxt,nzbr(1:np))
    nzg = sum(nzbr)
    if (nzg <0) then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#if defined(HAVE_ISO_FORTRAN_ENV)
    if (nrg > HUGE(1_psb_mpk_))  then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#endif
    if ((root_ == -1).or.(root_ == me)) then  
      if (info == psb_success_) call psb_realloc(nzg,glbia,info)
      if (info == psb_success_) call psb_realloc(nzg,glbja,info)    
      if (info == psb_success_) call glob_coo%allocate(nrg,ncg,nzg)
    else
      if (info == psb_success_) call psb_realloc(ione,glbia,info)
      if (info == psb_success_) call psb_realloc(ione,glbja,info)    
      if (info == psb_success_) call glob_coo%allocate(ione,ione,ione)
    end if
    
    if (info /= psb_success_) goto 9999

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1)
    if (root_ == -1) then 
      call mpi_allgatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(locia,ndx,psb_mpi_lpk_,&
           & glbia,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(locja,ndx,psb_mpi_lpk_,&
           & glbja,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
    else
      call mpi_gatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(locia,ndx,psb_mpi_lpk_,&
           & glbia,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(locja,ndx,psb_mpi_lpk_,&
           & glbja,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
     
    end if
    
    if (minfo /= psb_success_) then 
      info  = minfo
      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
      goto 9999
    end if    
    call loc_coo%free()
    deallocate(locia,locja, stat=info)
    !
    ! Is the code below safe? For very large cases
    ! the indices in glob_coo will overflow. But then,
    ! for very large cases it does not make sense to
    ! gather the matrix on a single procecss anyway...
    !
    if ((root_ == -1).or.(root_ == me)) then  
      glob_coo%ia(1:nzg) = glbia(1:nzg)
      glob_coo%ja(1:nzg) = glbja(1:nzg)
      call glob_coo%set_nzeros(nzg)
      if (present(dupl)) call glob_coo%set_dupl(dupl)
      call globa%mv_from(glob_coo)
    end if
    deallocate(glbia,glbja, stat=info)

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


subroutine  psb_lssp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
#if defined(HAVE_ISO_FORTRAN_ENV)
  use iso_fortran_env
#endif
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
  type(psb_lsspmat_type), intent(inout) :: globa
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: root, dupl
  logical, intent(in), optional   :: keepnum,keeploc

  type(psb_ls_coo_sparse_mat)     :: loc_coo, glob_coo
  integer(psb_lpk_) :: nrg, ncg, nzg
  integer(psb_ipk_) :: err_act, dupl_
  integer(psb_ipk_) :: ip,naggrm1,naggrp1, i, j, k, nzl
  logical :: keepnum_, keeploc_
  integer(psb_mpk_) :: ictxt,np,me
  integer(psb_mpk_) :: icomm, minfo, ndx, root_
  integer(psb_mpk_), allocatable :: nzbr(:), idisp(:)
  integer(psb_lpk_), allocatable :: lnzbr(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_gather'
  info=psb_success_

  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
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
  if (present(root)) then
    root_ = root
  else
    root_ = -1
  end if

  if ((root_ == -1).or.(root_ == me)) call globa%free()

  if (keepnum_) then 
    nrg = desc_a%get_global_rows()
    ncg = desc_a%get_global_rows()

    allocate(nzbr(np), idisp(np),lnzbr(np),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_; ierr(1) = 3*np
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
    lnzbr = nzbr
    nzg = sum(nzbr)
    if ((nzg < 0).or.(nzg /= sum(lnzbr))) then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#if defined(HAVE_ISO_FORTRAN_ENV)
    if ((nrg > HUGE(1_psb_mpk_)).or.(nzg > HUGE(1_psb_mpk_))&
         & .or.(sum(lnzbr) > HUGE(1_psb_mpk_)))  then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#endif
    if ((root_ == -1).or.(root_ == me)) then  
      if (info == psb_success_) call glob_coo%allocate(nrg,ncg,nzg)
    else
      if (info == psb_success_) call glob_coo%allocate(1_psb_lpk_,1_psb_lpk_,1_psb_lpk_)
    end if
    if (info /= psb_success_) goto 9999
    !
    !  PLS REVIEW AND ADD OVERFLOW ERROR CHECKING
    !

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1)
    if (root_ == -1) then 
      call mpi_allgatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(loc_coo%ia,ndx,psb_mpi_lpk_,&
           & glob_coo%ia,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(loc_coo%ja,ndx,psb_mpi_lpk_,&
           & glob_coo%ja,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
    else
      call mpi_gatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(loc_coo%ia,ndx,psb_mpi_lpk_,&
           & glob_coo%ia,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(loc_coo%ja,ndx,psb_mpi_lpk_,&
           & glob_coo%ja,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
    end if
    if (minfo /= psb_success_) then 
      info  = minfo
      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
      goto 9999
    end if
    call loc_coo%free()
    !
    ! Is the code below safe? For very large cases
    ! the indices in glob_coo will overflow. But then,
    ! for very large cases it does not make sense to
    ! gather the matrix on a single procecss anyway...
    !
    if ((root_ == -1).or.(root_ == me)) then  
      call glob_coo%set_nzeros(nzg)
      if (present(dupl)) call glob_coo%set_dupl(dupl)
      call globa%mv_from(glob_coo)
    end if

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

end subroutine psb_lssp_allgather

subroutine  psb_lslssp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
#if defined(HAVE_ISO_FORTRAN_ENV)
  use iso_fortran_env
#endif
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
  type(psb_lsspmat_type), intent(inout) :: loca
  type(psb_lsspmat_type), intent(inout) :: globa
  type(psb_desc_type), intent(in) :: desc_a
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: root, dupl
  logical, intent(in), optional   :: keepnum,keeploc

  type(psb_ls_coo_sparse_mat)     :: loc_coo, glob_coo
  integer(psb_lpk_) :: nrg, ncg, nzg
  integer(psb_ipk_) :: err_act, dupl_
  integer(psb_lpk_) :: ip,naggrm1,naggrp1, i, j, k, nzl
  logical :: keepnum_, keeploc_
  integer(psb_mpk_) :: ictxt,np,me
  integer(psb_mpk_) :: icomm, minfo, ndx, root_
  integer(psb_mpk_), allocatable :: nzbr(:), idisp(:)
  integer(psb_lpk_), allocatable :: lnzbr(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_gather'
  info=psb_success_

  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
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
  if (present(root)) then
    root_ = root
  else
    root_ = -1
  end if

  if ((root_ == -1).or.(root_ == me)) call globa%free()

  if (keepnum_) then 
    nrg = desc_a%get_global_rows()
    ncg = desc_a%get_global_rows()

    allocate(nzbr(np), idisp(np),lnzbr(np),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_; ierr(1) = 3*np
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
    lnzbr = nzbr
    nzg = sum(nzbr)
    if ((nzg < 0).or.(nzg /= sum(lnzbr))) then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#if defined(HAVE_ISO_FORTRAN_ENV)
    if ((nrg > HUGE(1_psb_mpk_)).or.(nzg > HUGE(1_psb_mpk_))&
         & .or.(sum(lnzbr) > HUGE(1_psb_mpk_)))  then
      info = psb_err_mpi_int_ovflw_
      call psb_errpush(info,name); goto 9999      
    end if
#endif
    if ((root_ == -1).or.(root_ == me)) then  
      if (info == psb_success_) call glob_coo%allocate(nrg,ncg,nzg)
    else
      if (info == psb_success_) call glob_coo%allocate(1_psb_lpk_,1_psb_lpk_,1_psb_lpk_)
    end if
    if (info /= psb_success_) goto 9999

    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1)

    if (root_ == -1) then  
      call mpi_allgatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(loc_coo%ia,ndx,psb_mpi_lpk_,&
           & glob_coo%ia,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_allgatherv(loc_coo%ja,ndx,psb_mpi_lpk_,&
           & glob_coo%ja,nzbr,idisp,&
           & psb_mpi_lpk_,icomm,minfo)
    else
      call mpi_gatherv(loc_coo%val,ndx,psb_mpi_r_spk_,&
           & glob_coo%val,nzbr,idisp,&
           & psb_mpi_r_spk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(loc_coo%ia,ndx,psb_mpi_lpk_,&
           & glob_coo%ia,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
      if (minfo == psb_success_) call &
           & mpi_gatherv(loc_coo%ja,ndx,psb_mpi_lpk_,&
           & glob_coo%ja,nzbr,idisp,&
           & psb_mpi_lpk_,root_,icomm,minfo)
    end if
    if (minfo /= psb_success_) then 
      info  = minfo
      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
      goto 9999
    end if    
    call loc_coo%free()
    !
    if ((root_ == -1).or.(root_ == me)) then  
      call glob_coo%set_nzeros(nzg)
      if (present(dupl)) call glob_coo%set_dupl(dupl)
      call globa%mv_from(glob_coo)
    end if

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
 
end subroutine psb_lslssp_allgather
