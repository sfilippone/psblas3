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
module psi_collective_mod
  use psi_penv_mod
  use psi_m_collective_mod
  use psi_e_collective_mod
  use psi_s_collective_mod
  use psi_d_collective_mod
  use psi_c_collective_mod
  use psi_z_collective_mod

  interface psb_bcast
    module procedure  psb_hbcasts, psb_hbcastv,&
         & psb_lbcasts, psb_lbcastv
  end interface psb_bcast


#if defined(SHORT_INTEGERS)
  interface psb_sum
    module procedure psb_i2sums, psb_i2sumv, psb_i2summ
  end interface psb_sum
#endif

contains

  subroutine psb_hbcasts(ictxt,dat,root,length)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer(psb_mpk_), intent(in), optional   :: root,length

    integer(psb_mpk_) :: iam, np, root_,length_,info, icomm

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    if (present(length)) then
      length_ = length
    else
      length_ = len(dat)
    endif

    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,length_,MPI_CHARACTER,root_,icomm,info)
#endif

  end subroutine psb_hbcasts

  subroutine psb_hbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_, icomm
    integer(psb_mpk_) :: length_,info, size_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ =  psb_root_
    endif
    length_ = len(dat)
    size_   = size(dat)

    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,length_*size_,MPI_CHARACTER,root_,icomm,info)
#endif

  end subroutine psb_hbcastv

  subroutine psb_lbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_,info, icomm

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,1,MPI_LOGICAL,root_,icomm,info)
#endif

  end subroutine psb_lbcasts

  subroutine psb_lallreduceand(ictxt,dat,rec)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    logical, intent(inout), optional :: rec

    integer(psb_mpk_) :: iam, np, info, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    if (present(rec)) then
      call mpi_allreduce(dat,rec,1,MPI_LOGICAL,MPI_LAND,icomm,info)
    else
      call mpi_allreduce(MPI_IN_PLACE,dat,1,MPI_LOGICAL,MPI_LAND,icomm,info)
    endif
#endif

end subroutine psb_lallreduceand


  subroutine psb_lbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)             :: ictxt
    logical, intent(inout)          :: dat(:)
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_,info, icomm

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,size(dat),MPI_LOGICAL,root_,icomm,info)
#endif

  end subroutine psb_lbcastv

#if defined(SHORT_INTEGERS)
  subroutine psb_i2sums(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_i2pk_) :: dat_
    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ictxt)
    if (root_ == -1) then
      call mpi_allreduce(dat,dat_,1,psb_mpi_i2pk_,mpi_sum,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif

#endif
  end subroutine psb_i2sums

  subroutine psb_i2sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_i2pk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ictxt)
    if (root_ == -1) then
      call psb_realloc(size(dat),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_i2pk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif
  end subroutine psb_i2sumv

  subroutine psb_i2summ(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_i2pk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ictxt)
    if (root_ == -1) then
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_i2pk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif
  end subroutine psb_i2summ

#endif

end module psi_collective_mod
