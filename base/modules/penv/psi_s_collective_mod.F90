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
module psi_s_collective_mod
  use psi_penv_mod

  interface psb_max
    module procedure psb_smaxs, psb_smaxv, psb_smaxm
  end interface

  interface psb_min
    module procedure psb_smins, psb_sminv, psb_sminm
  end interface psb_min

  interface psb_nrm2
    module procedure psb_s_nrm2s, psb_s_nrm2v
  end interface psb_nrm2

  interface psb_sum
    module procedure psb_ssums, psb_ssumv, psb_ssumm
  end interface

  interface psb_amx
    module procedure psb_samxs, psb_samxv, psb_samxm
  end interface

  interface psb_amn
    module procedure psb_samns, psb_samnv, psb_samnm
  end interface

  interface psb_bcast
    module procedure psb_sbcasts, psb_sbcastv, psb_sbcastm
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_sscan_sums, psb_sscan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_sexscan_sums, psb_sexscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_s_simple_a2av
  end interface psb_simple_a2av

  interface psb_simple_triad_a2av
    module procedure psb_s_e_simple_triad_a2av, psb_s_m_simple_triad_a2av
  end interface psb_simple_triad_a2av

contains 

  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Reduction operations
  !
  ! !!!!!!!!!!!!!!!!!!!!!!


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  MAX
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_smaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_max,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_max,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_smaxs

  subroutine psb_smaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_smaxv

  subroutine psb_smaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
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
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_smaxm

  !
  ! MIN: Minimum Value
  !


  subroutine psb_smins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_min,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_min,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_smins

  subroutine psb_sminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_sminv

  subroutine psb_sminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
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
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_sminm


  ! !!!!!!!!!!!!
  !
  ! Norm 2, only for reals
  !
  ! !!!!!!!!!!!!
  subroutine psb_s_nrm2s(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in) :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional  :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_snrm2_op,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_snrm2_op,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_s_nrm2s

  subroutine psb_s_nrm2v(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,&
           & mpi_snrm2_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,&
             & mpi_snrm2_op,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,&
             & mpi_snrm2_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_s_nrm2v


  !
  ! SUM
  !

  subroutine psb_ssums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_ssums

  subroutine psb_ssumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_ssumv

  subroutine psb_ssumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
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
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_ssumm

  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_samxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samx_op,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samx_op,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_samxs

  subroutine psb_samxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_samxv

  subroutine psb_samxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
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
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_samxm

  !
  ! AMN: Minimum Absolute Value
  !
  
  subroutine psb_samns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samn_op,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samn_op,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_samns

  subroutine psb_samnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
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
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_samnv

  subroutine psb_samnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
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
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_samnm

  !
  ! BCAST Broadcast
  !
  
  subroutine psb_sbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,1,psb_mpi_r_spk_,root_,icomm,info)

#endif    
  end subroutine psb_sbcasts

  subroutine psb_sbcastv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,size(dat),psb_mpi_r_spk_,root_,icomm,info)
#endif    
  end subroutine psb_sbcastv

  subroutine psb_sbcastm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_bcast(dat,size(dat),psb_mpi_r_spk_,root_,icomm,info)
#endif    
  end subroutine psb_sbcastm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_sscan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    real(psb_spk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_scan(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,icomm,minfo)
    dat = dat_
#endif    
  end subroutine psb_sscan_sums


  subroutine psb_sexscan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    real(psb_spk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: icomm, minfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_exscan(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,icomm,minfo)
    dat = dat_
#else
    dat = szero
#endif    
  end subroutine psb_sexscan_sums

  subroutine psb_sscan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_scan(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,icomm,minfo)
#endif
  end subroutine psb_sscan_sumv

  subroutine psb_sexscan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_exscan(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,icomm,minfo)
#else
    dat = szero
#endif
  end subroutine psb_sexscan_sumv

  subroutine psb_s_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ictxt,info)
    use psi_s_p2p_mod
    implicit none 
    real(psb_spk_), intent(in)  :: valsnd(:)
    real(psb_spk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ictxt
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: iam, np, i,j,k, ip, ipx, idx, sz

    call psb_info(ictxt,iam,np)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if

    do ip = 0, np-1
      sz = sdsz(ip+1) 
      if (sz > 0) then
        idx = bsdindx(ip+1)
        call psb_snd(ictxt,valsnd(idx+1:idx+sz),ip) 
      end if
    end do

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        idx = brvindx(ip+1)
        call psb_rcv(ictxt,valrcv(idx+1:idx+sz),ip) 
      end if
    end do

  end subroutine psb_s_simple_a2av

  subroutine psb_s_m_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    real(psb_spk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    real(psb_spk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ictxt
    integer(psb_ipk_), intent(out) :: info

    !Local variables
    integer(psb_ipk_)  :: iam, np, i,j,k, ip, ipx, idx, sz, counter
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret, icomm
    integer(psb_mpk_), allocatable :: prcid(:), rvhd(:,:)

    call psb_info(ictxt,iam,np)

    icomm = psb_get_mpi_comm(ictxt)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if
    allocate(prcid(np),rvhd(np,3))
    prcid = -1

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        prcid(ip+1) = psb_get_mpi_rank(ictxt,ip)
        idx = brvindx(ip+1)
        p2ptag =  psb_real_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_r_spk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,1),iret)
        p2ptag = psb_int_swap_tag
        call mpi_irecv(iarcv(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,2),iret)
        call mpi_irecv(jarcv(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,3),iret)
      end if
    Enddo


    do ip = 0, np-1
      sz = sdsz(ip+1) 
      if (sz > 0) then
        if (prcid(ip+1)<0) prcid(ip+1) = psb_get_mpi_rank(ictxt,ip)
        idx = bsdindx(ip+1)
        p2ptag =  psb_real_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_r_spk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
        p2ptag = psb_int_swap_tag
        call mpi_send(iasnd(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
        call mpi_send(jasnd(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
      end if
    Enddo

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        call mpi_wait(rvhd(ip+1,1),p2pstat,iret)
        call mpi_wait(rvhd(ip+1,2),p2pstat,iret)
        call mpi_wait(rvhd(ip+1,3),p2pstat,iret)
      end if
    Enddo

  end subroutine psb_s_m_simple_triad_a2av

  subroutine psb_s_e_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    real(psb_spk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    real(psb_spk_), intent(out) :: valrcv(:)
    integer(psb_epk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ictxt
    integer(psb_ipk_), intent(out) :: info

    !Local variables
    integer(psb_ipk_)  :: iam, np, i,j,k, ip, ipx, idx, sz, counter
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret, icomm
    integer(psb_mpk_), allocatable :: prcid(:), rvhd(:,:)

    call psb_info(ictxt,iam,np)

    icomm = psb_get_mpi_comm(ictxt)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if
    allocate(prcid(np),rvhd(np,3))
    prcid = -1

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        prcid(ip+1) = psb_get_mpi_rank(ictxt,ip)
        idx = brvindx(ip+1)
        p2ptag =  psb_real_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_r_spk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,1),iret)
        p2ptag = psb_int_swap_tag
        call mpi_irecv(iarcv(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,2),iret)
        call mpi_irecv(jarcv(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
             & p2ptag, icomm,rvhd(ip+1,3),iret)
      end if
    Enddo


    do ip = 0, np-1
      sz = sdsz(ip+1) 
      if (sz > 0) then
        if (prcid(ip+1)<0) prcid(ip+1) = psb_get_mpi_rank(ictxt,ip)
        idx = bsdindx(ip+1)
        p2ptag = psb_real_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_r_spk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
        p2ptag = psb_int_swap_tag
        call mpi_send(iasnd(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
        call mpi_send(jasnd(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
             & p2ptag, icomm,iret)
      end if
    Enddo

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        call mpi_wait(rvhd(ip+1,1),p2pstat,iret)
        call mpi_wait(rvhd(ip+1,2),p2pstat,iret)
        call mpi_wait(rvhd(ip+1,3),p2pstat,iret)
      end if
    Enddo

  end subroutine psb_s_e_simple_triad_a2av

  
end module psi_s_collective_mod
