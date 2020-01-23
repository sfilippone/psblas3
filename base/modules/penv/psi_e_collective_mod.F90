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
module psi_e_collective_mod
  use psi_penv_mod
  use psi_comm_buffers_mod

  interface psb_max
    module procedure psb_emaxs, psb_emaxv, psb_emaxm, &
         & psb_emaxs_ec, psb_emaxv_ec, psb_emaxm_ec
  end interface

  interface psb_min
    module procedure psb_emins, psb_eminv, psb_eminm, &
         & psb_emins_ec, psb_eminv_ec, psb_eminm_ec
  end interface psb_min


  interface psb_sum
    module procedure psb_esums, psb_esumv, psb_esumm, &
         & psb_esums_ec, psb_esumv_ec, psb_esumm_ec
  end interface

  interface psb_amx
    module procedure psb_eamxs, psb_eamxv, psb_eamxm, &
         & psb_eamxs_ec, psb_eamxv_ec, psb_eamxm_ec
  end interface

  interface psb_amn
    module procedure psb_eamns, psb_eamnv, psb_eamnm, &
         & psb_eamns_ec, psb_eamnv_ec, psb_eamnm_ec
  end interface


  interface psb_bcast
    module procedure psb_ebcasts, psb_ebcastv, psb_ebcastm, &
         & psb_ebcasts_ec, psb_ebcastv_ec, psb_ebcastm_ec
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_escan_sums, psb_escan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_eexscan_sums, psb_eexscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_e_simple_a2av
  end interface psb_simple_a2av

  interface psb_triad_a2av
    module procedure psb_e_e_triad_a2av, psb_e_m_triad_a2av
  end interface psb_triad_a2av


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

  subroutine psb_emaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_epk_,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_epk_,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_emaxs

  subroutine psb_emaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_emaxv

  subroutine psb_emaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_emaxm


  subroutine psb_emaxs_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_emaxs_ec

  subroutine psb_emaxv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_emaxv_ec

  subroutine psb_emaxm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_emaxm_ec


  !
  ! MIN: Minimum Value
  !


  subroutine psb_emins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_epk_,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_epk_,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_emins

  subroutine psb_eminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eminv

  subroutine psb_eminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eminm


  subroutine psb_emins_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_emins_ec

  subroutine psb_eminv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_eminv_ec

  subroutine psb_eminm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_eminm_ec




  !
  ! SUM
  !

  subroutine psb_esums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_epk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_epk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_esums

  subroutine psb_esumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_esumv

  subroutine psb_esumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_esumm

  subroutine psb_esums_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_esums_ec

  subroutine psb_esumv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_esumv_ec
  
  subroutine psb_esumm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_esumm_ec


  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_eamxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_epk_,mpi_eamx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_epk_,mpi_eamx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_eamxs

  subroutine psb_eamxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_eamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eamxv

  subroutine psb_eamxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_eamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eamxm


  subroutine psb_eamxs_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_eamxs_ec

  subroutine psb_eamxv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_eamxv_ec

  subroutine psb_eamxm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_eamxm_ec


  !
  ! AMN: Minimum Absolute Value
  !
  
  subroutine psb_eamns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_epk_,mpi_eamn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_epk_,mpi_eamn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_eamns

  subroutine psb_eamnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_eamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eamnv

  subroutine psb_eamnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_epk_,mpi_eamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_epk_,mpi_eamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_eamnm


  subroutine psb_eamns_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_eamns_ec

  subroutine psb_eamnv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_eamnv_ec

  subroutine psb_eamnm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_eamnm_ec


  !
  ! BCAST Broadcast
  !
  
  subroutine psb_ebcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    call mpi_bcast(dat,1,psb_mpi_epk_,root_,ictxt,info)

#endif    
  end subroutine psb_ebcasts

  subroutine psb_ebcastv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    
    call mpi_bcast(dat,size(dat),psb_mpi_epk_,root_,ictxt,info)
#endif    
  end subroutine psb_ebcastv

  subroutine psb_ebcastm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif

    call mpi_bcast(dat,size(dat),psb_mpi_epk_,root_,ictxt,info)
#endif    
  end subroutine psb_ebcastm


  subroutine psb_ebcasts_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_ebcasts_ec

  subroutine psb_ebcastv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_ebcastv_ec

  subroutine psb_ebcastm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_ebcastm_ec


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_escan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call mpi_scan(dat,dat_,1,psb_mpi_epk_,mpi_sum,ictxt,info)
    dat = dat_
#endif    
  end subroutine psb_escan_sums


  subroutine psb_eexscan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat
    integer(psb_epk_) :: dat_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call mpi_exscan(dat,dat_,1,psb_mpi_epk_,mpi_sum,ictxt,info)
    dat = dat_
#else
    dat = ezero
#endif    
  end subroutine psb_eexscan_sums

  subroutine psb_escan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_realloc(size(dat),dat_,iinfo)
    dat_ = dat
    if (iinfo == psb_success_) &
         & call mpi_scan(dat,dat_,size(dat),psb_mpi_epk_,mpi_sum,ictxt,info)
#endif
  end subroutine psb_escan_sumv

  subroutine psb_eexscan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_epk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_realloc(size(dat),dat_,iinfo)
    dat_ = dat
    if (iinfo == psb_success_) &
         & call mpi_exscan(dat,dat_,size(dat),psb_mpi_epk_,mpi_sum,ictxt,info)
#else
    dat = ezero
#endif
  end subroutine psb_eexscan_sumv

  subroutine psb_e_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ictxt,info)
    use psi_e_p2p_mod
    implicit none 
    integer(psb_epk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    integer(psb_ipk_), intent(in) :: ictxt
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

  end subroutine psb_e_simple_a2av

  subroutine psb_e_m_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_epk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_epk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    integer(psb_ipk_), intent(in) :: ictxt
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
        p2ptag =  psb_int8_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
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
        p2ptag =  psb_int8_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
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

  end subroutine psb_e_m_triad_a2av

  subroutine psb_e_e_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_epk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_epk_), intent(out) :: valrcv(:)
    integer(psb_epk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    integer(psb_ipk_), intent(in) :: ictxt
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
        p2ptag =  psb_int8_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_epk_,prcid(ip+1),&
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
        p2ptag = psb_int8_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_r_dpk_,prcid(ip+1),&
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

  end subroutine psb_e_e_triad_a2av

  
end module psi_e_collective_mod
