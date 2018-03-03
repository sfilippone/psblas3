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
module psi_m_reduce_mod
  use psi_penv_mod

  interface psb_max
    module procedure psb_mmaxs, psb_mmaxv, psb_mmaxm, &
         & psb_mmaxs_ec, psb_mmaxv_ec, psb_mmaxm_ec
  end interface

  interface psb_min
    module procedure psb_mmins, psb_mminv, psb_mminm, &
         & psb_mmins_ec, psb_mminv_ec, psb_mminm_ec
  end interface psb_min


  interface psb_sum
    module procedure psb_msums, psb_msumv, psb_msumm, &
         & psb_msums_ec, psb_msumv_ec, psb_msumm_ec
  end interface

  interface psb_amx
    module procedure psb_mamxs, psb_mamxv, psb_mamxm, &
         & psb_mamxs_ec, psb_mamxv_ec, psb_mamxm_ec
  end interface

  interface psb_amn
    module procedure psb_mamns, psb_mamnv, psb_mamnm, &
         & psb_mamns_ec, psb_mamnv_ec, psb_mamnm_ec
  end interface



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

  subroutine psb_mmaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_mpk_,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_mpk_,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_mmaxs

  subroutine psb_mmaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mmaxv

  subroutine psb_mmaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mmaxm


  subroutine psb_mmaxs_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_mmaxs_ec

  subroutine psb_mmaxv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_mmaxv_ec

  subroutine psb_mmaxm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_mmaxm_ec


  !
  ! MIN: Minimum Value
  !


  subroutine psb_mmins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_mpk_,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_mpk_,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_mmins

  subroutine psb_mminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mminv

  subroutine psb_mminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mminm


  subroutine psb_mmins_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_mmins_ec

  subroutine psb_mminv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_mminv_ec

  subroutine psb_mminm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_mminm_ec




  !
  ! SUM
  !

  subroutine psb_msums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_mpk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_mpk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_msums

  subroutine psb_msumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_msumv

  subroutine psb_msumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_msumm

  subroutine psb_msums_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_msums_ec

  subroutine psb_msumv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_msumv_ec
  
  subroutine psb_msumm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_msumm_ec


  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_mamxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_mpk_,mpi_mamx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_mpk_,mpi_mamx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_mamxs

  subroutine psb_mamxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_mamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mamxv

  subroutine psb_mamxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_mamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mamxm


  subroutine psb_mamxs_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_mamxs_ec

  subroutine psb_mamxv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_mamxv_ec

  subroutine psb_mamxm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_mamxm_ec


  !
  ! AMN: Minimum Absolute Value
  !
  
  subroutine psb_mamns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_mpk_,mpi_mamn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_mpk_,mpi_mamn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_mamns

  subroutine psb_mamnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_mamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mamnv

  subroutine psb_mamnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_mpk_,mpi_mamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_mpk_,mpi_mamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_mamnm


  subroutine psb_mamns_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_mamns_ec

  subroutine psb_mamnv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_mamnv_ec

  subroutine psb_mamnm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_mamnm_ec

end module psi_m_reduce_mod
