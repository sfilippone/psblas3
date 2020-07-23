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
module psi_c_collective_mod
  use psi_penv_mod
  use psi_comm_buffers_mod


  interface psb_sum
    module procedure psb_csums, psb_csumv, psb_csumm, &
         & psb_csums_ec, psb_csumv_ec, psb_csumm_ec
  end interface

  interface psb_amx
    module procedure psb_camxs, psb_camxv, psb_camxm, &
         & psb_camxs_ec, psb_camxv_ec, psb_camxm_ec
  end interface

  interface psb_amn
    module procedure psb_camns, psb_camnv, psb_camnm, &
         & psb_camns_ec, psb_camnv_ec, psb_camnm_ec
  end interface


  interface psb_bcast
    module procedure psb_cbcasts, psb_cbcastv, psb_cbcastm, &
         & psb_cbcasts_ec, psb_cbcastv_ec, psb_cbcastm_ec
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_cscan_sums, psb_cscan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_cexscan_sums, psb_cexscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_c_simple_a2av
  end interface psb_simple_a2av

  interface psb_simple_triad_a2av
    module procedure psb_c_e_simple_triad_a2av, psb_c_m_simple_triad_a2av
  end interface psb_simple_triad_a2av


contains 

  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Reduction operations
  !
  ! !!!!!!!!!!!!!!!!!!!!!!



  !
  ! SUM
  !

  subroutine psb_csums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_spk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_spk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_csums

  subroutine psb_csumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_csumv

  subroutine psb_csumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_csumm

  subroutine psb_csums_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csums_ec

  subroutine psb_csumv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csumv_ec
  
  subroutine psb_csumm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csumm_ec


  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_camxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_spk_,mpi_camx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_spk_,mpi_camx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_camxs

  subroutine psb_camxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_camx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_camxv

  subroutine psb_camxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_camx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_camxm


  subroutine psb_camxs_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxs_ec

  subroutine psb_camxv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxv_ec

  subroutine psb_camxm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxm_ec


  !
  ! AMN: Minimum Absolute Value
  !
  
  subroutine psb_camns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_) :: dat_
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
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_spk_,mpi_camn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_spk_,mpi_camn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_camns

  subroutine psb_camnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_camn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_camnv

  subroutine psb_camnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_spk_,mpi_camn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_camn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_camnm


  subroutine psb_camns_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camns_ec

  subroutine psb_camnv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camnv_ec

  subroutine psb_camnm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camnm_ec


  !
  ! BCAST Broadcast
  !
  
  subroutine psb_cbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
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
    call mpi_bcast(dat,1,psb_mpi_c_spk_,root_,ictxt,info)

#endif    
  end subroutine psb_cbcasts

  subroutine psb_cbcastv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
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
    
    call mpi_bcast(dat,size(dat),psb_mpi_c_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_cbcastv

  subroutine psb_cbcastm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
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

    call mpi_bcast(dat,size(dat),psb_mpi_c_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_cbcastm


  subroutine psb_cbcasts_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_cbcasts_ec

  subroutine psb_cbcastv_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_cbcastv_ec

  subroutine psb_cbcastm_ec(ictxt,dat,root)
    implicit none 
    integer(psb_epk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_epk_), intent(in), optional    :: root
    integer(psb_mpk_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_bcast(ictxt_,dat,root_)
    else
      call psb_bcast(ictxt_,dat)
    end if
  end subroutine psb_cbcastm_ec


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_cscan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    complex(psb_spk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: minfo, icomm


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_scan(dat,dat_,1,psb_mpi_c_spk_,mpi_sum,icomm,minfo)
    dat = dat_
#endif    
  end subroutine psb_cscan_sums


  subroutine psb_cexscan_sums(ictxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    complex(psb_spk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: icomm, minfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call mpi_exscan(dat,dat_,1,psb_mpi_c_spk_,mpi_sum,icomm,minfo)
    dat = dat_
#else
    dat = czero
#endif    
  end subroutine psb_cexscan_sums

  subroutine psb_cscan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_scan(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_sum,icomm,minfo)
#endif
  end subroutine psb_cscan_sumv

  subroutine psb_cexscan_sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    icomm = psb_get_mpi_comm(ictxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_exscan(dat,dat_,size(dat),psb_mpi_c_spk_,mpi_sum,icomm,minfo)
#else
    dat = czero
#endif
  end subroutine psb_cexscan_sumv

  subroutine psb_c_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ictxt,info)
    use psi_c_p2p_mod
    implicit none 
    complex(psb_spk_), intent(in)  :: valsnd(:)
    complex(psb_spk_), intent(out) :: valrcv(:)
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

  end subroutine psb_c_simple_a2av

  subroutine psb_c_m_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    complex(psb_spk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    complex(psb_spk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_complex_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_c_spk_,prcid(ip+1),&
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
        p2ptag =  psb_complex_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_c_spk_,prcid(ip+1),&
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

  end subroutine psb_c_m_simple_triad_a2av

  subroutine psb_c_e_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ictxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    complex(psb_spk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    complex(psb_spk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_complex_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_c_spk_,prcid(ip+1),&
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
        p2ptag = psb_complex_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_c_spk_,prcid(ip+1),&
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

  end subroutine psb_c_e_simple_triad_a2av

  
end module psi_c_collective_mod
