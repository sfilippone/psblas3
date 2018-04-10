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

module psi_p2p_mod
  use psi_penv_mod
  use psi_comm_buffers_mod

  interface psb_snd
    module procedure psb_isnds, psb_isndv, psb_isndm, &
         & psb_ssnds, psb_ssndv, psb_ssndm,&
         & psb_dsnds, psb_dsndv, psb_dsndm,&
         & psb_csnds, psb_csndv, psb_csndm,&
         & psb_zsnds, psb_zsndv, psb_zsndm,&
         & psb_lsnds, psb_lsndv, psb_lsndm,&
         & psb_hsnds
  end interface

  interface psb_rcv
    module procedure psb_ircvs, psb_ircvv, psb_ircvm, &
         & psb_srcvs, psb_srcvv, psb_srcvm,&
         & psb_drcvs, psb_drcvv, psb_drcvm,&
         & psb_crcvs, psb_crcvv, psb_crcvm,&
         & psb_zrcvs, psb_zrcvv, psb_zrcvm,&
         & psb_lrcvs, psb_lrcvv, psb_lrcvm,&
         & psb_hrcvs
  end interface


#if defined(LONG_INTEGERS)
  interface psb_snd
    module procedure psb_i4snds, psb_i4sndv, psb_i4sndm
  end interface

  interface psb_rcv
    module procedure psb_i4rcvs, psb_i4rcvv, psb_i4rcvm
  end interface
#endif

#if !defined(LONG_INTEGERS)
  interface psb_snd
    module procedure psb_i8snds, psb_i8sndv, psb_i8sndm
  end interface

  interface psb_rcv
    module procedure psb_i8rcvs, psb_i8rcvv, psb_i8rcvm
  end interface
#endif

#if defined(SHORT_INTEGERS)
  interface psb_snd
    module procedure psb_i2snds, psb_i2sndv, psb_i2sndm
  end interface

  interface psb_rcv
    module procedure psb_i2rcvs, psb_i2rcvv, psb_i2rcvm
  end interface
#endif


#if defined(LONG_INTEGERS)
  interface psb_snd
    module procedure psb_isnds_ic, psb_isndv_ic, psb_isndm_ic, &
         & psb_ssnds_ic, psb_ssndv_ic, psb_ssndm_ic,&
         & psb_dsnds_ic, psb_dsndv_ic, psb_dsndm_ic,&
         & psb_csnds_ic, psb_csndv_ic, psb_csndm_ic,&
         & psb_zsnds_ic, psb_zsndv_ic, psb_zsndm_ic,&
         & psb_lsnds_ic, psb_lsndv_ic, &
         & psb_lsndm_ic, psb_hsnds_ic
  end interface

  interface psb_rcv
    module procedure psb_ircvs_ic, psb_ircvv_ic, psb_ircvm_ic, &
         & psb_srcvs_ic, psb_srcvv_ic, psb_srcvm_ic,&
         & psb_drcvs_ic, psb_drcvv_ic, psb_drcvm_ic,&
         & psb_crcvs_ic, psb_crcvv_ic, psb_crcvm_ic,&
         & psb_zrcvs_ic, psb_zrcvv_ic, psb_zrcvm_ic,&
         & psb_lrcvs_ic, psb_lrcvv_ic, &
         & psb_lrcvm_ic, psb_hrcvs_ic
  end interface

#endif

  integer(psb_mpik_), private, parameter:: psb_int_tag      = 543987
  integer(psb_mpik_), private, parameter:: psb_real_tag     = psb_int_tag      + 1
  integer(psb_mpik_), private, parameter:: psb_double_tag   = psb_real_tag     + 1
  integer(psb_mpik_), private, parameter:: psb_complex_tag  = psb_double_tag   + 1
  integer(psb_mpik_), private, parameter:: psb_dcomplex_tag = psb_complex_tag  + 1
  integer(psb_mpik_), private, parameter:: psb_logical_tag  = psb_dcomplex_tag + 1
  integer(psb_mpik_), private, parameter:: psb_char_tag     = psb_logical_tag  + 1
  integer(psb_mpik_), private, parameter:: psb_int8_tag     = psb_char_tag     + 1
  integer(psb_mpik_), private, parameter:: psb_int2_tag     = psb_int8_tag     + 1
  integer(psb_mpik_), private, parameter:: psb_int4_tag     = psb_int2_tag     + 1

  integer(psb_mpik_),  parameter:: psb_int_swap_tag      = psb_int_tag      + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_real_swap_tag     = psb_real_tag     + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_double_swap_tag   = psb_double_tag   + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_complex_swap_tag  = psb_complex_tag  + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_dcomplex_swap_tag = psb_dcomplex_tag + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_logical_swap_tag  = psb_logical_tag  + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_char_swap_tag     = psb_char_tag     + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_int8_swap_tag     = psb_int8_tag     + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_int2_swap_tag     = psb_int2_tag     + psb_int_tag
  integer(psb_mpik_),  parameter:: psb_int4_swap_tag     = psb_int4_tag     + psb_int_tag


contains


  ! !!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Point-to-point SND
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_isnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_int_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_isnds

  subroutine psb_isndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_int_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_isndv

  subroutine psb_isndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_int_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_isndm

  subroutine psb_ssnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_real_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_ssnds

  subroutine psb_ssndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_real_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_ssndv

  subroutine psb_ssndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: i,j,k,m_,n_
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_real_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_ssndm


  subroutine psb_dsnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_double_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_dsnds

  subroutine psb_dsndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_double_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_dsndv

  subroutine psb_dsndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_double_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_dsndm


  subroutine psb_csnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_complex_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_csnds

  subroutine psb_csndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_complex_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_csndv

  subroutine psb_csndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_complex_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_csndm


  subroutine psb_zsnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_dcomplex_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_zsnds

  subroutine psb_zsndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_dcomplex_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_zsndv

  subroutine psb_zsndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_dcomplex_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_zsndm


  subroutine psb_lsnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    logical, allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_logical_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_lsnds

  subroutine psb_lsndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    logical, allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_logical_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_lsndv

  subroutine psb_lsndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    logical, allocatable :: dat_(:)
    integer(psb_mpik_) :: info
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_logical_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_lsndm


  subroutine psb_hsnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    character(len=*), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    character(len=1), allocatable :: dat_(:)
    integer(psb_mpik_) :: info, l, i
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    l = len(dat) 
    allocate(dat_(l), stat=info)
    do i=1, l
      dat_(i) = dat(i:i)
    end do
    call psi_snd(ictxt,psb_char_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_hsnds

#if defined(LONG_INTEGERS)
  subroutine psb_i4snds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i4snds

  subroutine psb_i4sndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_i4sndv

  subroutine psb_i4sndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_mpik_), intent(in), optional :: m
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i4sndm

#endif


#if !defined(LONG_INTEGERS)
  subroutine psb_i8snds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_int8_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i8snds

  subroutine psb_i8sndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_int8_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_i8sndv

  subroutine psb_i8sndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_int8_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i8sndm

#endif


#if defined(SHORT_INTEGERS)
  subroutine psb_i2snds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(in)  :: dat
    integer(psb_mpik_), intent(in)  :: dst
    integer(2), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_int2_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i2snds

  subroutine psb_i2sndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(in)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(2), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_int2_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_i2sndv

  subroutine psb_i2sndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(in)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    integer(2), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_ipk_) :: i,j,k,m_,n_

#if defined(SERIAL_MPI) 
#else
    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    end if
    n_ = size(dat,2)
    allocate(dat_(m_*n_), stat=info)
    k=1
    do j=1,n_
      do i=1, m_
        dat_(k) = dat(i,j)
        k = k + 1
      end do
    end do
    call psi_snd(ictxt,psb_int2_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_i2sndm

#endif

  ! !!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Point-to-point RCV
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_ircvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_ipk_integer,src,psb_int_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_ircvs

  subroutine psb_ircvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_ipk_integer,src,psb_int_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_ircvv

  subroutine psb_ircvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info, m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_ipk_integer,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_ipk_integer,src,psb_int_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_ircvm


  subroutine psb_srcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_r_spk_,src,psb_real_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_srcvs

  subroutine psb_srcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_r_spk_,src,psb_real_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_srcvv

  subroutine psb_srcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_mpik_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_r_spk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_real_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_r_spk_,src,psb_real_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_srcvm


  subroutine psb_drcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_r_dpk_,src,psb_double_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_drcvs

  subroutine psb_drcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_r_dpk_,src,psb_double_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_drcvv

  subroutine psb_drcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_), intent(in), optional :: m
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_r_dpk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_double_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_r_dpk_,src,&
           & psb_double_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_drcvm


  subroutine psb_crcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_c_spk_,src,psb_complex_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_crcvs

  subroutine psb_crcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_c_spk_,src,psb_complex_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_crcvv

  subroutine psb_crcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_c_spk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_complex_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_c_spk_,src,&
           & psb_complex_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_crcvm


  subroutine psb_zrcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_c_dpk_,src,psb_dcomplex_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_zrcvs

  subroutine psb_zrcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_c_dpk_,src,psb_dcomplex_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_zrcvv

  subroutine psb_zrcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_c_dpk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_dcomplex_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_c_dpk_,src,&
           & psb_dcomplex_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_zrcvm


  subroutine psb_lrcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_logical,src,psb_logical_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_lrcvs

  subroutine psb_lrcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_logical,src,psb_logical_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_lrcvv

  subroutine psb_lrcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    logical, intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_logical,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_logical_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_logical,src,&
           & psb_logical_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_lrcvm


  subroutine psb_hrcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    character(len=*), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    character(len=1), allocatable :: dat_(:)
    integer(psb_mpik_) :: info, l, i
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    l = len(dat) 
    allocate(dat_(l), stat=info)
    call mpi_recv(dat_,l,mpi_character,src,psb_char_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
    do i=1, l
      dat(i:i) = dat_(i) 
    end do
    deallocate(dat_)
#endif    
  end subroutine psb_hrcvs


#if defined(LONG_INTEGERS)

  subroutine psb_i4rcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_def_integer,src,psb_int4_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i4rcvs

  subroutine psb_i4rcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_def_integer,src,psb_int4_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_i4rcvv

  subroutine psb_i4rcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_mpik_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_), intent(in), optional :: m
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_mpik_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_def_integer,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int4_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_def_integer,src,&
           & psb_int4_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i4rcvm

#endif

#if !defined(LONG_INTEGERS)

  subroutine psb_i8rcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_lng_integer,src,psb_int8_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i8rcvs

  subroutine psb_i8rcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_lng_integer,src,psb_int8_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_i8rcvv

  subroutine psb_i8rcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_mpik_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_lng_integer,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int8_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_lng_integer,src,&
           & psb_int8_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i8rcvm

#endif

#if defined(SHORT_INTEGERS)

  subroutine psb_i2rcvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(out)  :: dat
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_def_integer2,src,psb_int2_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i2rcvs

  subroutine psb_i2rcvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(out)  :: dat(:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_) :: info 
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_def_integer2,src,psb_int2_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_i2rcvv

  subroutine psb_i2rcvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(2), intent(out)  :: dat(:,:)
    integer(psb_mpik_), intent(in)  :: src
    integer(psb_mpik_), intent(in), optional :: m
    integer(psb_mpik_) :: info , m_,n_, ld, mp_rcv_type
    integer(psb_ipk_) :: i,j,k
    integer(psb_mpik_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_def_integer2,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int2_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_def_integer2,src,&
           & psb_int2_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i2rcvm

#endif


! 
! Integer * 8 aliases. 
! 

#if defined(LONG_INTEGERS)
  ! !!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Point-to-point SND
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_isnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_isnds_ic

  subroutine psb_isndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_isndv_ic

  subroutine psb_isndm_ic(ictxt,dat,dst,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_isndm_ic

  subroutine psb_ssnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_ssnds_ic

  subroutine psb_ssndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_ssndv_ic

  subroutine psb_ssndm_ic(ictxt,dat,dst,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_ssndm_ic


  subroutine psb_dsnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_dsnds_ic

  subroutine psb_dsndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_dsndv_ic

  subroutine psb_dsndm_ic(ictxt,dat,dst,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_dsndm_ic


  subroutine psb_csnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_csnds_ic

  subroutine psb_csndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_csndv_ic

  subroutine psb_csndm_ic(ictxt,dat,dst,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_csndm_ic


  subroutine psb_zsnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_zsnds_ic

  subroutine psb_zsndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_zsndv_ic

  subroutine psb_zsndm_ic(ictxt,dat,dst,m)
    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_zsndm_ic


  subroutine psb_lsnds_ic(ictxt,dat,dst)
    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_lsnds_ic

  subroutine psb_lsndv_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(in)  :: dat(:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_lsndv_ic

  subroutine psb_lsndm_ic(ictxt,dat,dst,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(in)  :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_lsndm_ic


  subroutine psb_hsnds_ic(ictxt,dat,dst)

    integer(psb_ipk_), intent(in)  :: ictxt
    character(len=*), intent(in)  :: dat
    integer(psb_ipk_), intent(in)  :: dst
    
    integer(psb_mpik_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_hsnds_ic


  ! !!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Point-to-point RCV
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_ircvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ircvs_ic

  subroutine psb_ircvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ircvv_ic

  subroutine psb_ircvm_ic(ictxt,dat,src,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ircvm_ic

  subroutine psb_srcvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_srcvs_ic

  subroutine psb_srcvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_srcvv_ic

  subroutine psb_srcvm_ic(ictxt,dat,src,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_spk_), intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_srcvm_ic


  subroutine psb_drcvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_drcvs_ic

  subroutine psb_drcvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_drcvv_ic

  subroutine psb_drcvm_ic(ictxt,dat,src,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    real(psb_dpk_), intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_drcvm_ic


  subroutine psb_crcvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_crcvs_ic

  subroutine psb_crcvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_crcvv_ic

  subroutine psb_crcvm_ic(ictxt,dat,src,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_spk_), intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_crcvm_ic


  subroutine psb_zrcvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_zrcvs_ic

  subroutine psb_zrcvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_zrcvv_ic

  subroutine psb_zrcvm_ic(ictxt,dat,src,m)
    integer(psb_ipk_), intent(in)  :: ictxt
    complex(psb_dpk_), intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_zrcvm_ic


  subroutine psb_lrcvs_ic(ictxt,dat,src)
    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_lrcvs_ic

  subroutine psb_lrcvv_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(out) :: dat(:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_lrcvv_ic

  subroutine psb_lrcvm_ic(ictxt,dat,src,m)

    integer(psb_ipk_), intent(in)  :: ictxt
    logical, intent(out) :: dat(:,:)
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_lrcvm_ic


  subroutine psb_hrcvs_ic(ictxt,dat,src)

    integer(psb_ipk_), intent(in)  :: ictxt
    character(len=*), intent(out) :: dat
    integer(psb_ipk_), intent(in)  :: src
    
    integer(psb_mpik_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_hrcvs_ic


#endif


end module psi_p2p_mod
