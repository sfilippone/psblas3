
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



  integer, private, parameter:: psb_int_tag      = 543987
  integer, private, parameter:: psb_real_tag     = psb_int_tag      + 1
  integer, private, parameter:: psb_double_tag   = psb_real_tag     + 1
  integer, private, parameter:: psb_complex_tag  = psb_double_tag   + 1
  integer, private, parameter:: psb_dcomplex_tag = psb_complex_tag  + 1
  integer, private, parameter:: psb_logical_tag  = psb_dcomplex_tag + 1
  integer, private, parameter:: psb_char_tag     = psb_logical_tag  + 1
  integer, private, parameter:: psb_int8_tag     = psb_char_tag     + 1
  integer, private, parameter:: psb_int2_tag     = psb_int8_tag     + 1

  integer,  parameter:: psb_int_swap_tag      = psb_int_tag      + psb_int_tag
  integer,  parameter:: psb_real_swap_tag     = psb_real_tag     + psb_int_tag
  integer,  parameter:: psb_double_swap_tag   = psb_double_tag   + psb_int_tag
  integer,  parameter:: psb_complex_swap_tag  = psb_complex_tag  + psb_int_tag
  integer,  parameter:: psb_dcomplex_swap_tag = psb_dcomplex_tag + psb_int_tag
  integer,  parameter:: psb_logical_swap_tag  = psb_logical_tag  + psb_int_tag
  integer,  parameter:: psb_char_swap_tag     = psb_char_tag     + psb_int_tag
  integer,  parameter:: psb_int8_swap_tag     = psb_int8_tag     + psb_int_tag
  integer,  parameter:: psb_int2_swap_tag     = psb_int2_tag     + psb_int_tag


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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(in)  :: dat
    integer, intent(in)  :: dst
    integer(psb_int_k_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    integer(psb_int_k_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    integer(psb_int_k_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    real(psb_spk_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    real(psb_spk_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    real(psb_spk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    complex(psb_spk_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    complex(psb_spk_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    complex(psb_spk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    complex(psb_dpk_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    complex(psb_dpk_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    complex(psb_dpk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    logical, intent(in)  :: dat
    integer, intent(in)  :: dst
    logical, allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    logical, intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    logical, allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    logical, intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    logical, allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    character(len=*), intent(in)  :: dat
    integer, intent(in)  :: dst
    character(len=1), allocatable :: dat_(:)
    integer :: info, l, i
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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat
    integer, intent(in)  :: dst
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    integer(2), intent(in)  :: dat
    integer, intent(in)  :: dst
    integer(2), allocatable :: dat_(:)
    integer :: info 
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
    integer, intent(in)  :: ictxt
    integer(2), intent(in)  :: dat(:)
    integer, intent(in)  :: dst
    integer(2), allocatable :: dat_(:)
    integer :: info 

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
    integer, intent(in)  :: ictxt
    integer(2), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m
    integer(2), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_

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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_integer,src,psb_int_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    integer(psb_int_k_), allocatable :: dat_(:)
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_integer,src,psb_int_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(psb_int_k_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    integer(psb_int_k_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_integer,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_integer,src,psb_int_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_real,src,psb_real_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    real(psb_spk_), allocatable :: dat_(:)
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_real,src,psb_real_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    real(psb_spk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_real,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_real_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_real,src,psb_real_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_double_precision,src,psb_double_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_double_precision,src,psb_double_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_double_precision,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_double_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_double_precision,src,&
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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_complex,src,psb_complex_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    complex(psb_spk_), allocatable :: dat_(:)
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_complex,src,psb_complex_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    complex(psb_spk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_complex,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_complex_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_complex,src,&
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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_double_complex,src,psb_dcomplex_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    complex(psb_dpk_), allocatable :: dat_(:)
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_double_complex,src,psb_dcomplex_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    complex(psb_dpk_), allocatable :: dat_(:)
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_double_complex,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_dcomplex_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_double_complex,src,&
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
    integer, intent(in)  :: ictxt
    logical, intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
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
    integer, intent(in)  :: ictxt
    logical, intent(out)  :: dat(:)
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
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
    integer, intent(in)  :: ictxt
    logical, intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
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
    integer, intent(in)  :: ictxt
    character(len=*), intent(out)  :: dat
    integer, intent(in)  :: src
    character(len=1), allocatable :: dat_(:)
    integer :: info, l, i
    integer :: status(mpi_status_size)
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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_integer8,src,psb_int8_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_integer8,src,psb_int8_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(psb_long_int_k_), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_integer8,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int8_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_integer8,src,&
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
    integer, intent(in)  :: ictxt
    integer(2), intent(out)  :: dat
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,mpi_integer2,src,psb_int2_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(2), intent(out)  :: dat(:)
    integer, intent(in)  :: src
    integer :: info 
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),mpi_integer2,src,psb_int2_tag,ictxt,status,info)
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
    integer, intent(in)  :: ictxt
    integer(2), intent(out)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
    integer :: info ,i,j,k,m_,n_, ld, mp_rcv_type
    integer :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,mpi_integer2,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int2_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),mpi_integer2,src,&
           & psb_int2_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_i2rcvm

#endif

end module psi_p2p_mod
