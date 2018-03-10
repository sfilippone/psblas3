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

module psi_e_p2p_mod
  use psi_penv_mod
  use psi_comm_buffers_mod

  interface psb_snd
    module procedure psb_esnds, psb_esndv, psb_esndm, &
         & psb_esnds_ec, psb_esndv_ec, psb_esndm_ec 
  end interface

  interface psb_rcv
    module procedure psb_ercvs, psb_ercvv, psb_ercvm, &
         & psb_ercvs_ec, psb_ercvv_ec, psb_ercvm_ec 
  end interface

contains

  subroutine psb_esnds(ictxt,dat,dst)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ictxt,psb_int8_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_esnds

  subroutine psb_esndv(ictxt,dat,dst)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat(:)
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info 

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ictxt,psb_int8_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_esndv

  subroutine psb_esndm(ictxt,dat,dst,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat(:,:)
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: i,j,k,m_,n_
    integer(psb_mpk_) :: info 

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
  end subroutine psb_esndm

  subroutine psb_ercvs(ictxt,dat,src)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out)  :: dat
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_mpk_) :: info 
    integer(psb_mpk_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    call mpi_recv(dat,1,psb_mpi_epk_,src,psb_int8_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_ercvs

  subroutine psb_ercvv(ictxt,dat,src)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out)  :: dat(:)
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info 
    integer(psb_mpk_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    call mpi_recv(dat,size(dat),psb_mpi_epk_,src,psb_int8_tag,ictxt,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_ercvv

  subroutine psb_ercvm(ictxt,dat,src,m)
    use psi_comm_buffers_mod 

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out)  :: dat(:,:)
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_epk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_mpk_) :: i,j,k
    integer(psb_mpk_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_epk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int8_tag,ictxt,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      call mpi_recv(dat,size(dat),psb_mpi_epk_,src,psb_int8_tag,ictxt,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_ercvm


  subroutine psb_esnds_ec(ictxt,dat,dst)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat
    integer(psb_epk_), intent(in)  :: dst
    
    integer(psb_mpk_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_esnds_ec

  subroutine psb_esndv_ec(ictxt,dat,dst)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat(:)
    integer(psb_epk_), intent(in)  :: dst
    
    integer(psb_mpk_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_esndv_ec

  subroutine psb_esndm_ec(ictxt,dat,dst,m)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(in)  :: dat(:,:)
    integer(psb_epk_), intent(in)  :: dst
    
    integer(psb_mpk_) :: iictxt, idst 

    iictxt = ictxt
    idst   = dst 
    call psb_snd(iictxt, dat, idst)

  end subroutine psb_esndm_ec

  subroutine psb_ercvs_ec(ictxt,dat,src)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out) :: dat
    integer(psb_epk_), intent(in)  :: src
    
    integer(psb_mpk_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ercvs_ec

  subroutine psb_ercvv_ec(ictxt,dat,src)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out) :: dat(:)
    integer(psb_epk_), intent(in)  :: src
    
    integer(psb_mpk_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ercvv_ec

  subroutine psb_ercvm_ec(ictxt,dat,src,m)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out) :: dat(:,:)
    integer(psb_epk_), intent(in)  :: src
    
    integer(psb_mpk_) :: iictxt, isrc 

    iictxt = ictxt
    isrc   = src 
    call psb_rcv(iictxt, dat, isrc)

  end subroutine psb_ercvm_ec


end module psi_e_p2p_mod
