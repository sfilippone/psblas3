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

module psi_m_p2p_mod
  use psi_penv_mod

  interface psb_snd
    module procedure psb_msnds, psb_msndv, psb_msndm
  end interface

  interface psb_rcv
    module procedure psb_mrcvs, psb_mrcvv, psb_mrcvm
  end interface

contains

  subroutine psb_msnds(ctxt,dat,dst)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(in)  :: dat
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_mpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info 
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    allocate(dat_(1), stat=info)
    dat_(1) = dat
    call psi_snd(ctxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_msnds

  subroutine psb_msndv(ctxt,dat,dst)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(in)  :: dat(:)
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_mpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info

#if defined(SERIAL_MPI) 
#else
    allocate(dat_(size(dat)), stat=info)
    dat_(:) = dat(:)
    call psi_snd(ctxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    

  end subroutine psb_msndv

  subroutine psb_msndm(ctxt,dat,dst,m)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(in)  :: dat(:,:)
    integer(psb_mpk_), intent(in)  :: dst
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_mpk_), allocatable :: dat_(:)
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
    call psi_snd(ctxt,psb_int4_tag,dst,dat_,psb_mesg_queue)
#endif    
  end subroutine psb_msndm

  subroutine psb_mrcvs(ctxt,dat,src)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(out)  :: dat
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_mpk_) :: info, icomm
    integer(psb_mpk_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
    ! do nothing
#else
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_recv(dat,1,psb_mpi_mpk_,src,psb_int4_tag,icomm,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_mrcvs

  subroutine psb_mrcvv(ctxt,dat,src)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(out)  :: dat(:)
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_mpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info, icomm
    integer(psb_mpk_) :: status(mpi_status_size)
#if defined(SERIAL_MPI) 
#else
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_recv(dat,size(dat),psb_mpi_mpk_,src,psb_int4_tag,icomm,status,info)
    call psb_test_nodes(psb_mesg_queue)
#endif    

  end subroutine psb_mrcvv

  subroutine psb_mrcvm(ctxt,dat,src,m)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(out)  :: dat(:,:)
    integer(psb_mpk_), intent(in)  :: src
    integer(psb_ipk_), intent(in), optional :: m
    integer(psb_mpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: info ,m_,n_, ld, mp_rcv_type
    integer(psb_mpk_) :: i,j,k
    integer(psb_mpk_) :: status(mpi_status_size), icomm
#if defined(SERIAL_MPI) 
    ! What should we do here?? 
#else
    if (present(m)) then 
      m_ = m
      ld = size(dat,1)
      n_ = size(dat,2)
      call mpi_type_vector(n_,m_,ld,psb_mpi_mpk_,mp_rcv_type,info)
      if (info == mpi_success) call mpi_type_commit(mp_rcv_type,info)
      icomm = psb_get_mpi_comm(ctxt)
      if (info == mpi_success) call mpi_recv(dat,1,mp_rcv_type,src,&
           & psb_int4_tag,icomm,status,info)
      if (info == mpi_success) call mpi_type_free(mp_rcv_type,info)
    else
      icomm = psb_get_mpi_comm(ctxt) 
      call mpi_recv(dat,size(dat),psb_mpi_mpk_,src,psb_int4_tag,icomm,status,info)
    end if
    if (info /= mpi_success) then 
      write(psb_err_unit,*) 'Error in psb_recv', info
    end if
    call psb_test_nodes(psb_mesg_queue)
#endif    
  end subroutine psb_mrcvm

end module psi_m_p2p_mod
