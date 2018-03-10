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


module psi_bcast_mod
  use psb_const_mod
  use psi_penv_mod
  interface psb_bcast
    module procedure psb_ibcasts, psb_ibcastv, psb_ibcastm,&
         & psb_dbcasts, psb_dbcastv, psb_dbcastm,&
         & psb_zbcasts, psb_zbcastv, psb_zbcastm,&
         & psb_sbcasts, psb_sbcastv, psb_sbcastm,&
         & psb_cbcasts, psb_cbcastv, psb_cbcastm,&
         & psb_hbcasts, psb_hbcastv,&
         & psb_lbcasts, psb_lbcastv
  end interface psb_bcast

#if defined(LONG_INTEGERS)
  interface psb_bcast
    module procedure psb_ibcasts_ic, psb_ibcastv_ic, psb_ibcastm_ic,&
         & psb_dbcasts_ic, psb_dbcastv_ic, psb_dbcastm_ic,&
         & psb_zbcasts_ic, psb_zbcastv_ic, psb_zbcastm_ic,&
         & psb_sbcasts_ic, psb_sbcastv_ic, psb_sbcastm_ic,&
         & psb_cbcasts_ic, psb_cbcastv_ic, psb_cbcastm_ic,&
         & psb_hbcasts_ic, psb_hbcastv_ic, &
         & psb_lbcasts_ic, psb_lbcastv_ic
  end interface psb_bcast
#else 
  interface psb_bcast
    module procedure psb_i8bcasts, psb_i8bcastv, psb_i8bcastm
  end interface psb_bcast
#endif

contains

  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Broadcasts
  !
  ! !!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_ibcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    integer(psb_ipk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_ipk_,root_,ictxt,info)
#endif    
  end subroutine psb_ibcasts

  subroutine psb_ibcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional  :: root

    integer(psb_mpk_) :: iam, np, root_,  info
#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_ipk_,root_,ictxt,info)
#endif    
  end subroutine psb_ibcastv

  subroutine psb_ibcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_ipk_,root_,ictxt,info)
#endif    
  end subroutine psb_ibcastm


  subroutine psb_sbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_r_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_sbcasts


  subroutine psb_sbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_r_spk_,root_,ictxt,info)

#endif    
  end subroutine psb_sbcastv

  subroutine psb_sbcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_r_spk_,root_,ictxt,info)

#endif    
  end subroutine psb_sbcastm


  subroutine psb_dbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_r_dpk_,root_,ictxt,info)
#endif    
  end subroutine psb_dbcasts


  subroutine psb_dbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_r_dpk_,root_,ictxt,info)
#endif    
  end subroutine psb_dbcastv

  subroutine psb_dbcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_r_dpk_,root_,ictxt,info)
#endif    
  end subroutine psb_dbcastm

  subroutine psb_cbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_c_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_cbcasts

  subroutine psb_cbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_c_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_cbcastv

  subroutine psb_cbcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_c_spk_,root_,ictxt,info)
#endif    
  end subroutine psb_cbcastm

  subroutine psb_zbcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_c_dpk_,root_,ictxt,info)
#endif    
  end subroutine psb_zbcasts

  subroutine psb_zbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_c_dpk_,root_,ictxt,info) 
#endif    
  end subroutine psb_zbcastv

  subroutine psb_zbcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_c_dpk_,root_,ictxt,info)
#endif    
  end subroutine psb_zbcastm


  subroutine psb_hbcasts(ictxt,dat,root,length)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer(psb_mpk_), intent(in), optional   :: root,length

    integer(psb_mpk_) :: iam, np, root_,length_,info

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

    call mpi_bcast(dat,length_,MPI_CHARACTER,root_,ictxt,info)
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
    integer(psb_mpk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_,length_,info, size_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ =  psb_root_
    endif
    length_ = len(dat)
    size_   = size(dat) 

    call psb_info(ictxt,iam,np)

    call mpi_bcast(dat,length_*size_,MPI_CHARACTER,root_,ictxt,info)
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
    integer(psb_mpk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_,info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,MPI_LOGICAL,root_,ictxt,info)
#endif    

  end subroutine psb_lbcasts


  subroutine psb_lbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat(:)
    integer(psb_mpk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iam, np, root_,info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),MPI_LOGICAL,root_,ictxt,info)
#endif    

  end subroutine psb_lbcastv


#if !defined(LONG_INTEGERS)

  subroutine psb_i8bcasts(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)      :: ictxt
    integer(psb_epk_), intent(inout)   :: dat
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_lpk_,root_,ictxt,info)
#endif    
  end subroutine psb_i8bcasts

  subroutine psb_i8bcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    integer(psb_epk_), intent(inout) :: dat(:)
    integer(psb_mpk_), intent(in), optional  :: root

    integer(psb_mpk_) :: iam, np, root_,  info
#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_lpk_,root_,ictxt,info)
#endif    
  end subroutine psb_i8bcastv

  subroutine psb_i8bcastm(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)    :: ictxt
    integer(psb_epk_), intent(inout) :: dat(:,:)
    integer(psb_mpk_), intent(in), optional :: root

    integer(psb_mpk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_lpk_,root_,ictxt,info)
#endif    
  end subroutine psb_i8bcastm

#endif


#if defined(LONG_INTEGERS)

  subroutine psb_ibcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)      :: ictxt
    integer(psb_ipk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_ibcasts_ic

  subroutine psb_ibcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional  :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_ibcastv_ic

  subroutine psb_ibcastm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_ibcastm_ic


  subroutine psb_sbcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)      :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_sbcasts_ic


  subroutine psb_sbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_sbcastv_ic

  subroutine psb_sbcastm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_sbcastm_ic


  subroutine psb_dbcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)      :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_dbcasts_ic


  subroutine psb_dbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_dbcastv_ic

  subroutine psb_dbcastm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_dbcastm_ic

  subroutine psb_cbcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)      :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_cbcasts_ic

  subroutine psb_cbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_cbcastv_ic

  subroutine psb_cbcastm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_cbcastm_ic

  subroutine psb_zbcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)      :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_zbcasts_ic

  subroutine psb_zbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_zbcastv_ic

  subroutine psb_zbcastm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_zbcastm_ic


  subroutine psb_hbcasts_ic(ictxt,dat,root,length)
    implicit none 
    integer(psb_ipk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer(psb_ipk_), intent(in), optional   :: root,length

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_hbcasts_ic

  subroutine psb_hbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_hbcastv_ic

  subroutine psb_lbcasts_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_lbcasts_ic


  subroutine psb_lbcastv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat(:)
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_mpk_) :: iictxt, root_

    iictxt = ictxt 
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif
    call psb_bcast(iictxt,dat,root_)
  end subroutine psb_lbcastv_ic
#endif


end module psi_bcast_mod
