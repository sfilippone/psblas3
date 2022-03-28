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
module psi_m_collective_mod
  use psi_penv_mod
  use psb_desc_const_mod
  
  interface psb_max
    module procedure psb_mmaxs, psb_mmaxv, psb_mmaxm
  end interface

  interface psb_min
    module procedure psb_mmins, psb_mminv, psb_mminm
  end interface psb_min


  interface psb_sum
    module procedure psb_msums, psb_msumv, psb_msumm
  end interface

  interface psb_amx
    module procedure psb_mamxs, psb_mamxv, psb_mamxm
  end interface

  interface psb_amn
    module procedure psb_mamns, psb_mamnv, psb_mamnm
  end interface

  interface psb_bcast
    module procedure psb_mbcasts, psb_mbcastv, psb_mbcastm
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_mscan_sums, psb_mscan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_mexscan_sums, psb_mexscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_m_simple_a2av
  end interface psb_simple_a2av

  interface psb_simple_triad_a2av
    module procedure psb_m_e_simple_triad_a2av, psb_m_m_simple_triad_a2av
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

  subroutine psb_mmaxs(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        call mpi_allreduce(MPI_IN_PLACE,dat,1,psb_mpi_mpk_,mpi_max,icomm,info)
      else
        call mpi_reduce(MPI_IN_PLACE,dat,1,psb_mpi_mpk_,mpi_max,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_max,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_max,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_mmaxs

  subroutine psb_mmaxv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync


#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_max,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
    
#endif    
  end subroutine psb_mmaxv

  subroutine psb_mmaxm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync


#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_)&
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_max,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_max,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
   
#endif    
  end subroutine psb_mmaxm

  !
  ! MIN: Minimum Value
  !
  subroutine psb_mmins(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        call mpi_allreduce(MPI_IN_PLACE,dat,1,psb_mpi_mpk_,mpi_min,icomm,info)
      else
        call mpi_reduce(MPI_IN_PLACE,dat,1,psb_mpi_mpk_,mpi_min,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_min,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_min,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
      
#endif    
  end subroutine psb_mmins

  subroutine psb_mminv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_min,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_mminv

  subroutine psb_mminm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_)&
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_min,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,info)
        end if
      end if      
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_min,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_mminm



  !
  ! SUM
  !

  subroutine psb_msums(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync
    
#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        call mpi_allreduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_sum,icomm,info)
      else
        call mpi_reduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_sum,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_sum,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_sum,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_msums

  subroutine psb_msumv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_sum,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,root_,icomm,info)
        end if
      end if
    else
      if (collective_start) then
        if (root_ == -1) then 
          if (iinfo == psb_success_) &
               & call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,icomm,request,info)
        else
          if (iam == root_) then 
            call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
                 & psb_mpi_mpk_,mpi_sum,root_,icomm,request,info)
          else
            call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
                 & psb_mpi_mpk_,mpi_sum,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      endif
    end if
    
#endif    
  end subroutine psb_msumv

  subroutine psb_msumm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync
    
#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_sum,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,root_,icomm,info)
        end if
      end if
    else
      if (collective_start) then
        if (root_ == -1) then 
          if (iinfo == psb_success_) &
               & call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_sum,icomm,request,info)
        else
          if (iam == root_) then 
            call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
                 & psb_mpi_mpk_,mpi_sum,root_, icomm,request,info)
          else
            call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
                 & psb_mpi_mpk_,mpi_sum,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      endif
    end if
#endif    
  end subroutine psb_msumm

  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_mamxs(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        call mpi_allreduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_mamx_op,icomm,info)
      else
        call mpi_reduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_mamx_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_mamxs

  subroutine psb_mamxv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_mamx_op,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_mamxv

  subroutine psb_mamxm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_)&
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_mamx_op,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,info)
        end if
      endif
          else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamx_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_mamxm

  !
  ! AMN: Minimum Absolute Value
  !
  subroutine psb_mamns(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        call mpi_allreduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_mamn_op,icomm,info)
      else
        call mpi_reduce(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_mamn_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,1,&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_mamns

  subroutine psb_mamnv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_) &
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_mamn_op,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_mamnv

  subroutine psb_mamnm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      if (root_ == -1) then 
        if (iinfo == psb_success_)&
             & call mpi_allreduce(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_mamn_op,icomm,info)
      else
        if (iam == root_) then 
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,info)
        else
          call mpi_reduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,info)
        end if
      endif
          else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,icomm,request,info)
        else
          call mpi_ireduce(MPI_IN_PLACE,dat,size(dat),&
               & psb_mpi_mpk_,mpi_mamn_op,root_,icomm,request,info)
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_mamnm


  !
  ! BCAST Broadcast
  !  
  subroutine psb_mbcasts(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync
    

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_bcast(dat,1,psb_mpi_mpk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,1,psb_mpi_mpk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_mbcasts

  subroutine psb_mbcastv(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_bcast(dat,size(dat),psb_mpi_mpk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,size(dat),psb_mpi_mpk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
     
#endif    
  end subroutine psb_mbcastv

  subroutine psb_mbcastm(ctxt,dat,root,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_ipk_) :: iinfo
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_bcast(dat,size(dat),psb_mpi_mpk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,size(dat),psb_mpi_mpk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
      
#endif    
  end subroutine psb_mbcastm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_mscan_sums(ctxt,dat,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_mpk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_scan(MPI_IN_PLACE,dat,1,&
           & psb_mpi_mpk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iscan(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_sum,icomm,request,minfo)
      else if (collective_end) then
        call mpi_wait(request,status,minfo)
      end if
    end if
#endif    
  end subroutine psb_mscan_sums

  subroutine psb_mexscan_sums(ctxt,dat,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync


#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_exscan(MPI_IN_PLACE,dat,1,&
           & psb_mpi_mpk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_exscan(MPI_IN_PLACE,dat,1,&
             & psb_mpi_mpk_,mpi_sum,icomm,request,minfo)
      else if (collective_end) then
        call mpi_wait(request,status,minfo)
      end if
    end if
#else
    dat = mzero
#endif    
  end subroutine psb_mexscan_sums

  subroutine psb_mscan_sumv(ctxt,dat,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_scan(MPI_IN_PLACE,dat,size(dat),&
           & psb_mpi_mpk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_scan(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_sum,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif
  end subroutine psb_mscan_sumv

  subroutine psb_mexscan_sumv(ctxt,dat,mode,request)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_mpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    if (present(mode)) then
      collective_sync = .false.
      collective_start = iand(mode,psb_collective_start_) /= 0
      collective_end = iand(mode,psb_collective_end_) /= 0
      if (.not.present(request)) then
        collective_sync = .true.
        collective_start = .false.
        collective_end = .false.
      end if
    else
      collective_sync = .true.
      collective_start = .false.
      collective_end = .false.      
    end if
    if (collective_sync) then 
      call mpi_exscan(MPI_IN_PLACE,dat,size(dat),&
           & psb_mpi_mpk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iexscan(MPI_IN_PLACE,dat,size(dat),&
             & psb_mpi_mpk_,mpi_sum,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
    
#else
    dat = mzero
#endif
  end subroutine psb_mexscan_sumv

  subroutine psb_m_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ctxt,info)
    use psi_m_p2p_mod
    implicit none 
    integer(psb_mpk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: iam, np, i,j,k, ip, ipx, idx, sz

    call psb_info(ctxt,iam,np)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if

    do ip = 0, np-1
      sz = sdsz(ip+1) 
      if (sz > 0) then
        idx = bsdindx(ip+1)
        call psb_snd(ctxt,valsnd(idx+1:idx+sz),ip) 
      end if
    end do

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        idx = brvindx(ip+1)
        call psb_rcv(ctxt,valrcv(idx+1:idx+sz),ip) 
      end if
    end do

  end subroutine psb_m_simple_a2av

  subroutine psb_m_m_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_mpk_), intent(out) :: valrcv(:)
    integer(psb_mpk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_ipk_), intent(out) :: info

    !Local variables
    integer(psb_ipk_)  :: iam, np, i,j,k, ip, ipx, idx, sz, counter
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret, icomm
    integer(psb_mpk_), allocatable :: prcid(:), rvhd(:,:)

    call psb_info(ctxt,iam,np)

    icomm = psb_get_mpi_comm(ctxt)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if
    allocate(prcid(np),rvhd(np,3))
    prcid = -1

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        prcid(ip+1) = psb_get_mpi_rank(ctxt,ip)
        idx = brvindx(ip+1)
        p2ptag =  psb_int4_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
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
        if (prcid(ip+1)<0) prcid(ip+1) = psb_get_mpi_rank(ctxt,ip)
        idx = bsdindx(ip+1)
        p2ptag =  psb_int4_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
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

  end subroutine psb_m_m_simple_triad_a2av

  subroutine psb_m_e_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_mpk_), intent(out) :: valrcv(:)
    integer(psb_epk_), intent(out) :: iarcv(:), jarcv(:)
    integer(psb_mpk_), intent(in) :: bsdindx(:), brvindx(:), sdsz(:), rvsz(:)
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_ipk_), intent(out) :: info

    !Local variables
    integer(psb_ipk_)  :: iam, np, i,j,k, ip, ipx, idx, sz, counter
    integer(psb_mpk_) :: proc_to_comm, p2ptag, p2pstat(mpi_status_size), iret, icomm
    integer(psb_mpk_), allocatable :: prcid(:), rvhd(:,:)

    call psb_info(ctxt,iam,np)

    icomm = psb_get_mpi_comm(ctxt)

    if (min(size(bsdindx),size(brvindx),size(sdsz),size(rvsz))<np) then
      info = psb_err_internal_error_
      return
    end if
    allocate(prcid(np),rvhd(np,3))
    prcid = -1

    do ip = 0, np-1
      sz = rvsz(ip+1) 
      if (sz > 0) then
        prcid(ip+1) = psb_get_mpi_rank(ctxt,ip)
        idx = brvindx(ip+1)
        p2ptag =  psb_int4_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
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
        if (prcid(ip+1)<0) prcid(ip+1) = psb_get_mpi_rank(ctxt,ip)
        idx = bsdindx(ip+1)
        p2ptag = psb_int4_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_mpk_,prcid(ip+1),&
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

  end subroutine psb_m_e_simple_triad_a2av

  
end module psi_m_collective_mod
