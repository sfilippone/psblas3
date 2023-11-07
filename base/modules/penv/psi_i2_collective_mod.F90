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
module psi_i2_collective_mod
  use psi_penv_mod
  use psb_desc_const_mod
  
  interface psb_max
    module procedure psb_i2maxs, psb_i2maxv, psb_i2maxm
  end interface

  interface psb_min
    module procedure psb_i2mins, psb_i2minv, psb_i2minm
  end interface psb_min


  interface psb_gather
    module procedure psb_i2gather_s, psb_i2gather_v
  end interface psb_gather
  
  interface psb_gatherv
    module procedure psb_i2gatherv_v
  end interface

  interface psb_sum
    module procedure psb_i2sums, psb_i2sumv, psb_i2summ
  end interface

  interface psb_amx
    module procedure psb_i2amxs, psb_i2amxv, psb_i2amxm
  end interface

  interface psb_amn
    module procedure psb_i2amns, psb_i2amnv, psb_i2amnm
  end interface

  interface psb_bcast
    module procedure psb_i2bcasts, psb_i2bcastv, psb_i2bcastm
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_i2scan_sums, psb_i2scan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_i2exscan_sums, psb_i2exscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_i2_simple_a2av
  end interface psb_simple_a2av

  interface psb_simple_triad_a2av
    module procedure psb_i2_e_simple_triad_a2av, psb_i2_m_simple_triad_a2av
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

  subroutine psb_i2maxs(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    integer(psb_i2pk_) :: dat_

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
        call mpi_allreduce(mpi_in_place,dat,1,psb_mpi_i2pk_,mpi_max,icomm,info)
      else
        if (iam==root_) then
          call mpi_reduce(mpi_in_place,dat,1,psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,1,psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_max,icomm,request,info)
        else
          if (iam==root_) then
            call mpi_ireduce(mpi_in_place,dat,1,&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,1,&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2maxs

  subroutine psb_i2maxv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    integer(psb_i2pk_)  :: dat_(1) ! This is a dummy


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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_max,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
    
#endif    
  end subroutine psb_i2maxv

  subroutine psb_i2maxm(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    integer(psb_i2pk_)  :: dat_(1,1) ! this is a dummy


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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_max,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,root_,icomm,info)
        endif
      end if
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_max,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_max,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2maxm

  !
  ! MIN: Minimum Value
  !
  subroutine psb_i2mins(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,1,psb_mpi_i2pk_,mpi_min,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,1,psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,1,psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_min,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,1,&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,1,&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
      
#endif    
  end subroutine psb_i2mins

  subroutine psb_i2minv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_min,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_i2minv

  subroutine psb_i2minm(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_min,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,root_,icomm,info)
        end if
      end if      
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_min,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_min,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2minm



  !
  ! gather
  !
  subroutine psb_i2gather_s(ctxt,dat,resv,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat, resv(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    
#if defined(SERIAL_MPI)
    resv(0) = dat
#else
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
        call mpi_allgather(dat,1,psb_mpi_i2pk_,&
             & resv,1,psb_mpi_i2pk_,icomm,info)
      else
        call mpi_gather(dat,1,psb_mpi_i2pk_,&
             & resv,1,psb_mpi_i2pk_,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallgather(dat,1,psb_mpi_i2pk_,&
               & resv,1,psb_mpi_i2pk_,icomm,request,info)
        else
          call mpi_igather(dat,1,psb_mpi_i2pk_,&
               & resv,1,psb_mpi_i2pk_,root_,icomm,request,info)
        endif
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2gather_s

    subroutine psb_i2gather_v(ctxt,dat,resv,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:), resv(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    
#if defined(SERIAL_MPI)
    resv(0) = dat
#else
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
        call mpi_allgather(dat,size(dat),psb_mpi_i2pk_,&
             & resv,size(dat),psb_mpi_i2pk_,icomm,info)
      else
        call mpi_gather(dat,size(dat),psb_mpi_i2pk_,&
             & resv,size(dat),psb_mpi_i2pk_,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallgather(dat,size(dat),psb_mpi_i2pk_,&
               & resv,size(dat),psb_mpi_i2pk_,icomm,request,info)
        else
          call mpi_igather(dat,size(dat),psb_mpi_i2pk_,&
               & resv,size(dat),psb_mpi_i2pk_,root_,icomm,request,info)
        endif
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2gather_v

    subroutine psb_i2gatherv_v(ctxt,dat,resv,szs,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:), resv(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_), intent(in), optional    :: szs(:)
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info,i
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    integer(psb_mpk_), allocatable  :: displs(:)
    logical :: collective_start, collective_end, collective_sync
    
#if defined(SERIAL_MPI)
    resv(0) = dat
#else
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
        if (size(szs) < np) write(0,*) 'Error: bad input sizes'
        allocate(displs(np))
        displs(1) = 0
        do i=2, np
          displs(i) = displs(i-1) + szs(i-1)
        end do
        call mpi_allgatherv(dat,size(dat),psb_mpi_i2pk_,&
             & resv,szs,displs,psb_mpi_i2pk_,icomm,info)
      else
        if (iam == root_) then
          if (size(szs) < np) write(0,*) 'Error: bad input sizes'
          allocate(displs(np))
          displs(1) = 0
          do i=2, np
            displs(i) = displs(i-1) + szs(i-1)
          end do
        else
          allocate(displs(0))
        end if
        call mpi_gatherv(dat,size(dat),psb_mpi_i2pk_,&
             & resv,szs,displs,psb_mpi_i2pk_,root_,icomm,info)
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          if (size(szs) < np) write(0,*) 'Error: bad input sizes'
          allocate(displs(np))
          displs(1) = 0
          do i=2, np
            displs(i) = displs(i-1) + szs(i-1)
          end do
          call mpi_iallgatherv(dat,size(dat),psb_mpi_i2pk_,&
               & resv,szs,displs,psb_mpi_i2pk_,icomm,request,info)
        else
          if (iam == root_) then
            if (size(szs) < np) write(0,*) 'Error: bad input sizes'
            allocate(displs(np))
            displs(1) = 0
            do i=2, np
              displs(i) = displs(i-1) + szs(i-1)
            end do
          else
            allocate(displs(0))
          end if
          call mpi_igatherv(dat,size(dat),psb_mpi_i2pk_,&
               & resv,szs,displs,psb_mpi_i2pk_,root_,icomm,request,info)
        endif

      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2gatherv_v



  !
  ! SUM
  !

  subroutine psb_i2sums(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,1,&
             & psb_mpi_i2pk_,mpi_sum,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,1,&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_sum,icomm,request,info)
        else
          if(iam==root_) then                     
            call mpi_ireduce(mpi_in_place,dat,1,&
                 & psb_mpi_i2pk_,mpi_sum,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,1,&
                 & psb_mpi_i2pk_,mpi_sum,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2sums

  subroutine psb_i2sumv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
         call mpi_allreduce(mpi_in_place,dat,size(dat),&
              & psb_mpi_i2pk_,mpi_sum,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        end if
      end if
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_sum,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      endif
    end if
    
#endif    
  end subroutine psb_i2sumv

  subroutine psb_i2summ(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_sum,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,root_,icomm,info)
        end if
      end if
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_sum,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_sum,root_, icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_sum,root_, icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      endif
    end if
#endif    
  end subroutine psb_i2summ

  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_i2amxs(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,1,&
             & psb_mpi_i2pk_,mpi_i2amx_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amx_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,1,&
                 & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,1,&
                 & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_i2amxs

  subroutine psb_i2amxv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
              psb_mpi_i2pk_,mpi_i2amx_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_i2amxv

  subroutine psb_i2amxm(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_i2amx_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amx_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2amxm

  !
  ! AMN: Minimum Absolute Value
  !
  subroutine psb_i2amns(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,1,&
             & psb_mpi_i2pk_,mpi_i2amn_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,1,&
               & psb_mpi_i2pk_,mpi_i2amn_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,1,&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,1,&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_i2amns

  subroutine psb_i2amnv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_i2amn_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        end if
      endif
    else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amn_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if

#endif    
  end subroutine psb_i2amnv

  subroutine psb_i2amnm(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
        call mpi_allreduce(mpi_in_place,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_i2amn_op,icomm,info)
      else
        if(iam==root_) then 
          call mpi_reduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        else
          call mpi_reduce(dat,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,info)
        end if
      endif
          else
      if (collective_start) then
        if (root_ == -1) then 
          call mpi_iallreduce(mpi_in_place,dat,size(dat),&
               & psb_mpi_i2pk_,mpi_i2amn_op,icomm,request,info)
        else
          if(iam==root_) then 
            call mpi_ireduce(mpi_in_place,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          else
            call mpi_ireduce(dat,dat,size(dat),&
                 & psb_mpi_i2pk_,mpi_i2amn_op,root_,icomm,request,info)
          end if
        end if
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2amnm


  !
  ! BCAST Broadcast
  !  
  subroutine psb_i2bcasts(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np, info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
      call mpi_bcast(dat,1,psb_mpi_i2pk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,1,psb_mpi_i2pk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif    
  end subroutine psb_i2bcasts

  subroutine psb_i2bcastv(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
      call mpi_bcast(dat,size(dat),psb_mpi_i2pk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,size(dat),psb_mpi_i2pk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
     
#endif    
  end subroutine psb_i2bcastv

  subroutine psb_i2bcastm(ctxt,dat,root,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
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
      call mpi_bcast(dat,size(dat),psb_mpi_i2pk_,root_,icomm,info)
    else
      if (collective_start) then
        call mpi_ibcast(dat,size(dat),psb_mpi_i2pk_,root_,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
      
#endif    
  end subroutine psb_i2bcastm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_i2scan_sums(ctxt,dat,mode,request)
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
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_i2pk_) :: dat_
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
    dat_ = dat
    if (collective_sync) then 
      call mpi_scan(dat_,dat,1,&
           & psb_mpi_i2pk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iscan(dat_,dat,1,&
             & psb_mpi_i2pk_,mpi_sum,icomm,request,minfo)
      else if (collective_end) then
        call mpi_wait(request,status,minfo)
      end if
    end if
#endif    
  end subroutine psb_i2scan_sums

  subroutine psb_i2exscan_sums(ctxt,dat,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request
    integer(psb_i2pk_) :: dat_
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
    dat_ = dat
    if (collective_sync) then 
      call mpi_exscan(dat_,dat,1,&
           & psb_mpi_i2pk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iexscan(dat_,dat,1,&
             & psb_mpi_i2pk_,mpi_sum,icomm,request,minfo)
      else if (collective_end) then
        call mpi_wait(request,status,minfo)
      end if
    end if
#else
    dat = i2zero
#endif    
  end subroutine psb_i2exscan_sums

  subroutine psb_i2scan_sumv(ctxt,dat,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request

    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    integer(psb_i2pk_), allocatable :: dat_(:)
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
    dat_ = dat
    if (collective_sync) then 
      call mpi_scan(dat_,dat,size(dat),&
           & psb_mpi_i2pk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iscan(dat_,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_sum,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
#endif
  end subroutine psb_i2scan_sumv

  subroutine psb_i2exscan_sumv(ctxt,dat,mode,request)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    integer(psb_i2pk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: mode
    integer(psb_mpk_), intent(inout), optional :: request

    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo
    integer(psb_mpk_) :: icomm
    integer(psb_mpk_) :: status(mpi_status_size)
    logical :: collective_start, collective_end, collective_sync
    integer(psb_i2pk_), allocatable :: dat_(:)

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
    dat_ = dat
    if (collective_sync) then 
      call mpi_exscan(dat_,dat,size(dat),&
           & psb_mpi_i2pk_,mpi_sum,icomm,minfo)
    else
      if (collective_start) then
        call mpi_iexscan(dat_,dat,size(dat),&
             & psb_mpi_i2pk_,mpi_sum,icomm,request,info)
      else if (collective_end) then
        call mpi_wait(request,status,info)
      end if
    end if
    
#else
    dat = i2zero
#endif
  end subroutine psb_i2exscan_sumv

  subroutine psb_i2_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ctxt,info)
    use psi_i2_p2p_mod
    implicit none 
    integer(psb_i2pk_), intent(in)  :: valsnd(:)
    integer(psb_i2pk_), intent(out) :: valrcv(:)
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

  end subroutine psb_i2_simple_a2av

  subroutine psb_i2_m_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_i2pk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_i2pk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_int2_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_i2pk_,prcid(ip+1),&
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
        p2ptag =  psb_int2_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_i2pk_,prcid(ip+1),&
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

  end subroutine psb_i2_m_simple_triad_a2av

  subroutine psb_i2_e_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_i2pk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    integer(psb_i2pk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_int2_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_i2pk_,prcid(ip+1),&
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
        p2ptag = psb_int2_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_i2pk_,prcid(ip+1),&
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

  end subroutine psb_i2_e_simple_triad_a2av
  
end module psi_i2_collective_mod
