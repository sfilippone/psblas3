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
module psi_z_collective_mod
  use psi_penv_mod


  interface psb_sum
    module procedure psb_zsums, psb_zsumv, psb_zsumm
  end interface

  interface psb_amx
    module procedure psb_zamxs, psb_zamxv, psb_zamxm
  end interface

  interface psb_amn
    module procedure psb_zamns, psb_zamnv, psb_zamnm
  end interface

  interface psb_bcast
    module procedure psb_zbcasts, psb_zbcastv, psb_zbcastm
  end interface psb_bcast

  interface psb_scan_sum
    module procedure psb_zscan_sums, psb_zscan_sumv
  end interface psb_scan_sum

  interface psb_exscan_sum
    module procedure psb_zexscan_sums, psb_zexscan_sumv
  end interface psb_exscan_sum

  interface psb_simple_a2av
    module procedure psb_z_simple_a2av
  end interface psb_simple_a2av

  interface psb_simple_triad_a2av
    module procedure psb_z_e_simple_triad_a2av, psb_z_m_simple_triad_a2av
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

  subroutine psb_zsums(ctxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zsums

  subroutine psb_zsumv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zsumv

  subroutine psb_zsumm(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zsumm

  !
  ! AMX: Maximum Absolute Value
  !
  
  subroutine psb_zamxs(ctxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamx_op,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamx_op,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zamxs

  subroutine psb_zamxv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zamxv

  subroutine psb_zamxm(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zamxm

  !
  ! AMN: Minimum Absolute Value
  !
  
  subroutine psb_zamns(ctxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamn_op,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamn_op,root_,icomm,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zamns

  subroutine psb_zamnv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_) &
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,icomm,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zamnv

  subroutine psb_zamnm(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    icomm = psb_get_mpi_comm(ctxt)
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_ = dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,icomm,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,icomm,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_zamnm

  !
  ! BCAST Broadcast
  !
  
  subroutine psb_zbcasts(ctxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np, info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_bcast(dat,1,psb_mpi_c_dpk_,root_,icomm,info)

#endif    
  end subroutine psb_zbcasts

  subroutine psb_zbcastv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_bcast(dat,size(dat),psb_mpi_c_dpk_,root_,icomm,info)
#endif    
  end subroutine psb_zbcastv

  subroutine psb_zbcastm(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_

    integer(psb_mpk_) :: iam, np,  info, icomm
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ctxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    endif
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_bcast(dat,size(dat),psb_mpi_c_dpk_,root_,icomm,info)
#endif    
  end subroutine psb_zbcastm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SCAN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_zscan_sums(ctxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    complex(psb_dpk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_scan(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,icomm,minfo)
    dat = dat_
#endif    
  end subroutine psb_zscan_sums


  subroutine psb_zexscan_sums(ctxt,dat)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat
    complex(psb_dpk_) :: dat_
    integer(psb_ipk_) :: iam, np, info
    integer(psb_mpk_) :: icomm, minfo


#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    call mpi_exscan(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,icomm,minfo)
    dat = dat_
#else
    dat = zzero
#endif    
  end subroutine psb_zexscan_sums

  subroutine psb_zscan_sumv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_scan(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,icomm,minfo)
#endif
  end subroutine psb_zscan_sumv

  subroutine psb_zexscan_sumv(ctxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in)              :: ctxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpk_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_ipk_) :: iam, np,  info
    integer(psb_mpk_) :: minfo, icomm

#if !defined(SERIAL_MPI)
    call psb_info(ctxt,iam,np)
    icomm = psb_get_mpi_comm(ctxt)
    call psb_realloc(size(dat),dat_,info)
    dat_ = dat
    if (info == psb_success_) &
         & call mpi_exscan(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,icomm,minfo)
#else
    dat = zzero
#endif
  end subroutine psb_zexscan_sumv

  subroutine psb_z_simple_a2av(valsnd,sdsz,bsdindx,&
       & valrcv,rvsz,brvindx,ctxt,info)
    use psi_z_p2p_mod
    implicit none 
    complex(psb_dpk_), intent(in)  :: valsnd(:)
    complex(psb_dpk_), intent(out) :: valrcv(:)
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

  end subroutine psb_z_simple_a2av

  subroutine psb_z_m_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    complex(psb_dpk_), intent(in)  :: valsnd(:)
    integer(psb_mpk_), intent(in)  :: iasnd(:), jasnd(:)
    complex(psb_dpk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_dcomplex_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_c_dpk_,prcid(ip+1),&
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
        p2ptag =  psb_dcomplex_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_c_dpk_,prcid(ip+1),&
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

  end subroutine psb_z_m_simple_triad_a2av

  subroutine psb_z_e_simple_triad_a2av(valsnd,iasnd,jasnd,sdsz,bsdindx,&
       & valrcv,iarcv,jarcv,rvsz,brvindx,ctxt,info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    complex(psb_dpk_), intent(in)  :: valsnd(:)
    integer(psb_epk_), intent(in)  :: iasnd(:), jasnd(:)
    complex(psb_dpk_), intent(out) :: valrcv(:)
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
        p2ptag =  psb_dcomplex_tag
        call mpi_irecv(valrcv(idx+1:idx+sz),sz,&
             & psb_mpi_c_dpk_,prcid(ip+1),&
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
        p2ptag = psb_dcomplex_tag
        call mpi_send(valsnd(idx+1:idx+sz),sz,&
             & psb_mpi_c_dpk_,prcid(ip+1),&
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

  end subroutine psb_z_e_simple_triad_a2av

  
end module psi_z_collective_mod
