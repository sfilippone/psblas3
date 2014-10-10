!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psi_reduce_mod
  use psi_penv_mod
  interface psb_max
    module procedure psb_imaxs, psb_imaxv, psb_imaxm,&
         & psb_smaxs, psb_smaxv, psb_smaxm,&
         & psb_dmaxs, psb_dmaxv, psb_dmaxm
  end interface
#if defined(LONG_INTEGERS)
  interface psb_max 
    module procedure psb_i4maxs, psb_i4maxv, psb_i4maxm
  end interface
#endif
#if !defined(LONG_INTEGERS)
  interface psb_max 
    module procedure psb_i8maxs, psb_i8maxv, psb_i8maxm
  end interface
#endif

  interface psb_min
    module procedure psb_imins, psb_iminv, psb_iminm,&
         & psb_smins, psb_sminv, psb_sminm,&
         & psb_dmins, psb_dminv, psb_dminm
  end interface
#if !defined(LONG_INTEGERS)
  interface psb_min
    module procedure psb_i8mins, psb_i8minv, psb_i8minm
  end interface
#endif
#if defined(LONG_INTEGERS)
  interface psb_min
    module procedure psb_i4mins, psb_i4minv, psb_i4minm
  end interface
#endif


  interface psb_amx
    module procedure psb_iamxs, psb_iamxv, psb_iamxm,&
         & psb_samxs, psb_samxv, psb_samxm,&
         & psb_camxs, psb_camxv, psb_camxm,&
         & psb_damxs, psb_damxv, psb_damxm,&
         & psb_zamxs, psb_zamxv, psb_zamxm
  end interface
#if !defined(LONG_INTEGERS)
  interface psb_amx 
    module procedure psb_i8amxs, psb_i8amxv, psb_i8amxm
  end interface
#endif
#if defined(LONG_INTEGERS)
  interface psb_amx 
    module procedure psb_i4amxs, psb_i4amxv, psb_i4amxm
  end interface
#endif

  interface psb_amn
    module procedure psb_iamns, psb_iamnv, psb_iamnm,&
         & psb_samns, psb_samnv, psb_samnm,&
         & psb_camns, psb_camnv, psb_camnm,&
         & psb_damns, psb_damnv, psb_damnm,&
         & psb_zamns, psb_zamnv, psb_zamnm
  end interface
#if defined(LONG_INTEGERS)
  interface psb_amn 
    module procedure psb_i4amns, psb_i4amnv, psb_i4amnm
  end interface
#endif
#if !defined(LONG_INTEGERS)
  interface psb_amn 
    module procedure psb_i8amns, psb_i8amnv, psb_i8amnm
  end interface
#endif


  interface psb_sum
    module procedure psb_isums, psb_isumv, psb_isumm,&
         & psb_ssums, psb_ssumv, psb_ssumm,&
         & psb_csums, psb_csumv, psb_csumm,&
         & psb_dsums, psb_dsumv, psb_dsumm,&
         & psb_zsums, psb_zsumv, psb_zsumm
  end interface
#if defined(SHORT_INTEGERS)
  interface psb_sum
    module procedure psb_i2sums, psb_i2sumv, psb_i2summ
  end interface psb_sum
#endif
#if defined(LONG_INTEGERS)
  interface psb_sum
    module procedure psb_i4sums, psb_i4sumv, psb_i4summ
  end interface
#endif
#if !defined(LONG_INTEGERS)
  interface psb_sum
    module procedure psb_i8sums, psb_i8sumv, psb_i8summ
  end interface
#endif


  interface psb_nrm2
    module procedure psb_s_nrm2s, psb_s_nrm2v,&
         & psb_d_nrm2s, psb_d_nrm2v
  end interface



#if defined(LONG_INTEGERS)
  interface psb_max
    module procedure psb_imaxs_ic, psb_imaxv_ic, psb_imaxm_ic,&
         & psb_smaxs_ic, psb_smaxv_ic, psb_smaxm_ic,&
         & psb_dmaxs_ic, psb_dmaxv_ic, psb_dmaxm_ic
  end interface

  interface psb_min
    module procedure psb_imins_ic, psb_iminv_ic, psb_iminm_ic,&
         & psb_smins_ic, psb_sminv_ic, psb_sminm_ic,&
         & psb_dmins_ic, psb_dminv_ic, psb_dminm_ic
  end interface


  interface psb_amx
    module procedure psb_iamxs_ic, psb_iamxv_ic, psb_iamxm_ic,&
         & psb_samxs_ic, psb_samxv_ic, psb_samxm_ic,&
         & psb_camxs_ic, psb_camxv_ic, psb_camxm_ic,&
         & psb_damxs_ic, psb_damxv_ic, psb_damxm_ic,&
         & psb_zamxs_ic, psb_zamxv_ic, psb_zamxm_ic
  end interface

  interface psb_amn
    module procedure psb_iamns_ic, psb_iamnv_ic, psb_iamnm_ic,&
         & psb_samns_ic, psb_samnv_ic, psb_samnm_ic,&
         & psb_camns_ic, psb_camnv_ic, psb_camnm_ic,&
         & psb_damns_ic, psb_damnv_ic, psb_damnm_ic,&
         & psb_zamns_ic, psb_zamnv_ic, psb_zamnm_ic
  end interface


  interface psb_sum
    module procedure psb_isums_ic, psb_isumv_ic, psb_isumm_ic,&
         & psb_ssums_ic, psb_ssumv_ic, psb_ssumm_ic,&
         & psb_csums_ic, psb_csumv_ic, psb_csumm_ic,&
         & psb_dsums_ic, psb_dsumv_ic, psb_dsumm_ic,&
         & psb_zsums_ic, psb_zsumv_ic, psb_zsumm_ic
  end interface
#if defined(SHORT_INTEGERS)
  interface psb_sum
    module procedure psb_i2sums_ic, psb_i2sumv_ic, psb_i2summ_ic
  end interface psb_sum
#endif

  interface psb_nrm2
    module procedure psb_s_nrm2s_ic, psb_s_nrm2v_ic,&
         & psb_d_nrm2s_ic, psb_d_nrm2v_ic
  end interface
#endif


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

  subroutine psb_imaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_ipk_) :: dat_
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_imaxs

  subroutine psb_imaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_imaxv

  subroutine psb_imaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_imaxm

#if defined(LONG_INTEGERS)
  subroutine psb_i4maxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) ::  dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_i4maxs

  subroutine psb_i4maxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4maxv

  subroutine psb_i4maxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4maxm
  
#endif


#if !defined(LONG_INTEGERS)
  subroutine psb_i8maxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_) ::  dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_lng_integer,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_lng_integer,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_i8maxs

  subroutine psb_i8maxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8maxv

  subroutine psb_i8maxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8maxm

#endif


  subroutine psb_smaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_smaxs

  subroutine psb_smaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_smaxv

  subroutine psb_smaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_smaxm

  subroutine psb_dmaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_max,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_max,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_dmaxs

  subroutine psb_dmaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dmaxv

  subroutine psb_dmaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_max,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_max,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_max,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dmaxm


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  MIN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_imins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_ipk_) :: dat_
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_imins

  subroutine psb_iminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iminv

  subroutine psb_iminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iminm

#if defined(LONG_INTEGERS)
  subroutine psb_i4mins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) ::  dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_i4mins

  subroutine psb_i4minv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4minv

  subroutine psb_i4minm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4minm

#endif


#if !defined(LONG_INTEGERS)
  subroutine psb_i8mins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_) ::  dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_lng_integer,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_lng_integer,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_i8mins

  subroutine psb_i8minv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8minv

  subroutine psb_i8minm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8minm

#endif



  subroutine psb_smins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_smins

  subroutine psb_sminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_sminv

  subroutine psb_sminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_sminm

  subroutine psb_dmins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_min,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_min,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_dmins

  subroutine psb_dminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dminv

  subroutine psb_dminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_min,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_min,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_min,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dminm


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  AMX: maximum absolute value
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_iamxs(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_iamx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_iamx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_iamxs

  subroutine psb_iamxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_iamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_iamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_iamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iamxv

  subroutine psb_iamxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_iamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_iamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_iamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iamxm


#if defined(LONG_INTEGERS)
  subroutine psb_i4amxs(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer,mpi_i4amx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer,mpi_i4amx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i4amxs

  subroutine psb_i4amxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_i4amx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_i4amx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_i4amx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4amxv

  subroutine psb_i4amxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_i4amx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_i4amx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_i4amx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4amxm

#endif

#if !defined(LONG_INTEGERS)
  subroutine psb_i8amxs(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_lng_integer,mpi_i8amx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_lng_integer,mpi_i8amx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i8amxs

  subroutine psb_i8amxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_i8amx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_i8amx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_i8amx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8amxv

  subroutine psb_i8amxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_i8amx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_i8amx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_i8amx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8amxm

#endif



  subroutine psb_samxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_samxs

  subroutine psb_samxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_samxv

  subroutine psb_samxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_samxm

  subroutine psb_damxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_damx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_damx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_damxs

  subroutine psb_damxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_damx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_damx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_damxv

  subroutine psb_damxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_damx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_damxm


  subroutine psb_camxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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

  subroutine psb_zamxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamx_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamx_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zamxs

  subroutine psb_zamxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,&
           & mpi_zamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zamxv

  subroutine psb_zamxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamx_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zamxm



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  AMN: minimum absolute value
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_iamns(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_iamn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_iamn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_iamns

  subroutine psb_iamnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_iamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_iamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_iamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iamnv

  subroutine psb_iamnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_iamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_iamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_iamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_iamnm


#if defined(LONG_INTEGERS)
  subroutine psb_i4amns(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer,mpi_i4amn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer,mpi_i4amn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i4amns

  subroutine psb_i4amnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_i4amn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_i4amn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_i4amn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4amnv

  subroutine psb_i4amnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
    integer(psb_mpik_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_i4amn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_i4amn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_i4amn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4amnm

#endif

#if !defined(LONG_INTEGERS)
  subroutine psb_i8amns(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_lng_integer,mpi_i8amn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_lng_integer,mpi_i8amn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i8amns

  subroutine psb_i8amnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_i8amn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_i8amn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_i8amn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8amnv

  subroutine psb_i8amnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_i8amn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_i8amn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_i8amn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8amnm

#endif



  subroutine psb_samns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_samn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_samns

  subroutine psb_samnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_samnv

  subroutine psb_samnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_samn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_samnm

  subroutine psb_damns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_damn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_damn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_damns

  subroutine psb_damnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_damn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_damn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_damnv

  subroutine psb_damnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_damn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_damn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_damnm


  subroutine psb_camns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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

  subroutine psb_zamns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamn_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_zamn_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zamns

  subroutine psb_zamnv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,&
           & mpi_zamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zamnv

  subroutine psb_zamnm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_zamn_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zamnm


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SUM
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_isums(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_ipk_integer,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_isums

  subroutine psb_isumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)            :: ictxt
    integer(psb_ipk_), intent(inout)          :: dat(:)
    integer(psb_mpik_), intent(in), optional  :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_)&
           & call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_isumv

  subroutine psb_isumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_ipk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_ipk_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_ipk_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_ipk_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_isumm


#if defined(SHORT_INTEGERS)
  subroutine psb_i2sums(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(2) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer2,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer2,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i2sums

  subroutine psb_i2sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(2), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer2,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer2,mpi_sum,root_,ictxt,info)
      else
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer2,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i2sumv

  subroutine psb_i2summ(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(2), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer2,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer2,mpi_sum,root_,ictxt,info)
      else
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer2,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i2summ

#endif


#if defined(LONG_INTEGERS)
  subroutine psb_i4sums(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_def_integer,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_def_integer,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i4sums

  subroutine psb_i4sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_=dat
      if (info == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,info)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4sumv

  subroutine psb_i4summ(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_mpik_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_mpik_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_=dat
      if (info == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_def_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_def_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,info)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_def_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i4summ

#endif

#if !defined(LONG_INTEGERS)
  subroutine psb_i8sums(ictxt,dat,root)

#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo

#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_lng_integer,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_lng_integer,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif

#endif    
  end subroutine psb_i8sums

  subroutine psb_i8sumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8sumv

  subroutine psb_i8summ(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    integer(psb_long_int_k_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
      dat_=dat
      if (iinfo == psb_success_) call mpi_allreduce(dat_,dat,size(dat),&
           & psb_mpi_lng_integer,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_lng_integer,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_lng_integer,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_i8summ

#endif



  subroutine psb_ssums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_ssums

  subroutine psb_ssumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_ssumv

  subroutine psb_ssumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_ssumm

  subroutine psb_dsums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_dsums

  subroutine psb_dsumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dsumv

  subroutine psb_dsumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_dsumm


  subroutine psb_csums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_spk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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

  subroutine psb_zsums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_c_dpk_,mpi_sum,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_zsums

  subroutine psb_zsumv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,&
           & mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zsumv

  subroutine psb_zsumm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    complex(psb_dpk_), allocatable :: dat_(:,:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,ictxt,info)
      else
        call psb_realloc(1,1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_c_dpk_,mpi_sum,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_zsumm

  ! !!!!!!!!!!!!
  !
  ! Norm 2
  !
  ! !!!!!!!!!!!!
  subroutine psb_s_nrm2s(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)            :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional  :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_spk_,mpi_snrm2_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_spk_,mpi_snrm2_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_s_nrm2s

  subroutine psb_d_nrm2s(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)            :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_mpik_), intent(in), optional  :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_) :: dat_
    integer(psb_mpik_) :: iam, np, info
    integer(psb_ipk_) :: iinfo


#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_dnrm2_op,ictxt,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,psb_mpi_r_dpk_,mpi_dnrm2_op,root_,ictxt,info)
      if (iam == root_) dat = dat_
    endif
#endif    
  end subroutine psb_d_nrm2s

  subroutine psb_s_nrm2v(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional    :: root
    integer(psb_mpik_) :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_spk_,&
           & mpi_snrm2_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_spk_,&
             & mpi_snrm2_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_spk_,&
             & mpi_snrm2_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_s_nrm2v

  subroutine psb_d_nrm2v(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_), intent(in)            :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_mpik_), intent(in), optional  :: root
    integer(psb_mpik_) :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer(psb_mpik_) :: iam, np,  info
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
           & call mpi_allreduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
           & mpi_dnrm2_op,ictxt,info)
    else
      if (iam == root_) then 
        call psb_realloc(size(dat),dat_,iinfo)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),psb_mpi_r_dpk_,&
             & mpi_dnrm2_op,root_,ictxt,info)
      else
        call psb_realloc(1,dat_,iinfo)
        call mpi_reduce(dat,dat_,size(dat),psb_mpi_r_dpk_,&
             & mpi_dnrm2_op,root_,ictxt,info)
      end if
    endif
#endif    
  end subroutine psb_d_nrm2v

#if defined(LONG_INTEGERS)

  subroutine psb_imaxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if

  end subroutine psb_imaxs_ic

  subroutine psb_imaxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_imaxv_ic

  subroutine psb_imaxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_imaxm_ic


  subroutine psb_smaxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_smaxs_ic

  subroutine psb_smaxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_smaxv_ic

  subroutine psb_smaxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_smaxm_ic

  subroutine psb_dmaxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_dmaxs_ic

  subroutine psb_dmaxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_dmaxv_ic

  subroutine psb_dmaxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_max(ictxt_,dat,root_)
    else
      call psb_max(ictxt_,dat)
    end if
  end subroutine psb_dmaxm_ic


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  MIN
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_imins_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_imins_ic

  subroutine psb_iminv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_iminv_ic

  subroutine psb_iminm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_iminm_ic


  subroutine psb_smins_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_smins_ic

  subroutine psb_sminv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_sminv_ic

  subroutine psb_sminm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_sminm_ic

  subroutine psb_dmins_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_dmins_ic

  subroutine psb_dminv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_dminv_ic

  subroutine psb_dminm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_min(ictxt_,dat,root_)
    else
      call psb_min(ictxt_,dat)
    end if
  end subroutine psb_dminm_ic


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  AMX: maximum absolute value
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_iamxs_ic(ictxt,dat,root)

    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_iamxs_ic

  subroutine psb_iamxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_iamxv_ic

  subroutine psb_iamxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_iamxm_ic



  subroutine psb_samxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_samxs_ic

  subroutine psb_samxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_samxv_ic

  subroutine psb_samxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_samxm_ic

  subroutine psb_damxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_damxs_ic

  subroutine psb_damxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_damxv_ic

  subroutine psb_damxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_damxm_ic


  subroutine psb_camxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxs_ic

  subroutine psb_camxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxv_ic

  subroutine psb_camxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_camxm_ic

  subroutine psb_zamxs_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_zamxs_ic

  subroutine psb_zamxv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_zamxv_ic

  subroutine psb_zamxm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amx(ictxt_,dat,root_)
    else
      call psb_amx(ictxt_,dat)
    end if
  end subroutine psb_zamxm_ic



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  AMN: minimum absolute value
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_iamns_ic(ictxt,dat,root)

    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_iamns_ic

  subroutine psb_iamnv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_iamnv_ic

  subroutine psb_iamnm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_iamnm_ic




  subroutine psb_samns_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_samns_ic

  subroutine psb_samnv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_samnv_ic

  subroutine psb_samnm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_samnm_ic

  subroutine psb_damns_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_damns_ic

  subroutine psb_damnv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_damnv_ic

  subroutine psb_damnm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_damnm_ic


  subroutine psb_camns_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camns_ic

  subroutine psb_camnv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camnv_ic

  subroutine psb_camnm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_camnm_ic

  subroutine psb_zamns_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_zamns_ic

  subroutine psb_zamnv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_zamnv_ic

  subroutine psb_zamnm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_amn(ictxt_,dat,root_)
    else
      call psb_amn(ictxt_,dat)
    end if
  end subroutine psb_zamnm_ic


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !  SUM
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine psb_isums_ic(ictxt,dat,root)

    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_isums_ic

  subroutine psb_isumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_isumv_ic

  subroutine psb_isumm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(psb_ipk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_isumm_ic

#if defined(SHORT_INTEGERS) 
  subroutine psb_i2sums_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_i2sums_ic

  subroutine psb_i2sumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_i2sumv_ic

  subroutine psb_i2summ_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    integer(2), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_i2summ_ic

#endif




  subroutine psb_ssums_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_ssums_ic

  subroutine psb_ssumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_ssumv_ic

  subroutine psb_ssumm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_ssumm_ic

  subroutine psb_dsums_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_dsums_ic

  subroutine psb_dsumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_dsumv_ic

  subroutine psb_dsumm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_dsumm_ic


  subroutine psb_csums_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csums_ic

  subroutine psb_csumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csumv_ic

  subroutine psb_csumm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_csumm_ic

  subroutine psb_zsums_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_zsums_ic

  subroutine psb_zsumv_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_zsumv_ic

  subroutine psb_zsumm_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_sum(ictxt_,dat,root_)
    else
      call psb_sum(ictxt_,dat)
    end if
  end subroutine psb_zsumm_ic

  ! !!!!!!!!!!!!
  !
  ! Norm 2
  !
  ! !!!!!!!!!!!!
  subroutine psb_s_nrm2s_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)            :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional  :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_nrm2(ictxt_,dat,root_)
    else
      call psb_nrm2(ictxt_,dat)
    end if
  end subroutine psb_s_nrm2s_ic

  subroutine psb_d_nrm2s_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)            :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer(psb_ipk_), intent(in), optional  :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_nrm2(ictxt_,dat,root_)
    else
      call psb_nrm2(ictxt_,dat)
    end if
  end subroutine psb_d_nrm2s_ic

  subroutine psb_s_nrm2v_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional    :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_nrm2(ictxt_,dat,root_)
    else
      call psb_nrm2(ictxt_,dat)
    end if
  end subroutine psb_s_nrm2v_ic

  subroutine psb_d_nrm2v_ic(ictxt,dat,root)
    implicit none 
    integer(psb_ipk_), intent(in)            :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer(psb_ipk_), intent(in), optional  :: root
    integer(psb_mpik_) :: ictxt_, root_

    ictxt_ = ictxt
    if (present(root)) then 
      root_ = root
      call psb_nrm2(ictxt_,dat,root_)
    else
      call psb_nrm2(ictxt_,dat)
    end if
  end subroutine psb_d_nrm2v_ic

#endif
end module psi_reduce_mod
