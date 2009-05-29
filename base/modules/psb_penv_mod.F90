!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
#if defined(SERIAL_MPI)
! Provide a fake mpi module just to keep the compiler(s) happy.
module mpi
  integer, parameter :: mpi_success=0
  integer, parameter :: mpi_request_null=0
  integer, parameter :: mpi_status_size=1
  integer, parameter :: mpi_integer=1, mpi_double_precision=3
  integer, parameter :: mpi_double_complex=5 
  real(psb_dpk_), external :: mpi_wtime
end module mpi
#endif    

module psb_penv_mod
  
  use psb_const_mod
  use psb_blacs_mod 

  interface psb_init
    module procedure  psb_init
  end interface

  interface psb_exit
    module procedure  psb_exit
  end interface

  interface psb_abort
    module procedure  psb_abort
  end interface

  interface psb_info
    module procedure psb_info
  end interface
  
  interface psb_barrier
    module procedure  psb_barrier
  end interface

  interface psb_wtime
    module procedure  psb_wtime
  end interface

  interface psb_bcast
    module procedure psb_ibcasts, psb_ibcastv, psb_ibcastm,&
         & psb_dbcasts, psb_dbcastv, psb_dbcastm,&
         & psb_zbcasts, psb_zbcastv, psb_zbcastm,&
         & psb_sbcasts, psb_sbcastv, psb_sbcastm,&
         & psb_cbcasts, psb_cbcastv, psb_cbcastm,&
         & psb_hbcasts, psb_hbcastv, psb_lbcasts, psb_lbcastv
  end interface


  interface psb_snd
    module procedure psb_isnds, psb_isndv, psb_isndm,&
         & psb_ssnds, psb_ssndv, psb_ssndm,&
         & psb_csnds, psb_csndv, psb_csndm,&
         & psb_dsnds, psb_dsndv, psb_dsndm,&
         & psb_zsnds, psb_zsndv, psb_zsndm,&
         & psb_hsnds, psb_lsnds
  end interface

  interface psb_rcv
    module procedure psb_ircvs, psb_ircvv, psb_ircvm,&
         & psb_srcvs, psb_srcvv, psb_srcvm,&
         & psb_crcvs, psb_crcvv, psb_crcvm,&
         & psb_drcvs, psb_drcvv, psb_drcvm,&
         & psb_zrcvs, psb_zrcvv, psb_zrcvm,&
         & psb_hrcvs, psb_lrcvs
  end interface

  interface psb_max
    module procedure psb_imaxs, psb_imaxv, psb_imaxm,&
         & psb_i8maxs, &
         & psb_smaxs, psb_smaxv, psb_smaxm,&
         & psb_dmaxs, psb_dmaxv, psb_dmaxm
  end interface


  interface psb_min
    module procedure psb_imins, psb_iminv, psb_iminm,&
         & psb_i8mins, &
         & psb_smins, psb_sminv, psb_sminm,&
         & psb_dmins, psb_dminv, psb_dminm
  end interface


  interface psb_amx
    module procedure psb_iamxs, psb_iamxv, psb_iamxm,&
         & psb_i8amxs, &
         & psb_samxs, psb_samxv, psb_samxm,&
         & psb_camxs, psb_camxv, psb_camxm,&
         & psb_damxs, psb_damxv, psb_damxm,&
         & psb_zamxs, psb_zamxv, psb_zamxm
  end interface

  interface psb_amn
    module procedure psb_iamns, psb_iamnv, psb_iamnm,&
         & psb_i8amns, &
         & psb_samns, psb_samnv, psb_samnm,&
         & psb_camns, psb_camnv, psb_camnm,&
         & psb_damns, psb_damnv, psb_damnm,&
         & psb_zamns, psb_zamnv, psb_zamnm
  end interface

  interface psb_sum
    module procedure psb_isums, psb_isumv, psb_isumm,&
         & psb_i8sums, psb_i8sumv,&
         & psb_ssums, psb_ssumv, psb_ssumm,&
         & psb_csums, psb_csumv, psb_csumm,&
         & psb_dsums, psb_dsumv, psb_dsumm,&
         & psb_zsums, psb_zsumv, psb_zsumm
  end interface


#if defined(SERIAL_MPI)
  integer, private, save :: nctxt=0
#endif


#if defined(HAVE_KSENDID)
  interface 
    integer function krecvid(contxt,proc_to_comm,myrow)
      integer contxt,proc_to_comm,myrow
    end function krecvid
  end interface
  interface 
    integer function ksendid(contxt,proc_to_comm,myrow)
      integer contxt,proc_to_comm,myrow
    end function ksendid
  end interface
#endif  
  private psi_get_sizes
contains 

  subroutine psi_get_sizes()
    use psb_const_mod
    real(psb_dpk_) :: dv(2) 
    real(psb_spk_) :: sv(2) 
    integer        :: iv(2)
    integer(psb_long_int_k_) :: ilv(2)
    
    call psi_c_diffadd(sv(1),sv(2),psb_sizeof_sp)
    call psi_c_diffadd(dv(1),dv(2),psb_sizeof_dp)
    call psi_c_diffadd(iv(1),iv(2),psb_sizeof_int)
    call psi_c_diffadd(ilv(1),ilv(2),psb_sizeof_long_int)
    
  end subroutine psi_get_sizes
  
  subroutine psb_init(ictxt,np)
    use psb_const_mod
    use psb_error_mod
    integer, intent(out) :: ictxt
    integer, intent(in), optional :: np
    
    integer :: np_, npavail, iam, info
    character(len=20), parameter :: name='psb_init'
#if defined(SERIAL_MPI) 
    ictxt = nctxt
    nctxt = nctxt + 1
    np_   = 1
#else    
    call blacs_pinfo(iam, npavail)
    call blacs_get(izero, izero, ictxt)
    
    if (present(np)) then 
      np_ = max(1,min(np,npavail))
    else
      np_ = npavail
    endif
    
    call blacs_gridinit(ictxt, 'R', np_, ione)
#endif    
    if (present(np)) then 
      if (np_ < np) then 
        info = 2011
        call psb_errpush(info,name)
        call psb_error(ictxt)
      endif
    endif
    call psi_get_sizes()

  end subroutine psb_init

  subroutine psb_exit(ictxt,close)
    integer, intent(in) :: ictxt
    logical, intent(in), optional :: close
    logical  :: close_
    integer  :: nprow, npcol, myprow, mypcol

#if !defined(SERIAL_MPI)
    if (present(close)) then 
      close_ = close
    else
      close_ = .true.
    end if
    call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
    if ((myprow >=0).and.(mypcol>=0)) then
      call blacs_gridexit(ictxt)
    end if
    if (close_) call blacs_exit(0)
#endif
  end subroutine psb_exit


  subroutine psb_barrier(ictxt)
    integer, intent(in) :: ictxt

#if !defined(SERIAL_MPI)
    call blacs_barrier(ictxt,'All')
#endif    
    
  end subroutine psb_barrier

  function psb_wtime()
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    real(psb_dpk_) :: psb_wtime

    psb_wtime = mpi_wtime()
  end function psb_wtime

  subroutine psb_abort(ictxt)
    integer, intent(in) :: ictxt
    
#if defined(SERIAL_MPI) 
    stop
#else    
    call blacs_abort(ictxt,-1)
#endif    
    
  end subroutine psb_abort


  subroutine psb_info(ictxt,iam,np)

    integer, intent(in)  :: ictxt
    integer, intent(out) :: iam, np
    integer              :: nprow, npcol, myprow, mypcol
    
#if defined(SERIAL_MPI) 
    iam = 0
    np  = 1
#else    
    call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
    
    iam = myprow
    np  = nprow
#endif    
    
  end subroutine psb_info

    
  subroutine psb_ibcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    integer, intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_ibcasts

  subroutine psb_ibcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    integer, intent(inout) :: dat(:)
    integer, intent(in), optional  :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_ibcastv
    
  subroutine psb_ibcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    integer, intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_ibcastm


  subroutine psb_sbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_sbcasts


  subroutine psb_sbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_sbcastv
    
  subroutine psb_sbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_sbcastm
    


  subroutine psb_dbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_dbcasts


  subroutine psb_dbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_dbcastv
    
  subroutine psb_dbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_dbcastm
    

  subroutine psb_cbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_cbcasts

  subroutine psb_cbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_cbcastv
    
  subroutine psb_cbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_cbcastm

  subroutine psb_zbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_zbcasts

  subroutine psb_zbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_zbcastv
    
  subroutine psb_zbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
#endif    
  end subroutine psb_zbcastm


  subroutine psb_hbcasts(ictxt,dat,root,length)
#ifdef MPI_H
    include 'mpif.h'
#endif
#ifdef MPI_MOD
    use mpi
#endif
    integer, intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer, intent(in), optional   :: root,length

    integer  :: iam, np, root_,icomm,length_,info

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
    call psb_get_mpicomm(ictxt,icomm)
    
    call mpi_bcast(dat,length_,MPI_CHARACTER,root_,icomm,info)
#endif    

  end subroutine psb_hbcasts

  subroutine psb_hbcastv(ictxt,dat,root)
#ifdef MPI_H
    include 'mpif.h'
#endif
#ifdef MPI_MOD
    use mpi
#endif
    integer, intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat(:)
    integer, intent(in), optional   :: root

    integer  :: iam, np, root_,icomm,length_,info, size_

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ =  psb_root_
    endif
    length_ = len(dat)
    size_   = size(dat) 

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    
    call mpi_bcast(dat,length_*size_,MPI_CHARACTER,root_,icomm,info)
#endif    

  end subroutine psb_hbcastv

  subroutine psb_lbcasts(ictxt,dat,root)
#ifdef MPI_H
    include 'mpif.h'
#endif
#ifdef MPI_MOD
    use mpi
#endif
    integer, intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer, intent(in), optional   :: root

    integer  :: iam, np, root_,icomm,info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    call mpi_bcast(dat,1,MPI_LOGICAL,root_,icomm,info)
#endif    

  end subroutine psb_lbcasts


  subroutine psb_lbcastv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)             :: ictxt
    logical, intent(inout)          :: dat(:)
    integer, intent(in), optional   :: root

    integer  :: iam, np, root_,icomm,info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    call mpi_bcast(dat,size(dat),MPI_LOGICAL,root_,icomm,info)
#endif    

  end subroutine psb_lbcastv



  subroutine psb_imaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_, dat_
    integer :: iam, np, icomm,info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer,mpi_max,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer,mpi_max,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_imaxs

  subroutine psb_imaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_=dat
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_imaxv

  subroutine psb_imaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_=dat
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_imaxm
  
  subroutine psb_smaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_) :: dat_
    integer :: iam, np, icomm,info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_real,mpi_max,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_real,mpi_max,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_smaxs

  subroutine psb_smaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_real,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_real,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_real,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_smaxv

  subroutine psb_smaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_ = dat
      if (info ==0)&
           & call mpi_allreduce(dat_,dat,size(dat),mpi_real,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_real,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_real,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_smaxm
  
  subroutine psb_dmaxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_) :: dat_
    integer :: iam, np, icomm,info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_double_precision,mpi_max,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_double_precision,mpi_max,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_dmaxs
  subroutine psb_dmaxv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_dmaxv

  subroutine psb_dmaxm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_ = dat
      if (info ==0)&
           & call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_max,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_dmaxm
  
  
  subroutine psb_imins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_, dat_
    integer :: iam, np, icomm,info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer,mpi_min,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer,mpi_min,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_imins

  subroutine psb_iminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_=dat
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_iminv

  subroutine psb_iminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_=dat
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_=dat
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_iminm
  
  subroutine psb_smins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_) :: dat_
    integer :: iam, np, icomm,info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_real,mpi_min,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_real,mpi_min,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_smins

  subroutine psb_sminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_real,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_real,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_real,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_sminv

  subroutine psb_sminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_spk_), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_real,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_real,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_real,mpi_min,root_,icomm,info)
      end if
    endif
#endif    

  end subroutine psb_sminm
  
  subroutine psb_dmins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_) :: dat_
    integer :: iam, np, icomm,info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_double_precision,mpi_min,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_double_precision,mpi_min,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_dmins

  subroutine psb_dminv(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

#if !defined(SERIAL_MPI)
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_min,root_,icomm,info)
      end if
    endif
#endif    
  end subroutine psb_dminv

  subroutine psb_dminm(ictxt,dat,root)
    use psb_realloc_mod
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(psb_dpk_), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      dat_ = dat
      if (info ==0) &
           & call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,icomm,info)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_min,root_,icomm,info)
      end if
    endif
#endif    

  end subroutine psb_dminm
  
  
  subroutine psb_iamxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamx2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamxs

  subroutine psb_iamxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamxv

  subroutine psb_iamxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamxm

  subroutine psb_samxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamx2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samxs

  subroutine psb_samxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samxv

  subroutine psb_samxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samxm

  subroutine psb_damxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamx2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damxs

  subroutine psb_damxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damxv

  subroutine psb_damxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damxm


  subroutine psb_camxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamx2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camxs

  subroutine psb_camxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camxv

  subroutine psb_camxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2))) 
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camxm

  subroutine psb_zamxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamx2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamxs

  subroutine psb_zamxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamxv

  subroutine psb_zamxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2))) 
      call gamx2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamx2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamxm


  subroutine psb_iamns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamn2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamns

  subroutine psb_iamnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamnv

  subroutine psb_iamnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_iamnm


  subroutine psb_samns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamn2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samns

  subroutine psb_samnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samnv

  subroutine psb_samnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_samnm

  subroutine psb_damns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamn2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damns

  subroutine psb_damnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damnv

  subroutine psb_damnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_damnm


  subroutine psb_camns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamn2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camns

  subroutine psb_camnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camnv
  
  subroutine psb_camnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2))) 
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_camnm

  subroutine psb_zamns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      call gamn2d(ictxt,'A',dat,ria=ia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamns

  subroutine psb_zamnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia)))
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamnv

  subroutine psb_zamnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
    
#if defined(SERIAL_MPI) 
    if (present(ia)) then 
      ia = 0
    end if
#else
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (present(ia)) then 
      allocate(cia(size(ia,1),size(ia,2))) 
      call gamn2d(ictxt,'A',dat,ria=ia,cia=cia,rrt=root_) 
    else
      call gamn2d(ictxt,'A',dat,rrt=root_) 
    endif
#endif    
  end subroutine psb_zamnm



  subroutine psb_i8sumv(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)                    :: ictxt
    integer(psb_long_int_k_), intent(inout) :: dat(:)
    integer, intent(in), optional          :: root
    integer :: mpi_int8_type, info, icomm
    
    integer                 :: root_, iam, np, isz
    integer(psb_long_int_k_), allocatable  :: dat_(:)
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    mpi_int8_type = mpi_integer8
    isz = size(dat)
    allocate(dat_(isz),stat=info)
    if (root_ == -1) then 
      dat_=dat
      call mpi_allreduce(dat_,dat,isz,mpi_int8_type,mpi_sum,icomm,info)
    else
      if (iam==root_) then 
        dat_=dat
        call mpi_reduce(dat_,dat,isz,mpi_int8_type,mpi_sum,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,isz,mpi_int8_type,mpi_sum,root_,icomm,info)
      end if
    endif

  end subroutine psb_i8sumv


  subroutine psb_i8sums(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)                    :: ictxt
    integer(psb_long_int_k_), intent(inout) :: dat
    integer, intent(in), optional          :: root
    integer :: mpi_int8_type, info, icomm
    
    integer                 :: root_, iam, np
    integer(psb_long_int_k_) :: dat_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    mpi_int8_type = mpi_integer8
    if (root_ == -1) then 
      dat_=dat
      call mpi_allreduce(dat_,dat,1,mpi_int8_type,mpi_sum,icomm,info)
    else
      if (iam==root_) then 
        dat_=dat
        call mpi_reduce(dat_,dat,1,mpi_int8_type,mpi_sum,root_,icomm,info)
      else
        call mpi_reduce(dat,dat_,1,mpi_int8_type,mpi_sum,root_,icomm,info)
      end if
    endif

  end subroutine psb_i8sums

  subroutine psb_i8amx_mpi_user(inv, outv,len,type) 
    integer(psb_long_int_k_) :: inv(*),outv(*)
    integer :: len,type
    integer :: i
    
    do i=1, len
      if (abs(inv(i)) > abs(outv(i))) then 
        outv(i) = inv(i)
      end if
    end do
  end subroutine psb_i8amx_mpi_user
  
  subroutine psb_i8amn_mpi_user(inv, outv,len,type) 
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_long_int_k_) :: inv(*),outv(*)
    integer :: len,type
    integer :: i
    if (type /= mpi_integer8) then 
      write(0,*) 'Invalid type !!!'
    end if
    do i=1, len
      if (abs(inv(i)) < abs(outv(i))) then 
        outv(i) = inv(i)
      end if
    end do
  end subroutine psb_i8amn_mpi_user

  subroutine psb_i8amns(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    integer(psb_long_int_k_) :: dat_
    integer :: iam, np, icomm,info, i8amn
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    
    call mpi_op_create(psb_i8amn_mpi_user,.true.,i8amn,info)
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer8,i8amn,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer8,i8amn,root_,icomm,info)
      dat = dat_
    endif
    call mpi_op_free(i8amn,info)
#endif    
  end subroutine psb_i8amns

  subroutine psb_i8amxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    integer(psb_long_int_k_) :: dat_
    integer :: iam, np, icomm,info, i8amx
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    
    call mpi_op_create(psb_i8amx_mpi_user,.true.,i8amx,info)
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer8,i8amx,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer8,i8amx,root_,icomm,info)
      dat = dat_
    endif
    call mpi_op_free(i8amx,info)
#endif    
  end subroutine psb_i8amxs
  
  subroutine psb_i8mins(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    integer(psb_long_int_k_) :: dat_
    integer :: iam, np, icomm,info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer8,mpi_min,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer8,mpi_min,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_i8mins
  
  subroutine psb_i8maxs(ictxt,dat,root)
#ifdef MPI_MOD
    use mpi
#endif
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer, intent(in)              :: ictxt
    integer(psb_long_int_k_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    integer(psb_long_int_k_) :: dat_
    integer :: iam, np, icomm,info
    
#if !defined(SERIAL_MPI)

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer8,mpi_max,icomm,info)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer8,mpi_max,root_,icomm,info)
      dat = dat_
    endif
#endif    
  end subroutine psb_i8maxs

  subroutine psb_isums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_isums

  subroutine psb_isumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    
  end subroutine psb_isumv

  subroutine psb_isumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1

    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_isumm


  subroutine psb_ssums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_ssums

  subroutine psb_ssumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_ssumv

  subroutine psb_ssumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

#endif    
  end subroutine psb_ssumm

  subroutine psb_dsums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_dsums

  subroutine psb_dsumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_dsumv

  subroutine psb_dsumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

#endif    
  end subroutine psb_dsumm


  subroutine psb_csums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_csums

  subroutine psb_csumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_csumv

  subroutine psb_csumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

#endif    
  end subroutine psb_csumm

  subroutine psb_zsums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_zsums

  subroutine psb_zsumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 
#endif    

  end subroutine psb_zsumv

  subroutine psb_zsumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
#if !defined(SERIAL_MPI)
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

#endif    
  end subroutine psb_zsumm


  subroutine psb_hsnds(ictxt,dat,dst,length)
    use psb_error_mod
    integer, intent(in)           :: ictxt
    character(len=*), intent(in)  :: dat
    integer, intent(in)           :: dst
    integer, intent(in), optional :: length
    integer, allocatable          :: buffer(:)
    integer :: length_, i
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    if (present(length)) then 
      length_ = length
    else
      length_ = len(dat)
    endif
    allocate(buffer(length_))
    do i=1,length_
      buffer(i) = iachar(dat(i:i))
    end do
    
    call gesd2d(ictxt,buffer,dst,0) 
#endif
  end subroutine psb_hsnds

  subroutine psb_hrcvs(ictxt,dat,src,length)
    use psb_error_mod
    integer, intent(in)           :: ictxt
    character(len=*), intent(out)  :: dat
    integer, intent(in)           :: src
    integer, intent(in), optional :: length
    integer, allocatable          :: buffer(:)
    integer :: length_, i
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = ''
#else
    if (present(length)) then 
      length_ = length
    else
      length_ = len(dat)
    endif
    allocate(buffer(length_))
    
    call gerv2d(ictxt,buffer,src,0) 
    do i=1,length_
      dat(i:i) = achar(buffer(i))
    end do
#endif    

  end subroutine psb_hrcvs

  subroutine psb_lsnds(ictxt,dat,dst,length)
    use psb_error_mod
    integer, intent(in)           :: ictxt
    logical, intent(in)           :: dat
    integer, intent(in)           :: dst
    integer :: i
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    if (dat) then 
      i = 1
    else
      i = 0
    endif
    call gesd2d(ictxt,i,dst,0) 
#endif    

  end subroutine psb_lsnds

  subroutine psb_lrcvs(ictxt,dat,src,length)
    use psb_error_mod
    integer, intent(in)           :: ictxt
    logical, intent(out)          :: dat
    integer, intent(in)           :: src
    integer :: i
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = .false.
#else
    call gerv2d(ictxt,i,src,0) 

    dat = (i == 1) 
#endif    

  end subroutine psb_lrcvs


  subroutine psb_isnds(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat
    integer, intent(in)  :: dst
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_isnds

  subroutine psb_isndv(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat(:)
    integer, intent(in)  :: dst

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_isndv

  subroutine psb_isndm(ictxt,dat,dst,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0,m) 
#endif    

  end subroutine psb_isndm


  subroutine psb_ssnds(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_ssnds

  subroutine psb_ssndv(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_ssndv

  subroutine psb_ssndm(ictxt,dat,dst,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m


#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0,m) 
#endif    

  end subroutine psb_ssndm

  subroutine psb_dsnds(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_dsnds

  subroutine psb_dsndv(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_dsndv

  subroutine psb_dsndm(ictxt,dat,dst,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m


#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0,m) 
#endif    

  end subroutine psb_dsndm


  subroutine psb_csnds(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_csnds
  
  subroutine psb_csndv(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_csndv

  subroutine psb_csndm(ictxt,dat,dst,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0,m) 
#endif    

  end subroutine psb_csndm

  subroutine psb_zsnds(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat
    integer, intent(in)  :: dst
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_zsnds
  
  subroutine psb_zsndv(ictxt,dat,dst)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else
    call gesd2d(ictxt,dat,dst,0) 
#endif    

  end subroutine psb_zsndv

  subroutine psb_zsndm(ictxt,dat,dst,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst
    integer, intent(in), optional :: m

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >1) then 
      write(0,*) "Warning: process sending a message in serial mode (to itself)"
    endif
#else

    call gesd2d(ictxt,dat,dst,0,m) 
#endif    

  end subroutine psb_zsndm



  subroutine psb_ircvs(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in)  :: src
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_ircvs

  subroutine psb_ircvv(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in)  :: src

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_ircvv

  subroutine psb_ircvm(ictxt,dat,src,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m


#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else

    call gerv2d(ictxt,dat,src,0,m) 
#endif    

  end subroutine psb_ircvm


  subroutine psb_srcvs(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    integer, intent(in)  :: src
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else

    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_srcvs

  subroutine psb_srcvv(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_srcvv

  subroutine psb_srcvm(ictxt,dat,src,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m


#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0,m) 
#endif    

  end subroutine psb_srcvm

  subroutine psb_drcvs(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    integer, intent(in)  :: src
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else

    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_drcvs

  subroutine psb_drcvv(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_drcvv

  subroutine psb_drcvm(ictxt,dat,src,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m


#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0,m) 
#endif    

  end subroutine psb_drcvm


  subroutine psb_crcvs(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    integer, intent(in)  :: src
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else

    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_crcvs
  
  subroutine psb_crcvv(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_crcvv

  subroutine psb_crcvm(ictxt,dat,src,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0,m) 
#endif    

  end subroutine psb_crcvm

  subroutine psb_zrcvs(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    integer, intent(in)  :: src
    
#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else

    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_zrcvs
  
  subroutine psb_zrcvv(ictxt,dat,src)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0) 
#endif    

  end subroutine psb_zrcvv

  subroutine psb_zrcvm(ictxt,dat,src,m)
    use psb_error_mod
    integer, intent(in)  :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src
    integer, intent(in), optional :: m

#if defined(SERIAL_MPI) 
    if (psb_get_errverbosity() >0) then 
      write(0,*) "Warning: process receiving a message in serial mode (to itself)"
    endif
    dat = 0
#else
    call gerv2d(ictxt,dat,src,0,m) 
#endif    

  end subroutine psb_zrcvm


  subroutine psb_set_coher(ictxt,isvch)
    integer :: ictxt, isvch
    ! Ensure global repeatability for convergence checks.
#if !defined(HAVE_ESSL_BLACS)
    Call blacs_get(ictxt,15,isvch)
    Call blacs_set(ictxt,15,1)
#else
    ! Do nothing: ESSL does coherence by default,
    ! and does not handle req=16 
#endif
  end subroutine psb_set_coher
  subroutine psb_restore_coher(ictxt,isvch)
    integer :: ictxt, isvch
    ! Ensure global coherence for convergence checks.
#if !defined(HAVE_ESSL_BLACS)
    Call blacs_set(ictxt,15,isvch)
#else
    ! Do nothing: ESSL does coherence by default,
    ! and does not handle req=15
#endif
  end subroutine psb_restore_coher
  subroutine psb_get_mpicomm(ictxt,comm)
    integer :: ictxt, comm
#if !defined(SERIAL_MPI)
    call blacs_get(ictxt,10,comm)
#else    
    comm = ictxt
#endif
  end subroutine psb_get_mpicomm
  subroutine psb_get_rank(rank,ictxt,id)
    integer :: rank,ictxt, id
    integer :: blacs_pnum
#if defined(SERIAL_MPI)
    rank = 0
#else
    rank =  blacs_pnum(ictxt,id,0)
#endif    
  end subroutine psb_get_rank
  
#if (!defined(HAVE_KSENDID)) || defined(SERIAL_MPI)
  !
  ! Need these, as they are not in the ESSL implementation 
  ! of the BLACS. 
  !
  integer function krecvid(contxt,proc_to_comm,myrow)
    integer contxt,proc_to_comm,myrow
    
    krecvid=32766
    
    return
  end function krecvid
  integer function ksendid(contxt,proc_to_comm,myrow)
    integer contxt,proc_to_comm,myrow
    
    ksendid=32766
    
    return
  end function ksendid
#endif
  
  
end module psb_penv_mod
