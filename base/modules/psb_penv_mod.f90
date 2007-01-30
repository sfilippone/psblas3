!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
module psb_penv_mod


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
         & psb_hbcasts, psb_lbcasts
  end interface


  interface psb_snd
    module procedure psb_isnds, psb_isndv, psb_isndm,&
         & psb_dsnds, psb_dsndv, psb_dsndm,&
         & psb_zsnds, psb_zsndv, psb_zsndm,&
         & psb_hsnds, psb_lsnds
  end interface

  interface psb_rcv
    module procedure psb_ircvs, psb_ircvv, psb_ircvm,&
         & psb_drcvs, psb_drcvv, psb_drcvm,&
         & psb_zrcvs, psb_zrcvv, psb_zrcvm,&
         & psb_hrcvs, psb_lrcvs
  end interface

  interface psb_max
    module procedure psb_imaxs, psb_imaxv, psb_imaxm,&
         & psb_dmaxs, psb_dmaxv, psb_dmaxm
  end interface


  interface psb_min
    module procedure psb_imins, psb_iminv, psb_iminm,&
         & psb_dmins, psb_dminv, psb_dminm
  end interface


  interface psb_amx
    module procedure psb_iamxs, psb_iamxv, psb_iamxm,&
         & psb_damxs, psb_damxv, psb_damxm,&
         & psb_zamxs, psb_zamxv, psb_zamxm
  end interface

  interface psb_amn
    module procedure psb_iamns, psb_iamnv, psb_iamnm,&
         & psb_damns, psb_damnv, psb_damnm,&
         & psb_zamns, psb_zamnv, psb_zamnm
  end interface

  interface psb_sum
    module procedure psb_isums, psb_isumv, psb_isumm,&
         & psb_dsums, psb_dsumv, psb_dsumm,&
         & psb_zsums, psb_zsumv, psb_zsumm
  end interface

  
  interface gebs2d
    module procedure igebs2ds, igebs2dv, igebs2dm,&
         &           dgebs2ds, dgebs2dv, dgebs2dm,&
         &           zgebs2ds, zgebs2dv, zgebs2dm
  end interface

  interface gebr2d
    module procedure igebr2ds, igebr2dv, igebr2dm,&
         &           dgebr2ds, dgebr2dv, dgebr2dm,&
         &           zgebr2ds, zgebr2dv, zgebr2dm    
  end interface

  interface gesd2d
    module procedure igesd2ds, igesd2dv, igesd2dm,&
         &           dgesd2ds, dgesd2dv, dgesd2dm,&
         &           zgesd2ds, zgesd2dv, zgesd2dm
  end interface

  interface gerv2d
    module procedure igerv2ds, igerv2dv, igerv2dm,&
         &           dgerv2ds, dgerv2dv, dgerv2dm,&
         &           zgerv2ds, zgerv2dv, zgerv2dm
  end interface

  interface gsum2d
    module procedure igsum2ds, igsum2dv, igsum2dm,&
         &           dgsum2ds, dgsum2dv, dgsum2dm,&
         &           zgsum2ds, zgsum2dv, zgsum2dm
  end interface

  interface gamx2d
    module procedure igamx2ds, igamx2dv, igamx2dm,&
         &           dgamx2ds, dgamx2dv, dgamx2dm,&
         &           zgamx2ds, zgamx2dv, zgamx2dm
  end interface


  interface gamn2d
    module procedure igamn2ds, igamn2dv, igamn2dm,&
         &           dgamn2ds, dgamn2dv, dgamn2dm,&
         &           zgamn2ds, zgamn2dv, zgamn2dm
  end interface


contains 


  

  subroutine psb_init(ictxt,np)
    use psb_const_mod
    use psb_error_mod
    integer, intent(out) :: ictxt
    integer, intent(in), optional :: np
    
    integer :: np_, npavail, iam, info
    character(len=20), parameter :: name='psb_init'
    
    call blacs_pinfo(iam, npavail)
    call blacs_get(izero, izero, ictxt)
    
    if (present(np)) then 
      np_ = max(1,min(np,npavail))
    else
      np_ = npavail
    endif
    
    call blacs_gridinit(ictxt, 'R', np_, ione)
    
    if (present(np)) then 
      if (np_ < np) then 
        info = 2011
        call psb_errpush(info,name)
        call psb_error(ictxt)
      endif
    endif
        
  end subroutine psb_init

  subroutine psb_exit(ictxt,close)
    integer, intent(in) :: ictxt
    logical, intent(in), optional :: close
    logical  :: close_
    integer  :: nprow, npcol, myprow, mypcol
    
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
  end subroutine psb_exit


  subroutine psb_barrier(ictxt)
    integer, intent(in) :: ictxt
    
    call blacs_barrier(ictxt,'All')
    
  end subroutine psb_barrier

  function psb_wtime()
    real(kind(1.d0)) :: psb_wtime
    
    real(kind(1.d0)), external :: mpi_wtime
    
    psb_wtime = mpi_wtime()
  end function psb_wtime

  subroutine psb_abort(ictxt)
    integer, intent(in) :: ictxt
    
    call blacs_abort(ictxt,-1)
    
  end subroutine psb_abort


  subroutine psb_info(ictxt,iam,np)

    integer, intent(in)  :: ictxt
    integer, intent(out) :: iam, np
    integer              :: nprow, npcol, myprow, mypcol
    
    call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
    
    iam = myprow
    np  = nprow
    
  end subroutine psb_info

    
  subroutine psb_ibcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    integer, intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_ibcasts

  subroutine psb_ibcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    integer, intent(inout) :: dat(:)
    integer, intent(in), optional  :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_ibcastv
    
  subroutine psb_ibcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    integer, intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_ibcastm


  subroutine psb_dbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_dbcasts


  subroutine psb_dbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(kind(1.d0)), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_dbcastv
    
  subroutine psb_dbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    real(kind(1.d0)), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_dbcastm
    

  subroutine psb_zbcasts(ictxt,dat,root)
    integer, intent(in)      :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_zbcasts

  subroutine psb_zbcastv(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(kind(1.d0)), intent(inout) :: dat(:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_zbcastv
    
  subroutine psb_zbcastm(ictxt,dat,root)
    integer, intent(in)    :: ictxt
    complex(kind(1.d0)), intent(inout) :: dat(:,:)
    integer, intent(in), optional :: root

    integer  :: iam, np, root_

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    
    if (iam==root_) then 
      call gebs2d(ictxt,'A',dat)
    else
      call gebr2d(ictxt,'A',dat,rrt=root_)
    endif
  end subroutine psb_zbcastm


  subroutine psb_hbcasts(ictxt,dat,root,length)
    use mpi
    integer, intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer, intent(in), optional   :: root,length

    integer  :: iam, np, root_,icomm,length_,info

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif
    if (present(length)) then
      length_ = length
    else
      length_ = len(dat)
    endif

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    
    call mpi_bcast(dat,length_,MPI_CHARACTER,root_,icomm,info)

  end subroutine psb_hbcasts

  subroutine psb_lbcasts(ictxt,dat,root)
    use mpi
    integer, intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer, intent(in), optional   :: root

    integer  :: iam, np, root_,icomm,info

    if (present(root)) then
      root_ = root
    else
      root_ = 0
    endif

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)
    call mpi_bcast(dat,1,MPI_LOGICAL,root_,icomm,info)

  end subroutine psb_lbcasts



  subroutine psb_imaxs(ictxt,dat,root)
    use mpi
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_, dat_
    integer :: iam, np, icomm
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer,mpi_max,icomm)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer,mpi_max,root_,icomm)
      dat = dat_
    endif
  end subroutine psb_imaxs
  subroutine psb_imaxv(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_max,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_max,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_max,root_,icomm)
      end if
    endif
  end subroutine psb_imaxv
  subroutine psb_imaxm(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_max,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_max,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_max,root_,icomm)
      end if
    endif
  end subroutine psb_imaxm
  
  subroutine psb_dmaxs(ictxt,dat,root)
    use mpi
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)) :: dat_
    integer :: iam, np, icomm
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_double_precision,mpi_max,icomm)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_double_precision,mpi_max,root_,icomm)
      dat = dat_
    endif
  end subroutine psb_dmaxs
  subroutine psb_dmaxv(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

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
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_max,root_,icomm)
      end if
    endif
  end subroutine psb_dmaxv
  subroutine psb_dmaxm(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

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
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_max,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_max,root_,icomm)
      end if
    endif
  end subroutine psb_dmaxm
  
  
  subroutine psb_imins(ictxt,dat,root)
    use mpi
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_, dat_
    integer :: iam, np, icomm
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_integer,mpi_min,icomm)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_integer,mpi_min,root_,icomm)
      dat = dat_
    endif
  end subroutine psb_imins
  subroutine psb_iminv(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat),dat_,info)
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_min,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_min,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_min,root_,icomm)
      end if
    endif
  end subroutine psb_iminv
  subroutine psb_iminm(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    integer, allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call psb_realloc(size(dat,1),size(dat,2),dat_,info)
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_integer,mpi_min,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        call mpi_reduce(dat_,dat,size(dat),mpi_integer,mpi_min,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_integer,mpi_min,root_,icomm)
      end if
    endif
  end subroutine psb_iminm
  
  subroutine psb_dmins(ictxt,dat,root)
    use mpi
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)) :: dat_
    integer :: iam, np, icomm
    

    call psb_info(ictxt,iam,np)
    call psb_get_mpicomm(ictxt,icomm)

    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    if (root_ == -1) then 
      call mpi_allreduce(dat,dat_,1,mpi_double_precision,mpi_min,icomm)
      dat = dat_
    else
      call mpi_reduce(dat,dat_,1,mpi_double_precision,mpi_min,root_,icomm)
      dat = dat_
    endif
  end subroutine psb_dmins
  subroutine psb_dminv(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)), allocatable :: dat_(:)
    integer :: iam, np, icomm, info
    

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
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_min,root_,icomm)
      end if
    endif
  end subroutine psb_dminv
  subroutine psb_dminm(ictxt,dat,root)
    use mpi
    use psb_realloc_mod
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer :: root_
    real(kind(1.d0)), allocatable :: dat_(:,:)
    integer :: iam, np, icomm, info
    

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
      if (info ==0) call mpi_allreduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,icomm)
    else
      if (iam==root_) then 
        call psb_realloc(size(dat,1),size(dat,2),dat_,info)
        dat_ = dat
        call mpi_reduce(dat_,dat,size(dat),mpi_double_precision,mpi_min,root_,icomm)
      else
        call mpi_reduce(dat,dat_,size(dat),mpi_double_precision,mpi_min,root_,icomm)
      end if
    endif
  end subroutine psb_dminm
  
  
  subroutine psb_iamxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_iamxs

  subroutine psb_iamxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_iamxv

  subroutine psb_iamxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_iamxm


  subroutine psb_damxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_damxs

  subroutine psb_damxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_damxv

  subroutine psb_damxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_damxm


  subroutine psb_zamxs(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_zamxs

  subroutine psb_zamxv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_zamxv

  subroutine psb_zamxm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_zamxm



  subroutine psb_iamns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_iamns

  subroutine psb_iamnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_iamnv

  subroutine psb_iamnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_iamnm


  subroutine psb_damns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_damns

  subroutine psb_damnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_damnv

  subroutine psb_damnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_damnm


  subroutine psb_zamns(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia
    
    integer   :: root_
    
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
  end subroutine psb_zamns

  subroutine psb_zamnv(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:)
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
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
  end subroutine psb_zamnv

  subroutine psb_zamnm(ictxt,dat,root,ia)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    integer, intent(inout), optional :: ia(:,:)
    
    integer   :: root_
    integer, allocatable :: cia(:,:)
    
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
  end subroutine psb_zamnm



  subroutine psb_isums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_isums

  subroutine psb_isumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 
  end subroutine psb_isumv

  subroutine psb_isumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1

    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_isumm


  subroutine psb_dsums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_dsums

  subroutine psb_dsumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_dsumv

  subroutine psb_dsumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_dsumm


  subroutine psb_zsums(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_zsums

  subroutine psb_zsumv(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    integer, allocatable :: cia(:)
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif

    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_zsumv

  subroutine psb_zsumm(ictxt,dat,root)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in), optional    :: root
    
    integer   :: root_
    
    if (present(root)) then 
      root_ = root
    else
      root_ = -1
    endif
    
    call gsum2d(ictxt,'A',dat,rrt=root_) 

  end subroutine psb_zsumm



  subroutine psb_hsnds(ictxt,dat,dst,length)
    integer, intent(in)           :: ictxt
    character(len=*), intent(in)  :: dat
    integer, intent(in)           :: dst
    integer, intent(in), optional :: length
    integer, allocatable          :: buffer(:)
    integer :: length_, i
    
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

  end subroutine psb_hsnds

  subroutine psb_hrcvs(ictxt,dat,src,length)
    integer, intent(in)           :: ictxt
    character(len=*), intent(out)  :: dat
    integer, intent(in)           :: src
    integer, intent(in), optional :: length
    integer, allocatable          :: buffer(:)
    integer :: length_, i
    
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

  end subroutine psb_hrcvs

  subroutine psb_lsnds(ictxt,dat,dst,length)
    integer, intent(in)           :: ictxt
    logical, intent(in)           :: dat
    integer, intent(in)           :: dst
    integer :: i
    
    if (dat) then 
      i = 1
    else
      i = 0
    endif
    call gesd2d(ictxt,i,dst,0) 

  end subroutine psb_lsnds

  subroutine psb_lrcvs(ictxt,dat,src,length)
    integer, intent(in)           :: ictxt
    logical, intent(out)          :: dat
    integer, intent(in)           :: src
    integer :: i
    
    call gerv2d(ictxt,i,src,0) 

    dat = (i == 1) 

  end subroutine psb_lrcvs


  subroutine psb_isnds(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat
    integer, intent(in)  :: dst
    

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_isnds

  subroutine psb_isndv(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat(:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_isndv

  subroutine psb_isndm(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    integer, intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_isndm


  subroutine psb_dsnds(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(in)  :: dat
    integer, intent(in)  :: dst
    

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_dsnds

  subroutine psb_dsndv(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_dsndv

  subroutine psb_dsndm(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_dsndm


  subroutine psb_zsnds(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(in)  :: dat
    integer, intent(in)  :: dst
    

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_zsnds
  
  subroutine psb_zsndv(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(in)  :: dat(:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_zsndv

  subroutine psb_zsndm(ictxt,dat,dst)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(in)  :: dat(:,:)
    integer, intent(in)  :: dst

    call gesd2d(ictxt,dat,dst,0) 

  end subroutine psb_zsndm



  subroutine psb_ircvs(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat
    integer, intent(in)  :: src
    

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_ircvs

  subroutine psb_ircvv(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat(:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_ircvv

  subroutine psb_ircvm(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    integer, intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_ircvm


  subroutine psb_drcvs(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in)  :: src
    

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_drcvs

  subroutine psb_drcvv(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_drcvv

  subroutine psb_drcvm(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_drcvm


  subroutine psb_zrcvs(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    integer, intent(in)  :: src
    

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_zrcvs
  
  subroutine psb_zrcvv(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_zrcvv

  subroutine psb_zrcvm(ictxt,dat,src)
    integer, intent(in)  :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    integer, intent(in)  :: src

    call gerv2d(ictxt,dat,src,0) 

  end subroutine psb_zrcvm







  
  subroutine igebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    character             :: top_ 
    
    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine igebs2ds

  subroutine igebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    

    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(ictxt,scope,top_,size(dat,1),1,dat,size(dat,1))
    
  end subroutine igebs2dv

  subroutine igebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat(:,:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine igebs2dm



  subroutine dgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine dgebs2ds

  subroutine dgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine dgebs2dv

  subroutine dgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine dgebs2dm



  subroutine zgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine zgebs2ds

  subroutine zgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine zgebs2dv

  subroutine zgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine zgebs2dm





  subroutine dgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol
    
    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call dgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgebr2ds

  subroutine dgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgebr2dv

  subroutine dgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgebr2dm




  subroutine zgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call zgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgebr2ds

  subroutine zgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgebr2dv

  subroutine zgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgebr2dm



  subroutine igebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call igebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igebr2ds

  subroutine igebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igebr2dv

  subroutine igebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igebr2dm



  subroutine dgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgesd2ds


  subroutine dgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgesd2dv

  subroutine dgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call dgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgesd2dm


  subroutine igesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igesd2ds


  subroutine igesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igesd2dv

  subroutine igesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call igesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igesd2dm



  subroutine zgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgesd2ds


  subroutine zgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgesd2dv

  subroutine zgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call zgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgesd2dm



  subroutine dgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgerv2ds


  subroutine dgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgerv2dv

  subroutine dgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call dgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgerv2dm


  subroutine igerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igerv2ds


  subroutine igerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igerv2dv

  subroutine igerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface


    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call igerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igerv2dm



  subroutine zgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgerv2ds


  subroutine zgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgerv2dv

  subroutine zgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface


    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call zgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgerv2dm



  subroutine dgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call dgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgsum2ds

  subroutine dgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgsum2dv

  subroutine dgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgsum2dm



  subroutine igsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call igsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igsum2ds

  subroutine igsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igsum2dv

  subroutine igsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igsum2dm



  subroutine zgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call zgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgsum2ds

  subroutine zgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgsum2dv

  subroutine zgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgsum2dm




  subroutine dgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call dgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamx2ds


  subroutine dgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).and.present(cia)) then 
      call dgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamx2dv

  subroutine dgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamx2dm



  subroutine igamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        integer, intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call igamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamx2ds


  subroutine igamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        integer, intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamx2dv

  subroutine igamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        integer, intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine igamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamx2dm

  

  subroutine zgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call zgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamx2ds


  subroutine zgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamx2dv

  subroutine zgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamx2dm
  

  subroutine dgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call dgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamn2ds


  subroutine dgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamn2dv

  subroutine dgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamn2dm



  subroutine igamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        integer, intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call igamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamn2ds


  subroutine igamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        integer, intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamn2dv

  subroutine igamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        integer, intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine igamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamn2dm

  

  subroutine zgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call zgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamn2ds


  subroutine zgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamn2dv

  subroutine zgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamn2dm
  

end module psb_penv_mod
