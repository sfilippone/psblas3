

module psi_bcast_mod
  use psb_const_mod
  use psi_penv_mod
  interface psb_bcast
    module procedure psb_ibcasts, psb_ibcastv, psb_ibcastm,&
         & psb_dbcasts, psb_dbcastv, psb_dbcastm,&
         & psb_zbcasts, psb_zbcastv, psb_zbcastm,&
         & psb_sbcasts, psb_sbcastv, psb_sbcastm,&
         & psb_cbcasts, psb_cbcastv, psb_cbcastm,&
         & psb_hbcasts, psb_hbcastv, psb_lbcasts, psb_lbcastv
  end interface

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
    integer(psb_ipk_), intent(in)      :: ictxt
    integer(psb_ipk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,psb_mpi_integer,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional  :: root

    integer(psb_ipk_) :: iam, np, root_,  info
#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_integer,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    integer(psb_ipk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),psb_mpi_integer,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)      :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,mpi_real,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_real,root_,ictxt,info)

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
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_real,root_,ictxt,info)

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
    integer(psb_ipk_), intent(in)      :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,mpi_double_precision,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_double_precision,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    real(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_double_precision,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)      :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,mpi_complex,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_complex,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_spk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_complex,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)      :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,1,mpi_double_complex,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_double_complex,root_,ictxt,info) 
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
    integer(psb_ipk_), intent(in)    :: ictxt
    complex(psb_dpk_), intent(inout) :: dat(:,:)
    integer(psb_ipk_), intent(in), optional :: root

    integer(psb_ipk_) :: iam, np, root_,  info

#if !defined(SERIAL_MPI)
    if (present(root)) then
      root_ = root
    else
      root_ = psb_root_
    endif

    call psb_info(ictxt,iam,np)
    call mpi_bcast(dat,size(dat),mpi_double_complex,root_,ictxt,info)
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
    integer(psb_ipk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat
    integer(psb_ipk_), intent(in), optional   :: root,length

    integer(psb_ipk_) :: iam, np, root_,length_,info

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
    integer(psb_ipk_), intent(in)             :: ictxt
    character(len=*), intent(inout) :: dat(:)
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_ipk_) :: iam, np, root_,length_,info, size_

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
    integer(psb_ipk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_ipk_) :: iam, np, root_,info

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
    integer(psb_ipk_), intent(in)             :: ictxt
    logical, intent(inout)          :: dat(:)
    integer(psb_ipk_), intent(in), optional   :: root

    integer(psb_ipk_) :: iam, np, root_,info

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

end module psi_bcast_mod
