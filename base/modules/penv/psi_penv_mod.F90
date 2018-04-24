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
module psi_penv_mod
  use psb_const_mod
  use psi_comm_buffers_mod, only : psb_buffer_queue

  interface psb_init
    module procedure  psb_init_mpik
  end interface

  interface psb_exit
    module procedure  psb_exit_mpik
  end interface

  interface psb_abort
    module procedure  psb_abort_mpik
  end interface

  interface psb_info
    module procedure psb_info_mpik
  end interface

  interface psb_barrier
    module procedure  psb_barrier_mpik
  end interface
  
  interface psb_init
    module procedure  psb_init_epk
  end interface

  interface psb_exit
    module procedure  psb_exit_epk
  end interface

  interface psb_abort
    module procedure  psb_abort_epk
  end interface

  interface psb_info
    module procedure psb_info_epk
  end interface

  interface psb_barrier
    module procedure  psb_barrier_epk
  end interface

  interface psb_wtime
    module procedure  psb_wtime
  end interface


#if defined(SERIAL_MPI)
  integer(psb_mpk_), private, save :: nctxt=0

#else 

  integer(psb_mpk_), save :: mpi_iamx_op, mpi_iamn_op
  integer(psb_mpk_), save :: mpi_mamx_op, mpi_mamn_op
  integer(psb_mpk_), save :: mpi_eamx_op, mpi_eamn_op
  integer(psb_mpk_), save :: mpi_samx_op, mpi_samn_op
  integer(psb_mpk_), save :: mpi_damx_op, mpi_damn_op
  integer(psb_mpk_), save :: mpi_camx_op, mpi_camn_op
  integer(psb_mpk_), save :: mpi_zamx_op, mpi_zamn_op
  integer(psb_mpk_), save :: mpi_snrm2_op, mpi_dnrm2_op

  type(psb_buffer_queue), save :: psb_mesg_queue 

#endif

  private :: psi_get_sizes,  psi_register_mpi_extras
  private :: psi_iamx_op, psi_iamn_op 
  private :: psi_mamx_op, psi_mamn_op 
  private :: psi_eamx_op, psi_eamn_op 
  private :: psi_samx_op, psi_samn_op 
  private :: psi_damx_op, psi_damn_op 
  private :: psi_camx_op, psi_camn_op 
  private :: psi_zamx_op, psi_zamn_op 
  private :: psi_snrm2_op, psi_dnrm2_op 


contains
  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Environment handling 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!

  subroutine psi_get_sizes()
    use psb_const_mod
    real(psb_dpk_) :: dv(2) 
    real(psb_spk_) :: sv(2) 
    integer(psb_i2pk_):: i2v(2)
    integer(psb_mpk_) :: mv(2)
    integer(psb_ipk_) :: iv(2)
    integer(psb_lpk_) :: lv(2)
    integer(psb_epk_) :: ev(2)

    call psi_c_diffadd(sv(1),sv(2),psb_sizeof_sp)
    call psi_c_diffadd(dv(1),dv(2),psb_sizeof_dp)
    call psi_c_diffadd(i2v(1),i2v(2),psb_sizeof_i2p)
    call psi_c_diffadd(mv(1),mv(2),psb_sizeof_mp)
    call psi_c_diffadd(iv(1),iv(2),psb_sizeof_ip)
    call psi_c_diffadd(lv(1),lv(2),psb_sizeof_lp)
    call psi_c_diffadd(ev(1),ev(2),psb_sizeof_ep)

  end subroutine psi_get_sizes

  subroutine  psi_register_mpi_extras(info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_) :: info
    
    info = 0
#if 0
    if (info == 0) call mpi_type_create_f90_integer(psb_ipk_, psb_mpi_ipk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_lpk_, psb_mpi_lpk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_mpk_, psb_mpi_mpk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_epk_, psb_mpi_lpk_ ,info)
    if (info == 0) call mpi_type_create_f90_real(psb_spk_p_,psb_spk_r_, psb_mpi_r_spk_,info)
    if (info == 0) call mpi_type_create_f90_real(psb_dpk_p_,psb_dpk_r_, psb_mpi_r_dpk_,info)
    if (info == 0) call mpi_type_create_f90_complex(psb_spk_p_,psb_spk_r_, psb_mpi_c_spk_,info)
    if (info == 0) call mpi_type_create_f90_complex(psb_dpk_p_,psb_dpk_r_, psb_mpi_c_dpk_,info)
#else
#if defined(IPK4) && defined(LPK4)
    psb_mpi_ipk_ = mpi_integer4
    psb_mpi_lpk_ = mpi_integer4
#elif defined(IPK4) && defined(LPK8)
    psb_mpi_ipk_ = mpi_integer4
    psb_mpi_lpk_ = mpi_integer8
#elif defined(IPK8) && defined(LPK8)
    psb_mpi_ipk_ = mpi_integer8
    psb_mpi_lpk_ = mpi_integer8
#else
    ! This should never happen
    write(psb_err_unit,*) 'Warning: an impossible IPK/LPK combination.'
    write(psb_err_unit,*) 'Something went wrong at configuration time.'
    psb_mpi_ipk_ = -1
    psb_mpi_lpk_ = -1
#endif
    psb_mpi_i2pk_ = mpi_integer2
    psb_mpi_mpk_  = mpi_integer4
    psb_mpi_epk_  = mpi_integer8
    psb_mpi_r_spk_  = mpi_real
    psb_mpi_r_dpk_  = mpi_double_precision
    psb_mpi_c_spk_  = mpi_complex
    psb_mpi_c_dpk_  = mpi_double_complex
#endif

#if defined(SERIAL_MPI)
#else 
    if (info == 0) call mpi_op_create(psi_mamx_op,.true.,mpi_mamx_op,info)
    if (info == 0) call mpi_op_create(psi_mamn_op,.true.,mpi_mamn_op,info)
    if (info == 0) call mpi_op_create(psi_eamx_op,.true.,mpi_eamx_op,info)
    if (info == 0) call mpi_op_create(psi_eamn_op,.true.,mpi_eamn_op,info)
    if (info == 0) call mpi_op_create(psi_samx_op,.true.,mpi_samx_op,info)
    if (info == 0) call mpi_op_create(psi_samn_op,.true.,mpi_samn_op,info)
    if (info == 0) call mpi_op_create(psi_damx_op,.true.,mpi_damx_op,info)
    if (info == 0) call mpi_op_create(psi_damn_op,.true.,mpi_damn_op,info)
    if (info == 0) call mpi_op_create(psi_camx_op,.true.,mpi_camx_op,info)
    if (info == 0) call mpi_op_create(psi_camn_op,.true.,mpi_camn_op,info)
    if (info == 0) call mpi_op_create(psi_zamx_op,.true.,mpi_zamx_op,info)
    if (info == 0) call mpi_op_create(psi_zamn_op,.true.,mpi_zamn_op,info)
    if (info == 0) call mpi_op_create(psi_snrm2_op,.true.,mpi_snrm2_op,info)
    if (info == 0) call mpi_op_create(psi_dnrm2_op,.true.,mpi_dnrm2_op,info)
#endif

  end subroutine psi_register_mpi_extras

  subroutine psb_init_epk(ictxt,np,basectxt,ids)
    integer(psb_epk_), intent(out) :: ictxt
    integer(psb_epk_), intent(in), optional :: np, basectxt, ids(:)

    integer(psb_mpk_) :: iictxt
    integer(psb_mpk_) :: inp, ibasectxt
    integer(psb_mpk_), allocatable :: ids_(:)

    if (present(ids)) then 
      allocate(ids_(size(ids)))
      ids_ = ids
    else
      allocate(ids_(0))
    end if
    if (present(np).and.present(basectxt)) then 
      inp       = np
      ibasectxt = basectxt
      call psb_init(iictxt,np=inp,basectxt=ibasectxt,ids=ids_)
    else if (present(np)) then 
      inp       = np
      call psb_init(iictxt,np=inp,ids=ids_)
    else if (present(basectxt)) then 
      ibasectxt = basectxt
      call psb_init(iictxt,basectxt=ibasectxt,ids=ids_)
    else
      call psb_init(iictxt,ids=ids_)
    end if
    ictxt = iictxt
  end subroutine psb_init_epk

  subroutine psb_exit_epk(ictxt,close)
    integer(psb_epk_), intent(inout) :: ictxt
    logical, intent(in), optional :: close
    integer(psb_mpk_) :: iictxt
    
    iictxt = ictxt
    call psb_exit(iictxt, close)
  end subroutine psb_exit_epk

  subroutine psb_barrier_epk(ictxt)
    integer(psb_epk_), intent(in) :: ictxt
    integer(psb_mpk_) :: iictxt
    
    iictxt = ictxt
    call psb_barrier(iictxt)
  end subroutine psb_barrier_epk

  subroutine psb_abort_epk(ictxt,errc)
    integer(psb_epk_), intent(in) :: ictxt
    integer(psb_epk_), intent(in), optional :: errc
    integer(psb_mpk_) :: iictxt, ierrc

    iictxt = ictxt
    if (present(errc)) then 
      ierrc = errc
      call psb_abort(iictxt,ierrc)
    else
      call psb_abort(iictxt)
    end if
  end subroutine psb_abort_epk
  
  subroutine psb_info_epk(ictxt,iam,np)

    integer(psb_epk_), intent(in)  :: ictxt
    integer(psb_epk_), intent(out) :: iam, np

    !
    ! Simple caching scheme, keep track
    ! of the last CTXT encountered. 
    !
    integer(psb_mpk_), save :: lctxt=-1, lam, lnp
    if (ictxt /= lctxt) then
      lctxt = ictxt
      call psb_info(lctxt,lam,lnp)
    end if
    iam = lam
    np  = lnp
  end subroutine psb_info_epk
  
  subroutine psb_init_mpik(ictxt,np,basectxt,ids)
    use psi_comm_buffers_mod 
    use psb_const_mod
    use psb_error_mod
    use psb_mat_mod
    use psb_vect_mod
! !$    use psb_rsb_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(out) :: ictxt
    integer(psb_mpk_), intent(in), optional :: np, basectxt, ids(:)

    integer(psb_mpk_) :: i, isnullcomm
    integer(psb_mpk_), allocatable :: iids(:) 
    logical :: initialized    
    integer(psb_mpk_) :: np_, npavail, iam, info, basecomm, basegroup, newgroup
    character(len=20), parameter :: name='psb_init'
    integer(psb_ipk_) :: iinfo
    !    
    call psb_set_debug_unit(psb_err_unit)

#if defined(SERIAL_MPI) 
    ictxt = nctxt
    nctxt = nctxt + 1

    call psi_register_mpi_extras(info)
    call psi_get_sizes()

#else    
    call mpi_initialized(initialized,info)
    if ((.not.initialized).or.(info /= mpi_success)) then 
      if (info == mpi_success) call mpi_init(info) 
      if (info /= mpi_success) then
        write(psb_err_unit,*) 'Error in initalizing MPI, bailing out',info 
        stop 
      end if
    end if

    if (present(basectxt)) then 
      basecomm = basectxt
    else
      basecomm = mpi_comm_world
    end if

    if (present(np)) then 
      if (np < 1) then 
        iinfo=psb_err_initerror_neugh_procs_
        call psb_errpush(iinfo,name)
        call psb_error()
        ictxt = mpi_comm_null
        return
      endif
      call mpi_comm_size(basecomm,np_,info)
      if (np_ < np) then 
        iinfo=psb_err_initerror_neugh_procs_
        call psb_errpush(iinfo,name)
        call psb_error()
        ictxt = mpi_comm_null
        return
      endif
      call mpi_comm_group(basecomm,basegroup,info)
      if (present(ids)) then 
        if (size(ids)<np) then 
          write(psb_err_unit,*) 'Error in init: too few ids in input'
          ictxt = mpi_comm_null
          return
        end if
        do i=1, np 
          if ((ids(i)<0).or.(ids(i)>np_)) then 
            write(psb_err_unit,*) 'Error in init: invalid rank in input'
            ictxt = mpi_comm_null
            return
          end if
        end do
        call mpi_group_incl(basegroup,np,ids,newgroup,info)
        if (info /= mpi_success) then 
          ictxt = mpi_comm_null 
          return
        endif
      else
        allocate(iids(np),stat=info)
        if (info /= 0) then 
          ictxt = mpi_comm_null
          return
        endif
        do i=1, np
          iids(i) = i-1
        end do
        call mpi_group_incl(basegroup,np,iids,newgroup,info)
        if (info /= mpi_success) then 
          ictxt = mpi_comm_null 
          return
        endif
        deallocate(iids)
      end if
      call mpi_comm_create(basecomm,newgroup,ictxt,info)

    else
      if (basecomm /= mpi_comm_null) then 
        call mpi_comm_dup(basecomm,ictxt,info)
      else 
        ictxt = mpi_comm_null
      end if
    endif
    call psi_register_mpi_extras(info)
    call psi_get_sizes()
    if (ictxt == mpi_comm_null) return 
#endif
    call psb_init_vect_defaults()
    call psb_init_mat_defaults()
    ! !$    call psb_rsb_init(info)
    ! !$    if (info.ne.psb_rsb_const_success) then 
    ! !$      if (info.eq.psb_rsb_const_not_available) then 
    ! !$        info=psb_success_ ! rsb is not present
    ! !$      else
    ! !$        ! rsb failed to initialize, and we issue an internal error.
    ! !$        ! or shall we tolerate this ?
    ! !$        info=psb_err_internal_error_
    ! !$        call psb_errpush(info,name)
    ! !$        call psb_error(ictxt)
    ! !$      endif
    ! !$    endif

  end subroutine psb_init_mpik

  subroutine psb_exit_mpik(ictxt,close)
    use psi_comm_buffers_mod 
! !$    use psb_rsb_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(inout) :: ictxt
    logical, intent(in), optional :: close
    logical  :: close_
    integer(psb_mpk_) :: info
    character(len=20), parameter :: name='psb_exit'

    info = 0
    if (present(close)) then 
      close_ = close
    else
      close_ = .true.
    end if
! !$    if (close_) call psb_rsb_exit(info)
! !$    if (info.ne.psb_rsb_const_success) then 
! !$      if (info.eq.psb_rsb_const_not_available) then 
! !$        info=psb_success_ ! rsb is not present
! !$      else
! !$        info=psb_err_internal_error_ ! rsb failed to exit, and we issue an internal error. or  shall we tolerate this ?
! !$        call psb_errpush(info,name)
! !$        call psb_error(ictxt)
! !$      endif
! !$    endif
#if defined(SERIAL_MPI)
    ! Under serial mode, CLOSE has no effect, but reclaim
    ! the used ICTXT number. 
    nctxt = max(0, nctxt - 1)    
#else 
    if (close_) then 
      call psb_close_all_context(psb_mesg_queue)
    else
      call psb_close_context(psb_mesg_queue,ictxt)
    end if
    if ((ictxt /= mpi_comm_null).and.(ictxt /= mpi_comm_world)) then 
      call mpi_comm_Free(ictxt,info)
    end if

    if (close_) call mpi_finalize(info)
#endif


  end subroutine psb_exit_mpik


  subroutine psb_barrier_mpik(ictxt)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_), intent(in) :: ictxt

    integer(psb_mpk_) :: info
#if !defined(SERIAL_MPI)
    if (ictxt /= mpi_comm_null) then 
      call mpi_barrier(ictxt, info)
    end if
#endif    

  end subroutine psb_barrier_mpik

  function psb_wtime()
    use psb_const_mod
!    use mpi_constants
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    real(psb_dpk_) :: psb_wtime

    psb_wtime = mpi_wtime()
  end function psb_wtime

  subroutine psb_abort_mpik(ictxt,errc)
    use psi_comm_buffers_mod 

    integer(psb_mpk_), intent(in) :: ictxt
    integer(psb_mpk_), intent(in), optional :: errc
    
    integer(psb_mpk_) :: code, info 

#if defined(SERIAL_MPI) 
    stop 
#else    
    if (present(errc)) then 
      code = errc
    else
      code = -1 
    endif

    call mpi_abort(ictxt,code,info)
#endif    

  end subroutine psb_abort_mpik


  subroutine psb_info_mpik(ictxt,iam,np)
    use psi_comm_buffers_mod 
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_mpk_), intent(out) :: iam, np
    integer(psb_mpk_) :: info
    !
    ! Simple caching scheme, keep track
    ! of the last CTXT encountered. 
    !  
    integer(psb_mpk_), save :: lctxt=-1, lam, lnp
    
#if defined(SERIAL_MPI) 
    iam = 0
    np  = 1
#else    
    iam = -1
    np  = -1
    if (ictxt == lctxt) then
      iam = lam
      np  = lnp
    else
      if (ictxt /= mpi_comm_null) then 
        call mpi_comm_size(ictxt,np,info) 
        if (info /= mpi_success) np = -1 
        call mpi_comm_rank(ictxt,iam,info) 
        if (info /= mpi_success) iam = -1 
      end if
      lctxt = ictxt
      lam   = iam
      lnp   = np
    end if
#endif    

  end subroutine psb_info_mpik


  subroutine psb_get_mpicomm(ictxt,comm)
    integer(psb_mpk_) :: ictxt, comm

    comm = ictxt
  end subroutine psb_get_mpicomm

  subroutine psb_get_rank(rank,ictxt,id)
    integer(psb_mpk_) :: rank,ictxt,id

    rank = id
  end subroutine psb_get_rank


  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Base binary  operations
  !
  ! Note: len & type are always default integer.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!
  subroutine psi_mamx_op(inv, outv,len,type) 
    integer(psb_mpk_) :: inv(*),outv(*)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) > abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_mamx_op

  subroutine psi_mamn_op(inv, outv,len,type) 
    integer(psb_mpk_) :: inv(*),outv(*)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) < abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_mamn_op

  subroutine psi_eamx_op(inv, outv,len,type) 
    integer(psb_epk_) :: inv(*),outv(*)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) > abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_eamx_op

  subroutine psi_eamn_op(inv, outv,len,type) 
    integer(psb_epk_) :: inv(*),outv(*)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) < abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_eamn_op

  subroutine psi_samx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_samx_op

  subroutine psi_samn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_samn_op

  subroutine psi_damx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_damx_op

  subroutine psi_damn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_damn_op

  subroutine psi_camx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_spk_), intent(in)    :: vin(len)
    complex(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_camx_op

  subroutine psi_camn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_spk_), intent(in)    :: vin(len)
    complex(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_camn_op

  subroutine psi_zamx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_dpk_), intent(in)    :: vin(len)
    complex(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_zamx_op

  subroutine psi_zamn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_dpk_), intent(in)    :: vin(len)
    complex(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_zamn_op

  subroutine psi_snrm2_op(vin,vinout,len,itype)
    implicit none 
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    real(psb_spk_) :: w, z
    do i=1, len
      w = max( vin(i), vinout(i) )
      z = min( vin(i), vinout(i) )
      if ( z == szero ) then
        vinout(i) = w
      else
        vinout(i) = w*sqrt( sone+( z / w )**2 )
      end if
    end do
  end subroutine psi_snrm2_op

  subroutine psi_dnrm2_op(vin,vinout,len,itype)
    implicit none 
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    real(psb_dpk_) :: w, z
    do i=1, len
      w = max( vin(i), vinout(i) )
      z = min( vin(i), vinout(i) )
      if ( z == dzero ) then
        vinout(i) = w
      else
        vinout(i) = w*sqrt( done+( z / w )**2 )
      end if
    end do
  end subroutine psi_dnrm2_op

end module psi_penv_mod
