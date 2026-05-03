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
!         software without specific prior written permission.
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
!
! package: psb_l_base_vect_mod
!
! This module contains the definition of the psb_l_base_vect type which
! is a container for dense vectors.
!  This is encapsulated instead of being just a simple array to allow for
!  more complicated situations, such as GPU programming, where the memory
!  area we are interested in is not easily accessible from the host/Fortran
!  side. It is also meant to be encapsulated in an outer type, to allow
!  runtime switching as per the STATE design pattern, similar to the
!  sparse matrix types.
!
!
submodule (psb_l_base_vect_mod) psb_l_base_vect_impl

    use psi_serial_mod
    use psb_realloc_mod
    use psb_string_mod
contains

  !
  ! Constructors.
  !

  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  module function constructor(x) result(this)
    integer(psb_lpk_)   :: x(:)
    type(psb_l_base_vect_type) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x,kind=psb_ipk_),info)
  end function constructor
  


  !> Function  constructor:
  !! \brief     Constructor from size
  !!  \param    n   Size of vector to be built.
  !!
  module function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_l_base_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(n,info)

  end function size_const
  

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_l_base_vect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  module subroutine l_base_bld_x(x,this,scratch)
    implicit none
    integer(psb_lpk_), intent(in) :: this(:)
    class(psb_l_base_vect_type), intent(inout) :: x
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: i
    
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    call psb_realloc(size(this),x%v,info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_vect_bld')
      return
    end if
#if defined (PSB_OPENMP)
    !$omp parallel do private(i)
    do i = 1, size(this)
      x%v(i) = this(i)
    end do
#else
    x%v(:)  = this(:)
#endif
  end subroutine l_base_bld_x
  

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_mn:
  !! \memberof  psb_l_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  module subroutine l_base_bld_mn(x,n,scratch)
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_l_base_vect_type), intent(inout) :: x
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    call psb_realloc(n,x%v,info)
    call x%asb(n,info,scratch=scratch_)

  end subroutine l_base_bld_mn
  

  !> Function  bld_en:
  !! \memberof  psb_l_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  module subroutine l_base_bld_en(x,n,scratch)
    implicit none
    integer(psb_epk_), intent(in) :: n
    class(psb_l_base_vect_type), intent(inout) :: x
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: info

    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    call psb_realloc(n,x%v,info)
    call x%asb(n,info,scratch=scratch_)

  end subroutine l_base_bld_en
  

  !> Function  base_all:
  !! \memberof  psb_l_base_vect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  module subroutine l_base_all(n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)               :: n
    class(psb_l_base_vect_type), intent(out)    :: x
    integer(psb_ipk_), intent(out)              :: info

    call psb_realloc(n,x%v,info)
    if (try_newins) then
      call psb_realloc(n,x%iv,info)
      call x%set_ncfs(0)
    end if

  end subroutine l_base_all
  

  !> Function  base_mold:
  !! \memberof  psb_l_base_vect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  module subroutine l_base_mold(x, y, info)
    implicit none
    class(psb_l_base_vect_type), intent(in)   :: x
    class(psb_l_base_vect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info

    allocate(psb_l_base_vect_type :: y, stat=info)

  end subroutine l_base_mold
  

  module subroutine l_base_reinit(x, info,clear)
    implicit none
    class(psb_l_base_vect_type), intent(inout)    :: x
    integer(psb_ipk_), intent(out)              :: info
    logical, intent(in), optional               :: clear
    logical :: clear_

    info =  0
    if (present(clear)) then
      clear_ = clear
    else
      clear_ = .true.
    end if

    if (allocated(x%v)) then 
      if (x%is_dev()) call x%sync()
      if (clear_) x%v(:) = lzero
      call x%set_host()
      call x%set_upd()
    end if

  end subroutine l_base_reinit
  

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_ins:
  !! \memberof  psb_l_base_vect_type
  !! \brief Insert coefficients.
  !!
  !!
  !!         Given  a list of N pairs
  !!           (IRL(i),VAL(i))
  !!         record a new coefficient in X such that
  !!            X(IRL(1:N)) = VAL(1:N).
  !!
  !!         - the update operation will perform either
  !!               X(IRL(1:n)) = VAL(1:N)
  !!           or
  !!               X(IRL(1:n)) = X(IRL(1:n))+VAL(1:N)
  !!           according to the value of DUPLICATE.
  !!
  !!
  !!  \param n     number of pairs in input
  !!  \param irl(:)  the input row indices
  !!  \param val(:)  the input coefficients
  !!  \param dupl    how to treat duplicate entries
  !!  \param info  return code
  !!
  !
  module subroutine l_base_ins_a(n,irl,val,dupl,x,maxr,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl, maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz, dupl_, ncfs_, k

    info = 0
    if (psb_errstatus_fatal()) return

    if (try_newins) then 
      if (x%is_bld()) then
        ncfs_ = x%get_ncfs()
        isz = ncfs_ + n
        call psb_ensure_size(isz,x%v,info)
        call psb_ensure_size(isz,x%iv,info)
        k   = ncfs_
        do i = 1, n
          !loop over all val's rows
          ! row actual block row
          if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
            k  = k + 1 
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(k)  = val(i)
            x%iv(k) = irl(i)
          end if
        enddo
        call x%set_ncfs(k)      

      else if (x%is_upd()) then

        dupl_ = x%get_dupl()
        if (.not.allocated(x%v)) then
          info = psb_err_invalid_vect_state_
        else if (n > min(size(irl),size(val))) then
          info = psb_err_invalid_input_
        else
          isz = size(x%v)
          select case(dupl_)
          case(psb_dupl_ovwrt_)
            do i = 1, n
              !loop over all val's rows
              ! row actual block row
              if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
                ! this row belongs to me
                ! copy i-th row of block val in x
                x%v(irl(i)) = val(i)
              end if
            enddo

          case(psb_dupl_add_)

            do i = 1, n
              !loop over all val's rows
              if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
                ! this row belongs to me
                ! copy i-th row of block val in x
                x%v(irl(i)) = x%v(irl(i)) +  val(i)
              end if
            enddo

          case default
            info = 321
            ! !$      call psb_errpush(info,name)
            ! !$      goto 9999
          end select
        end if
      else
        info = psb_err_invalid_vect_state_
      end if
    else
      if (.not.allocated(x%v)) then
        info = psb_err_invalid_vect_state_
      else if (n > min(size(irl),size(val))) then
        info = psb_err_invalid_input_

      else
        isz = size(x%v)
        select case(dupl)
        case(psb_dupl_ovwrt_)
          do i = 1, n
            !loop over all val's rows
            ! row actual block row
            if ((1 <= irl(i)).and.(irl(i) <= isz)) then
              ! this row belongs to me
              ! copy i-th row of block val in x
              x%v(irl(i)) = val(i)
            end if
          enddo

        case(psb_dupl_add_)

          do i = 1, n
            !loop over all val's rows
            if ((1 <= irl(i)).and.(irl(i) <= isz)) then
              ! this row belongs to me
              ! copy i-th row of block val in x
              x%v(irl(i)) = x%v(irl(i)) +  val(i)
            end if
          enddo

        case default
          info = 321
          ! !$      call psb_errpush(info,name)
          ! !$      goto 9999
        end select
      end if
    end if
    call x%set_host()
    if (info /= 0) then
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine l_base_ins_a
  

  module subroutine l_base_ins_v(n,irl,val,dupl,x,maxr,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl, maxr
    class(psb_i_base_vect_type), intent(inout)  :: irl
    class(psb_l_base_vect_type), intent(inout)  :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: isz

    info = 0
    if (psb_errstatus_fatal()) return

    if (irl%is_dev()) call irl%sync()
    if (val%is_dev()) call val%sync()
    if (x%is_dev())   call x%sync()
    call x%ins(n,irl%v,val%v,dupl,maxr,info)

    if (info /= 0) then
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine l_base_ins_v
  


  !
  !> Function  base_zero
  !! \memberof  psb_l_base_vect_type
  !! \brief Zero out contents
  !!
  !
  module subroutine l_base_zero(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout)    :: x

    if (allocated(x%v)) then
      !$omp workshare
      x%v(:)=lzero
      !$omp end workshare
    end if
    call x%set_host()
  end subroutine l_base_zero
  


  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_l_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  module subroutine l_base_asb_m(n, x, info, scratch)
    implicit none
    integer(psb_mpk_), intent(in)              :: n
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: i, ncfs, xvsz
    integer(psb_lpk_), allocatable ::  vv(:)

    info = 0
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    if (try_newins) then 
      if (x%is_bld()) then
        ncfs = x%get_ncfs()
        xvsz = psb_size(x%v)
        call psb_realloc(n,vv,info)
        vv(:) = lzero
        select case(x%get_dupl())
        case(psb_dupl_add_)
          do i=1,ncfs
            vv(x%iv(i)) = vv(x%iv(i)) + x%v(i)
          end do
        case(psb_dupl_ovwrt_)
          do i=1,ncfs
            vv(x%iv(i)) = x%v(i)
          end do
        case(psb_dupl_err_)
          do i=1,ncfs
            if (vv(x%iv(i)).ne. lzero) then
              call psb_errpush(psb_err_duplicate_coo,'vect-asb')
              return
            else
              vv(x%iv(i)) = x%v(i)
            end if
          end do
        case default
          write(psb_err_unit,*) 'Error in vect_asb: unsafe dupl',x%get_dupl()
          info =-7
        end select
        call psb_move_alloc(vv,x%v,info)
        if (allocated(x%iv)) deallocate(x%iv,stat=info)
      else if (x%is_upd().or.x%is_asb().or.scratch_) then
        if (x%get_nrows() < n) &
             & call psb_realloc(n,x%v,info)
        if (info /= 0) &
             & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
      else
        info = psb_err_invalid_vect_state_        
        call psb_errpush(info,'vect_asb')
      end if
    else
      if (x%get_nrows() < n) &
           & call psb_realloc(n,x%v,info)
      if (info /= 0) &
           & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
    end if
    call x%set_host()
    call x%set_asb()
    call x%sync()
  end subroutine l_base_asb_m
  

  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_l_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  module subroutine l_base_asb_e(n, x, info, scratch)
    implicit none
    integer(psb_epk_), intent(in)              :: n
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: i, ncfs, xvsz
    integer(psb_lpk_), allocatable ::  vv(:)

    info = 0
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    if (try_newins) then
      if (info /= 0) &
           & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb unhandled')    
      if (x%is_bld()) then
        call psb_realloc(n,vv,info)
        vv(:) =  lzero
        select case(x%get_dupl())
        case(psb_dupl_add_)
          do i=1,x%get_ncfs()
            vv(x%iv(i)) = vv(x%iv(i)) + x%v(i)
          end do
        case(psb_dupl_ovwrt_)
          do i=1,x%get_ncfs()
            vv(x%iv(i)) = x%v(i)
          end do
        case(psb_dupl_err_)
          do i=1,x%get_ncfs()
            if (vv(x%iv(i)).ne. lzero) then
              call psb_errpush(psb_err_duplicate_coo,'vect_asb')
              return
            else
              vv(x%iv(i)) = x%v(i)
            end if
          end do
        case default
          write(psb_err_unit,*) 'Error in vect_asb: unsafe dupl',x%get_dupl()
          info =-7
        end select
        call psb_move_alloc(vv,x%v,info)
        if (allocated(x%iv)) deallocate(x%iv,stat=info)
      else if (x%is_upd().or.x%is_asb().or.scratch_) then
        if (x%get_nrows() < n) &
             & call psb_realloc(n,x%v,info)
        if (info /= 0) &
             & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
      else
        info = psb_err_invalid_vect_state_        
        call psb_errpush(info,'vect_asb')
      end if
    else
      if (x%get_nrows() < n) &
           & call psb_realloc(n,x%v,info)
      if (info /= 0) &
           & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
    end if
    call x%set_host()
    call x%set_asb()
    call x%sync()
  end subroutine l_base_asb_e
  

  !
  !> Function  base_free:
  !! \memberof  psb_l_base_vect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  module subroutine l_base_free(x, info)
    implicit none
    class(psb_l_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if ((info == 0).and.allocated(x%combuf)) call x%free_buffer(info)
    if ((info == 0).and.allocated(x%comid)) call x%free_comid(info)
    if ((info == 0).and.allocated(x%iv)) deallocate(x%iv, stat=info)    
    if (info /= 0) call &
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')
    call x%set_null()
  end subroutine l_base_free
  

  !
  !> Function  base_free_buffer:
  !! \memberof  psb_l_base_vect_type
  !! \brief Free aux buffer
  !!
  !!  \param info  return code
  !!
  !
  module subroutine l_base_free_buffer(x,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%combuf)) &
         &  deallocate(x%combuf,stat=info)
  end subroutine l_base_free_buffer
  

  !
  !> Function  base_maybe_free_buffer:
  !! \memberof  psb_l_base_vect_type
  !! \brief Conditionally Free aux buffer.
  !!        In some derived classes, e.g. GPU,
  !!        does not really frees to avoid  runtime
  !!        costs
  !!
  !!  \param info  return code
  !!
  !
  module subroutine l_base_maybe_free_buffer(x,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (psb_get_maybe_free_buffer())&
         &  call x%free_buffer(info)

  end subroutine l_base_maybe_free_buffer
  

  !
  !> Function  base_free_comid:
  !! \memberof  psb_l_base_vect_type
  !! \brief Free aux MPI communication id buffer
  !!
  !!  \param info  return code
  !!
  !
  module subroutine l_base_free_comid(x,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%comid)) &
         &  deallocate(x%comid,stat=info)
  end subroutine l_base_free_comid
  
  
  module function l_base_get_ncfs(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%ncfs
  end function l_base_get_ncfs
  

  module function l_base_get_dupl(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%dupl
  end function l_base_get_dupl
  

  module function l_base_get_state(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%bldstate 
  end function l_base_get_state
  

  module function l_base_is_null(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_null_)
  end function l_base_is_null
  

  module function l_base_is_bld(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_bld_)
  end function l_base_is_bld
  

  module function l_base_is_upd(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_upd_)
  end function l_base_is_upd
  

  module function l_base_is_asb(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_asb_)
  end function l_base_is_asb
  

  module subroutine  l_base_set_ncfs(n,x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%ncfs = n
  end subroutine l_base_set_ncfs
  

  module subroutine  l_base_set_dupl(n,x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%dupl = n
  end subroutine l_base_set_dupl
  

  module subroutine  l_base_set_state(n,x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%bldstate = n
  end subroutine l_base_set_state
  

  module subroutine  l_base_set_null(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

    x%bldstate = psb_vect_null_
  end subroutine l_base_set_null
  

  module subroutine  l_base_set_bld(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

    x%bldstate = psb_vect_bld_
  end subroutine l_base_set_bld
  

  module subroutine  l_base_set_upd(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

    x%bldstate = psb_vect_upd_
  end subroutine l_base_set_upd
  

  module subroutine  l_base_set_asb(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

    x%bldstate = psb_vect_asb_
  end subroutine l_base_set_asb
  

  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_sync:
  !! \memberof  psb_l_base_vect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  module subroutine l_base_sync(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

  end subroutine l_base_sync
  

  !
  !> Function  base_set_host:
  !! \memberof  psb_l_base_vect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  module subroutine l_base_set_host(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

  end subroutine l_base_set_host
  

  !
  !> Function  base_set_dev:
  !! \memberof  psb_l_base_vect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  module subroutine l_base_set_dev(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

  end subroutine l_base_set_dev
  

  !
  !> Function  base_set_sync:
  !! \memberof  psb_l_base_vect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  module subroutine l_base_set_sync(x)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x

  end subroutine l_base_set_sync
  

  !
  !> Function  base_is_dev:
  !! \memberof  psb_l_base_vect_type
  !! \brief Is  vector on external device    .
  !!
  !
  module function l_base_is_dev(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical  :: res

    res = .false.
  end function l_base_is_dev
  

  !
  !> Function  base_is_host
  !! \memberof  psb_l_base_vect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  module function l_base_is_host(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function l_base_is_host
  

  !
  !> Function  base_is_sync
  !! \memberof  psb_l_base_vect_type
  !! \brief Is  vector on sync               .
  !!
  !
  module function l_base_is_sync(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function l_base_is_sync
  

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  module subroutine l_base_cpy(x, y)
    implicit none
    class(psb_l_base_vect_type), intent(in)   :: x
    class(psb_l_base_vect_type), intent(out)  :: y

    if (allocated(x%v)) call y%bld(x%v)
    call y%set_state(x%get_state())
    call y%set_dupl(x%get_dupl())
    call y%set_ncfs(x%get_ncfs())
    if (allocated(x%iv)) y%iv = x%iv
  end subroutine l_base_cpy
  

  !
  ! Size info.
  !
  !
  !> Function  base_get_nrows
  !! \memberof  psb_l_base_vect_type
  !! \brief  Number of entries
  !!
  !
  module function l_base_get_nrows(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v)

  end function l_base_get_nrows
  

  !
  !> Function  base_get_sizeof
  !! \memberof  psb_l_base_vect_type
  !! \brief  Size in bytes
  !!
  !
  module function l_base_sizeof(x) result(res)
    implicit none
    class(psb_l_base_vect_type), intent(in) :: x
    integer(psb_epk_) :: res

    ! Force 8-byte integers.
    res = (1_psb_epk_ * psb_sizeof_lp) * x%get_nrows()

  end function l_base_sizeof
  

  !
  !> Function  base_get_fmt
  !! \memberof  psb_l_base_vect_type
  !! \brief  Format
  !!
  !
  module function l_base_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'BASE'
  end function l_base_get_fmt
  


  !
  !
  !
  !> Function  base_get_vect
  !! \memberof  psb_l_base_vect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  module function  l_base_get_vect(x,n) result(res)
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_lpk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n
    ! Local variables
    integer(psb_ipk_) :: isz, i

    if (.not.allocated(x%v)) return
    if (.not.x%is_host()) call x%sync()
    isz = x%get_nrows()
    if (present(n)) isz = max(0,min(isz,n))
    allocate(res(isz),stat=info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_get_vect')
      return
    end if
    if (.false.) then 
      res(1:isz) = x%v(1:isz)
    else
      !$omp parallel do private(i)
      do i=1, isz
        res(i) = x%v(i)
      end do
    end if
    
  end function l_base_get_vect
  

  !
  ! Reset all values
  !
  !
  !> Function  base_set_scal
  !! \memberof  psb_l_base_vect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  module subroutine l_base_set_scal(x,val,first,last)
    implicit none
    class(psb_l_base_vect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: first_, last_, i

    first_=1
    last_=size(x%v)
    if (present(first)) first_ = max(1,first)
    if (present(last))  last_  = min(last,last_)

    if (x%is_dev()) call x%sync()
#if defined(PSB_OPENMP)
    !$omp parallel do private(i)
    do i = first_, last_        
      x%v(i) = val
    end do
#else
    x%v(first_:last_) = val
#endif
    call x%set_host()

  end subroutine l_base_set_scal
  


  !
  !> Function  base_set_vect
  !! \memberof  psb_l_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  module subroutine l_base_set_vect(x,val,first,last)
    implicit none
    class(psb_l_base_vect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: first_, last_, i,  info

    if (.not.allocated(x%v)) then
      call psb_realloc(size(val),x%v,info)
    end if
    
    first_                     = 1
    if (present(first)) first_ = max(1,first)
    last_                      = min(psb_size(x%v),first_+size(val)-1)
    if (present(last))  last_  = min(last,last_)

    if (x%is_dev()) call x%sync()

#if defined(PSB_OPENMP)
      !$omp parallel do private(i)
      do i  = first_, last_
        x%v(i) = val(i-first_+1)
      end do
#else
      x%v(first_:last_) = val(1:last_-first_+1)
#endif
    call x%set_host()

  end subroutine l_base_set_vect
  

  module subroutine l_base_check_addr(x)
    class(psb_l_base_vect_type), intent(inout) :: x

    write(0,*) 'Check addr: base version, do nothing'

  end subroutine l_base_check_addr
  



  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_gthab
  !! \memberof  psb_l_base_vect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  module subroutine l_base_gthab(n,idx,alpha,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_base_vect_type) :: x

    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,alpha,x%v,beta,y)

  end subroutine l_base_gthab
  
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_l_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  module subroutine l_base_gthzv_x(i,n,idx,x,y)
    implicit none
    integer(psb_ipk_) :: i
    integer(psb_mpk_) :: n
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) ::  y(:)
    class(psb_l_base_vect_type) :: x

    if (idx%is_dev()) call idx%sync()
    call x%gth(n,idx%v(i:),y)

  end subroutine l_base_gthzv_x
  

  !
  ! New comm internals impl.
  !
  module subroutine l_base_gthzbuf(i,n,idx,x)
    implicit none
    integer(psb_ipk_) :: i
    integer(psb_mpk_) :: n
    class(psb_i_base_vect_type) :: idx
    class(psb_l_base_vect_type) :: x

    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if
    if (idx%is_dev()) call idx%sync()
    if (x%is_dev()) call x%sync()
    call x%gth(n,idx%v(i:),x%combuf(i:))

  end subroutine l_base_gthzbuf
  
  !
  !> Function  base_device_wait:
  !! \memberof  psb_l_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  module subroutine l_base_device_wait()
    implicit none

  end subroutine l_base_device_wait
  

  module function l_base_use_buffer() result(res)
    logical :: res

    res = .true.
  end function l_base_use_buffer
  

  module subroutine l_base_new_buffer(n,x,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,x%combuf,info)
  end subroutine l_base_new_buffer
  

  module subroutine l_base_new_comid(n,x,info)
    implicit none
    class(psb_l_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,2_psb_ipk_,x%comid,info)
  end subroutine l_base_new_comid
  


  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_l_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  module subroutine l_base_gthzv(n,idx,x,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_base_vect_type) :: x

    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,x%v,y)

  end subroutine l_base_gthzv
  

  !
  ! Scatter:
  ! Y(IDX(:)) = beta*Y(IDX(:)) + X(:)
  !
  !
  !> Function  base_sctb
  !! \memberof  psb_l_base_vect_type
  !! \brief scatter into a class(base_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  module subroutine l_base_sctb(n,idx,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_base_vect_type) :: y

    if (y%is_dev()) call y%sync()
    call psi_sct(n,idx,x,beta,y%v)
    call y%set_host()

  end subroutine l_base_sctb
  

  module subroutine l_base_sctb_x(i,n,idx,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_base_vect_type) :: y

    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),x,beta)
    call y%set_host()

  end subroutine l_base_sctb_x
  

  module subroutine l_base_sctb_buf(i,n,idx,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) :: beta
    class(psb_l_base_vect_type) :: y


    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    if (y%is_dev()) call y%sync()
    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),y%combuf(i:),beta)
    call y%set_host()

  end subroutine l_base_sctb_buf
  


end submodule psb_l_base_vect_impl


submodule (psb_l_base_multivect_mod) psb_l_base_multivect_impl
  use psb_realloc_mod
  use psi_serial_mod
  use psb_string_mod
contains

  !
  ! Constructors.
  !

  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  module function constructor(x) result(this)
    integer(psb_lpk_)   :: x(:,:)
    type(psb_l_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x,dim=1,kind=psb_ipk_),&
         & size(x,dim=2,kind=psb_ipk_),info)
  end function constructor
  


  !> Function  constructor:
  !! \brief     Constructor from size
  !!  \param    n   Size of vector to be built.
  !!
  module function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_l_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(m,n,info)

  end function size_const
  

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_l_base_multivect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  module subroutine l_base_mlv_bld_x(x,this)
    integer(psb_lpk_), intent(in) :: this(:,:)
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this,1),size(this,2),x%v,info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_mlv_vect_bld')
      return
    end if
    x%v(:,:)  = this(:,:)

  end subroutine l_base_mlv_bld_x
  

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_l_base_multivect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  module subroutine l_base_mlv_bld_n(x,m,n,scratch)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info
    logical, intent(in), optional        :: scratch

    call psb_realloc(m,n,x%v,info)
    call x%asb(m,n,info,scratch=scratch)

  end subroutine l_base_mlv_bld_n
  

  !> Function  base_mlv_all:
  !! \memberof  psb_l_base_multivect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  module subroutine l_base_mlv_all(m,n, x, info)
    implicit none
    integer(psb_ipk_), intent(in)               :: m,n
    class(psb_l_base_multivect_type), intent(out) :: x
    integer(psb_ipk_), intent(out)              :: info

    call psb_realloc(m,n,x%v,info)
    if (try_newins) then
      call psb_realloc(n,x%iv,info)
      call x%set_ncfs(0)
    end if

  end subroutine l_base_mlv_all
  

  !> Function  base_mlv_mold:
  !! \memberof  psb_l_base_multivect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  module subroutine l_base_mlv_mold(x, y, info)
    implicit none
    class(psb_l_base_multivect_type), intent(in)   :: x
    class(psb_l_base_multivect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info

    allocate(psb_l_base_multivect_type :: y, stat=info)

  end subroutine l_base_mlv_mold
  

  module subroutine l_base_mlv_reinit(x, info)
    implicit none
    class(psb_l_base_multivect_type), intent(out)    :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) then 
      call x%sync()
      x%v(:,:) = lzero
      call x%set_host()
      call x%set_upd()
    end if

  end subroutine l_base_mlv_reinit
  

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_mlv_ins:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Insert coefficients.
  !!
  !!
  !!         Given  a list of N pairs
  !!           (IRL(i),VAL(i))
  !!         record a new coefficient in X such that
  !!            X(IRL(1:N)) = VAL(1:N).
  !!
  !!         - the update operation will perform either
  !!               X(IRL(1:n)) = VAL(1:N)
  !!           or
  !!               X(IRL(1:n)) = X(IRL(1:n))+VAL(1:N)
  !!           according to the value of DUPLICATE.
  !!
  !!
  !!  \param n     number of pairs in input
  !!  \param irl(:)  the input row indices
  !!  \param val(:)  the input coefficients
  !!  \param dupl    how to treat duplicate entries
  !!  \param info  return code
  !!
  !
  module subroutine l_base_mlv_ins(n,irl,val,dupl,x,maxr,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl,maxr
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_lpk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz, nc, dupl_, ncfs_, k

    info = 0
    if (psb_errstatus_fatal()) return

    if (try_newins) then 
      if (x%is_bld()) then
        nc = size(x%v,2)
        ncfs_ = x%get_ncfs()
        isz = ncfs_ + n
        call psb_realloc(isz,nc,x%v,info)
        call psb_ensure_size(isz,x%iv,info)
        k   = ncfs_
        do i = 1, n
          !loop over all val's rows
          ! row actual block row
          if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
            k  = k + 1 
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(k,:)  = val(i,:)
            x%iv(k) = irl(i)
          end if
        enddo
        call x%set_ncfs(k)      

      else if (x%is_upd()) then

        dupl_ = x%get_dupl()
        if (.not.allocated(x%v)) then
          info = psb_err_invalid_vect_state_
        else if (n > min(size(irl),size(val))) then
          info = psb_err_invalid_input_
        else
          isz = size(x%v,1)
          nc = size(x%v,2)
          select case(dupl_)
          case(psb_dupl_ovwrt_)
            do i = 1, n
              !loop over all val's rows
              ! row actual block row
              if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
                ! this row belongs to me
                ! copy i-th row of block val in x
                x%v(irl(i),:) = val(i,:)
              end if
            enddo

          case(psb_dupl_add_)

            do i = 1, n
              !loop over all val's rows
              if ((1 <= irl(i)).and.(irl(i) <= maxr)) then
                ! this row belongs to me
                ! copy i-th row of block val in x
                x%v(irl(i),:) = x%v(irl(i),:) +  val(i,:)
              end if
            enddo

          case default
            info = 321
            ! !$      call psb_errpush(info,name)
            ! !$      goto 9999
          end select
        end if
      else
        info = psb_err_invalid_vect_state_
      end if
    else
      if (.not.allocated(x%v)) then
        info = psb_err_invalid_vect_state_
      else if (n > min(size(irl),size(val))) then
        info = psb_err_invalid_input_

      else
        isz = size(x%v,1)
        nc = size(x%v,2)
        select case(dupl)
        case(psb_dupl_ovwrt_)
          do i = 1, n
            !loop over all val's rows
            ! row actual block row
            if ((1 <= irl(i)).and.(irl(i) <= isz)) then
              ! this row belongs to me
              ! copy i-th row of block val in x
              x%v(irl(i),:) = val(i,:)
            end if
          enddo

        case(psb_dupl_add_)

          do i = 1, n
            !loop over all val's rows
            if ((1 <= irl(i)).and.(irl(i) <= isz)) then
              ! this row belongs to me
              ! copy i-th row of block val in x
              x%v(irl(i),:) = x%v(irl(i),:) +  val(i,:)
            end if
          enddo

        case default
          info = 321
          ! !$      call psb_errpush(info,name)
          ! !$      goto 9999
        end select
      end if
    end if
    call x%set_host()
    if (info /= 0) then
      call psb_errpush(info,'base_mlv_vect_ins')
      return
    end if

  end subroutine l_base_mlv_ins
  

  !
  !> Function  base_mlv_zero
  !! \memberof  psb_l_base_multivect_type
  !! \brief Zero out contents
  !!
  !
  module subroutine l_base_mlv_zero(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) x%v=lzero
    call x%set_host()

  end subroutine l_base_mlv_zero
  


  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_mlv_asb:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  module subroutine l_base_mlv_asb(m,n, x, info, scratch)
    implicit none
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    logical, intent(in), optional        :: scratch

    logical :: scratch_
    integer(psb_ipk_) :: i, ncfs, xvsz
    integer(psb_lpk_), allocatable ::  vv(:,:)

    info = 0
    if (present(scratch)) then
      scratch_ = scratch
    else
      scratch_ = .false.
    end if
    if (try_newins) then 
      if (x%is_bld()) then
        ncfs = x%get_ncfs()
        xvsz = psb_size(x%v)
        call psb_realloc(m,n,vv,info)
        vv(:,:) = lzero
        select case(x%get_dupl())
        case(psb_dupl_add_)
          do i=1,ncfs
            vv(x%iv(i),:) = vv(x%iv(i),:) + x%v(i,:)
          end do
        case(psb_dupl_ovwrt_)
          do i=1,ncfs
            vv(x%iv(i),:) = x%v(i,:)
          end do
        case(psb_dupl_err_)
          do i=1,ncfs
            if (any(vv(x%iv(i),:).ne.lzero)) then
              info = psb_err_duplicate_coo
              call psb_errpush(info,'mvect-asb')
              return
            else
              vv(x%iv(i),:) = x%v(i,:)
            end if
          end do
        case default
          write(psb_err_unit,*) 'Error in mvect_asb: unsafe dupl',x%get_dupl()
          info =-7
        end select
        call psb_move_alloc(vv,x%v,info)
        if (allocated(x%iv)) deallocate(x%iv,stat=info)
      else if (x%is_upd().or.x%is_asb().or.scratch_) then
        if ((x%get_nrows() < m).or.(x%get_ncols()<n)) &
             & call psb_realloc(m,n,x%v,info)
        if (info /= 0) then
          info = psb_err_alloc_dealloc_
          call psb_errpush(psb_err_alloc_dealloc_,'mvect_asb')
        end if
      else
        info = psb_err_invalid_vect_state_        
        call psb_errpush(info,'vect_asb')
      end if
    else
      if ((x%get_nrows() < m).or.(x%get_ncols()<n)) &
           & call psb_realloc(m,n,x%v,info)
      if (info /= 0) then
          info = psb_err_alloc_dealloc_
          call psb_errpush(psb_err_alloc_dealloc_,'mvect_asb')
        end if
    end if
    call x%set_host()
    call x%set_asb()
    call x%sync()
  end subroutine l_base_mlv_asb
  


  !
  !> Function  base_mlv_free:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  module subroutine l_base_mlv_free(x, info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info /= 0) call &
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')

  end subroutine l_base_mlv_free
  

  module function l_base_mlv_get_ncfs(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%ncfs
  end function l_base_mlv_get_ncfs
  

  module function l_base_mlv_get_dupl(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%dupl
  end function l_base_mlv_get_dupl
  

  module function l_base_mlv_get_state(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res
    res = x%bldstate 
  end function l_base_mlv_get_state
  

  module function l_base_mlv_is_null(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_null_)
  end function l_base_mlv_is_null
  

  module function l_base_mlv_is_bld(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_bld_)
  end function l_base_mlv_is_bld
  

  module function l_base_mlv_is_upd(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_upd_)
  end function l_base_mlv_is_upd
  

  module function l_base_mlv_is_asb(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical :: res
    res = (x%bldstate == psb_vect_asb_)
  end function l_base_mlv_is_asb
  

  module subroutine  l_base_mlv_set_ncfs(n,x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%ncfs = n
  end subroutine l_base_mlv_set_ncfs
  

  module subroutine  l_base_mlv_set_dupl(n,x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%dupl = n
  end subroutine l_base_mlv_set_dupl
  

  module subroutine  l_base_mlv_set_state(n,x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in) :: n
    x%bldstate = n
  end subroutine l_base_mlv_set_state
  

  module subroutine  l_base_mlv_set_null(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

    x%bldstate = psb_vect_null_
  end subroutine l_base_mlv_set_null
  

  module subroutine  l_base_mlv_set_bld(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

    x%bldstate = psb_vect_bld_
  end subroutine l_base_mlv_set_bld
  

  module subroutine  l_base_mlv_set_upd(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

    x%bldstate = psb_vect_upd_
  end subroutine l_base_mlv_set_upd
  

  module subroutine  l_base_mlv_set_asb(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

    x%bldstate = psb_vect_asb_
  end subroutine l_base_mlv_set_asb
  


  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_mlv_sync:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  module subroutine l_base_mlv_sync(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

  end subroutine l_base_mlv_sync
  

  !
  !> Function  base_mlv_set_host:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  module subroutine l_base_mlv_set_host(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

  end subroutine l_base_mlv_set_host
  

  !
  !> Function  base_mlv_set_dev:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  module subroutine l_base_mlv_set_dev(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

  end subroutine l_base_mlv_set_dev
  

  !
  !> Function  base_mlv_set_sync:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  module subroutine l_base_mlv_set_sync(x)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x

  end subroutine l_base_mlv_set_sync
  

  !
  !> Function  base_mlv_is_dev:
  !! \memberof  psb_l_base_multivect_type
  !! \brief Is  vector on external device    .
  !!
  !
  module function l_base_mlv_is_dev(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .false.
  end function l_base_mlv_is_dev
  

  !
  !> Function  base_mlv_is_host
  !! \memberof  psb_l_base_multivect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  module function l_base_mlv_is_host(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function l_base_mlv_is_host
  

  !
  !> Function  base_mlv_is_sync
  !! \memberof  psb_l_base_multivect_type
  !! \brief Is  vector on sync               .
  !!
  !
  module function l_base_mlv_is_sync(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function l_base_mlv_is_sync
  

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  module subroutine l_base_mlv_cpy(x, y)
    implicit none
    class(psb_l_base_multivect_type), intent(in)   :: x
    class(psb_l_base_multivect_type), intent(out)  :: y

    if (allocated(x%v)) call y%bld(x%v)
    call y%set_state(x%get_state())
    call y%set_dupl(x%get_dupl())
    call y%set_ncfs(x%get_ncfs())
    if (allocated(x%iv)) y%iv = x%iv
  end subroutine l_base_mlv_cpy
  


  !
  ! Size info.
  !
  !
  !> Function  base_mlv_get_nrows
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Number of entries
  !!
  !
  module function l_base_mlv_get_nrows(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,1)

  end function l_base_mlv_get_nrows
  

  module function l_base_mlv_get_ncols(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,2)

  end function l_base_mlv_get_ncols
  

  !
  !> Function  base_mlv_get_sizeof
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Size in bytesa
  !!
  !
  module function l_base_mlv_sizeof(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res

    ! Force 8-byte integers.
    res = (1_psb_epk_ * psb_sizeof_lp) * x%get_nrows() * x%get_ncols()

  end function l_base_mlv_sizeof
  

  !
  !> Function  base_mlv_get_fmt
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Format
  !!
  !
  module function l_base_mlv_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'BASE'
  end function l_base_mlv_get_fmt
  


  !
  !
  !
  !> Function  base_mlv_get_vect
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  module function  l_base_mlv_get_vect(x) result(res)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_lpk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info,m,n
    m = x%get_nrows()
    n = x%get_ncols()
    if (.not.allocated(x%v)) return
    call x%sync()
    allocate(res(m,n),stat=info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_mlv_get_vect')
      return
    end if
    res(1:m,1:n) = x%v(1:m,1:n)
  end function l_base_mlv_get_vect
  

  !
  ! Reset all values
  !
  !
  !> Function  base_mlv_set_scal
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  module subroutine l_base_mlv_set_scal(x,val)
    implicit none
    class(psb_l_base_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val

    integer(psb_ipk_) :: info
    x%v = val

  end subroutine l_base_mlv_set_scal
  

  !
  !> Function  base_mlv_set_vect
  !! \memberof  psb_l_base_multivect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  module subroutine l_base_mlv_set_vect(x,val)
    implicit none
    class(psb_l_base_multivect_type), intent(inout)  :: x
    integer(psb_lpk_), intent(in) :: val(:,:)
    integer(psb_ipk_) :: nr, nc
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      nr = min(size(x%v,1),size(val,1))
      nc = min(size(x%v,2),size(val,2))

      x%v(1:nr,1:nc) = val(1:nr,1:nc)
    else
      x%v = val
    end if

  end subroutine l_base_mlv_set_vect
  


  module function l_base_mlv_use_buffer() result(res)
    implicit none
    logical :: res

    res = .true.
  end function l_base_mlv_use_buffer

  module subroutine l_base_mlv_new_buffer(n,x,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    integer(psb_ipk_)               :: nc
    nc = x%get_ncols()
    call psb_realloc(n*nc,x%combuf,info)
  end subroutine l_base_mlv_new_buffer
  

  module subroutine l_base_mlv_new_comid(n,x,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,2_psb_ipk_,x%comid,info)
  end subroutine l_base_mlv_new_comid
  


  module subroutine l_base_mlv_maybe_free_buffer(x,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info


    info = 0
    if (psb_get_maybe_free_buffer())&
         &  call x%free_buffer(info)

  end subroutine l_base_mlv_maybe_free_buffer
  

  module subroutine l_base_mlv_free_buffer(x,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%combuf)) &
         &  deallocate(x%combuf,stat=info)
  end subroutine l_base_mlv_free_buffer
  

  module subroutine l_base_mlv_free_comid(x,info)
    implicit none
    class(psb_l_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%comid)) &
         &  deallocate(x%comid,stat=info)
  end subroutine l_base_mlv_free_comid
  


  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_mlv_gthab
  !! \memberof  psb_l_base_multivect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  module subroutine l_base_mlv_gthab(n,idx,alpha,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: alpha, beta, y(:)
    class(psb_l_base_multivect_type) :: x
    integer(psb_mpk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)
    call psi_gth(n,nc,idx,alpha,x%v,beta,y)

  end subroutine l_base_mlv_gthab
  
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_l_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  module subroutine l_base_mlv_gthzv_x(i,n,idx,x,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) ::  y(:)
    class(psb_l_base_multivect_type) :: x

    if (x%is_dev()) call x%sync()
    call x%gth(n,idx%v(i:),y)

  end subroutine l_base_mlv_gthzv_x
  

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_l_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  module subroutine l_base_mlv_gthzv(n,idx,x,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) ::  y(:)
    class(psb_l_base_multivect_type) :: x
    integer(psb_mpk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)

    call psi_gth(n,nc,idx,x%v,y)

  end subroutine l_base_mlv_gthzv
  
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_l_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  module subroutine l_base_mlv_gthzm(n,idx,x,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) ::  y(:,:)
    class(psb_l_base_multivect_type) :: x
    integer(psb_mpk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)

    call psi_gth(n,nc,idx,x%v,y)

  end subroutine l_base_mlv_gthzm
  

  !
  ! New comm internals impl.
  !
  module subroutine l_base_mlv_gthzbuf(i,ixb,n,idx,x)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i, ixb
    class(psb_i_base_vect_type) :: idx
    class(psb_l_base_multivect_type) :: x
    integer(psb_ipk_) :: nc

    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if
    if (idx%is_dev()) call idx%sync()
    if (x%is_dev()) call x%sync()
    nc = x%get_ncols()
    call x%gth(n,idx%v(i:),x%combuf(ixb:))

  end subroutine l_base_mlv_gthzbuf
  

  !
  ! Scatter:
  ! Y(IDX(:),:) = beta*Y(IDX(:),:) + X(:)
  !
  !
  !> Function  base_mlv_sctb
  !! \memberof  psb_l_base_multivect_type
  !! \brief scatter into a class(base_mlv_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  module subroutine l_base_mlv_sctb(n,idx,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: beta, x(:)
    class(psb_l_base_multivect_type) :: y
    integer(psb_mpk_) :: nc

    if (y%is_dev()) call y%sync()
    nc = psb_size(y%v,2_psb_ipk_)
    call psi_sct(n,nc,idx,x,beta,y%v)
    call y%set_host()

  end subroutine l_base_mlv_sctb
  

  module subroutine l_base_mlv_sctbr2(n,idx,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: idx(:)
    integer(psb_lpk_) :: beta, x(:,:)
    class(psb_l_base_multivect_type) :: y
    integer(psb_mpk_) :: nc

    if (y%is_dev()) call y%sync()
    nc = y%get_ncols()
    call psi_sct(n,nc,idx,x,beta,y%v)
    call y%set_host()

  end subroutine l_base_mlv_sctbr2
  

  module subroutine l_base_mlv_sctb_x(i,n,idx,x,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i
    class(psb_i_base_vect_type) :: idx
    integer( psb_lpk_) :: beta, x(:)
    class(psb_l_base_multivect_type) :: y

    call y%sct(n,idx%v(i:),x,beta)

  end subroutine l_base_mlv_sctb_x
  

  module subroutine l_base_mlv_sctb_buf(i,iyb,n,idx,beta,y)
    implicit none
    integer(psb_mpk_) :: n
    integer(psb_ipk_) :: i, iyb
    class(psb_i_base_vect_type) :: idx
    integer(psb_lpk_) :: beta
    class(psb_l_base_multivect_type) :: y
    integer(psb_ipk_) :: nc

    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    if (y%is_dev()) call y%sync()
    if (idx%is_dev()) call idx%sync()
    nc = y%get_ncols()
    call y%sct(n,idx%v(i:),y%combuf(iyb:),beta)
    call y%set_host()

  end subroutine l_base_mlv_sctb_buf
  

  !
  !> Function  base_device_wait:
  !! \memberof  psb_l_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  module subroutine l_base_mlv_device_wait()
    implicit none

  end subroutine l_base_mlv_device_wait
  

end submodule psb_l_base_multivect_impl
