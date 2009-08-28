
module psbn_d_mat_mod

  use psbn_d_base_mat_mod
  
  type :: psbn_d_sparse_mat

    class(psbn_d_base_sparse_mat), allocatable  :: a 
    
  contains
    
    procedure, pass(a) :: get_nrows
    procedure, pass(a) :: get_ncols
    procedure, pass(a) :: get_nzeros
    procedure, pass(a) :: get_size
    procedure, pass(a) :: get_state

    procedure, pass(a) :: get_dupl
    procedure, pass(a) :: is_null
    procedure, pass(a) :: is_bld
    procedure, pass(a) :: is_upd
    procedure, pass(a) :: is_asb
    procedure, pass(a) :: is_sorted
    procedure, pass(a) :: is_upper
    procedure, pass(a) :: is_lower
    procedure, pass(a) :: is_triangle
    procedure, pass(a) :: is_unit
    procedure, pass(a) :: get_neigh
    procedure, pass(a) :: allocate_mn
    procedure, pass(a) :: allocate_mnnz
    procedure, pass(a) :: reallocate_nz
    procedure, pass(a) :: free
    generic,   public  :: allocate => allocate_mn, allocate_mnnz
    generic,   public  :: reallocate => reallocate_nz

    procedure, pass(a) :: d_csmv
    procedure, pass(a) :: d_csmm
    generic, public    :: psbn_csmm => d_csmm, d_csmv

    procedure, pass(a) :: d_cssv
    procedure, pass(a) :: d_cssm
    generic, public    :: psbn_cssm => d_cssm, d_cssv
    
  end type psbn_d_sparse_mat

  private :: get_nrows, get_ncols, get_nzeros, get_size, &
       & get_state, get_dupl, is_null, is_bld, is_upd, &
       & is_asb, is_sorted, is_upper, is_lower, is_triangle, &
       & is_unit, get_neigh, allocate_mn, allocate_mnnz, &
       & reallocate_nz, free, d_csmv, d_csmm, d_cssv, d_cssm 

contains 

  function get_dupl(a) result(res)
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psbn_invalid_
    end if
  end function get_dupl
 
 
  function get_state(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psbn_spmat_null_
    end if
  end function get_state
 
  function get_nrows(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function get_nrows

  function get_ncols(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function get_ncols

  function is_triangle(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function is_triangle

  function is_unit(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function is_unit

  function is_upper(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function is_upper

  function is_lower(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function is_lower

  function is_null(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res
    
    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function is_null

  function is_bld(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function is_bld

  function is_upd(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function is_upd

  function is_asb(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function is_asb

  function is_sorted(a) result(res)
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function is_sorted


  function get_nzeros(a) result(res)
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    res = a%a%get_nzeros()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end function get_nzeros

  function get_size(a) result(res)
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act
    character(len=20)  :: name='get_size'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    
    res = a%a%get_size()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end function get_size


  subroutine get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a   
    integer, intent(in)                :: idx 
    integer, intent(out)               :: n   
    integer, allocatable, intent(out)  :: neigh(:)
    integer, intent(out)               :: info
    integer, optional, intent(in)      :: lev 
    
    Integer :: err_act
    character(len=20)  :: name='get_neigh'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%get_neigh(idx,neigh,n,info,lev)

    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine get_neigh

  subroutine  allocate_mn(m,n,a,type,mold) 
    use psb_error_mod
    use psb_string_mod
    integer, intent(in) :: m,n
    class(psbn_d_sparse_mat), intent(inout) :: a
    character(len=*), intent(in), optional :: type
    class(psbn_d_base_sparse_mat), intent(in), optional :: mold

    Integer :: err_act, info
    character(len=20)  :: name='allocate_mn'
    character(len=8)   :: type_
    logical, parameter :: debug=.false.


    call psb_erractionsave(err_act)
    info = 0 
    if (allocated(a%a)) then 
      call a%a%free()
      deallocate(a%a)
    end if

    if (present(mold)) then 
      allocate(a%a, source=mold, stat=info)

    else

      if (present(type)) then 
        type_ = psb_toupper(type)
      else
        type_ = 'COO'
      end if

      select case(type) 
      case('COO')
        allocate(psbn_d_coo_sparse_mat :: a%a, stat=info)
! Add here a few other data structures inplemented by default.

!!$      case('CSR') 
!!$        allocate(psbn_d_csr_sparse_mat :: a%a, stat=info)

      case default
        allocate(psbn_d_coo_sparse_mat :: a%a, stat=info)
      end select

    end if

    if (info /= 0) then 
      info = 4010
      goto 9999
    end if
    
    call a%a%allocate(m,n)
          
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine allocate_mn

  subroutine  allocate_mnnz(m,n,nz,a,type,mold) 
    use psb_error_mod
    use psb_string_mod
    integer, intent(in) :: m,n,nz
    class(psbn_d_sparse_mat), intent(inout) :: a
    character(len=*), intent(in), optional :: type
    class(psbn_d_base_sparse_mat), intent(in), optional :: mold

    Integer :: err_act, info
    character(len=20)  :: name='allocate_mnnz'
    character(len=8)   :: type_
    logical, parameter :: debug=.false.


    call psb_erractionsave(err_act)
    info = 0 
    if (allocated(a%a)) then 
      call a%a%free()
      deallocate(a%a)
    end if

    if (present(mold)) then 
      allocate(a%a, source=mold, stat=info)

    else

      if (present(type)) then 
        type_ = psb_toupper(type)
      else
        type_ = 'COO'
      end if

      select case(type) 
      case('COO')
        allocate(psbn_d_coo_sparse_mat :: a%a, stat=info)
! Add here a few other data structures inplemented by default.

!!$      case('CSR') 
!!$        allocate(psbn_d_csr_sparse_mat :: a%a, stat=info)

      case default
        allocate(psbn_d_coo_sparse_mat :: a%a, stat=info)
      end select

    end if

    if (info /= 0) then 
      info = 4010
      goto 9999
    end if
    
    call a%a%allocate(m,n,nz)
          
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine allocate_mnnz

  subroutine  reallocate_nz(nz,a) 
    use psb_error_mod
    integer, intent(in) :: nz
    class(psbn_d_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='reallocate_nz'
    logical, parameter :: debug=.false.

    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%reallocate(nz)

    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine reallocate_nz

  subroutine  free(a) 
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='free'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%free()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine free


  subroutine d_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_csmm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%psbn_csmm(alpha,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_csmm

  subroutine d_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_csmv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%psbn_csmm(alpha,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_csmv

  subroutine d_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_cssm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    
    call a%a%psbn_cssm(alpha,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_cssm

  subroutine d_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_cssv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%psbn_cssm(alpha,x,beta,y,info,trans) 


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_cssv

end module psbn_d_mat_mod

