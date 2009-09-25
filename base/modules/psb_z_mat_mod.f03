module psb_z_mat_mod

  use psb_z_base_mat_mod
  use psb_z_csr_mat_mod

  type :: psb_z_sparse_mat

    class(psb_z_base_sparse_mat), allocatable  :: a 

  contains
    ! Setters
    procedure, pass(a) :: set_nrows
    procedure, pass(a) :: set_ncols
    procedure, pass(a) :: set_dupl
    procedure, pass(a) :: set_state
    procedure, pass(a) :: set_null
    procedure, pass(a) :: set_bld
    procedure, pass(a) :: set_upd
    procedure, pass(a) :: set_asb
    procedure, pass(a) :: set_sorted
    procedure, pass(a) :: set_upper
    procedure, pass(a) :: set_lower
    procedure, pass(a) :: set_triangle
    procedure, pass(a) :: set_unit
    ! Getters
    procedure, pass(a) :: get_nrows
    procedure, pass(a) :: get_ncols
    procedure, pass(a) :: get_nzeros
    procedure, pass(a) :: get_nz_row
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
    procedure, pass(a) :: get_fmt => sparse_get_fmt
    procedure, pass(a) :: sizeof => z_sizeof


    ! Memory/data management 
    procedure, pass(a) :: csall
    procedure, pass(a) :: free
    procedure, pass(a) :: trim
    procedure, pass(a) :: csput 
    procedure, pass(a) :: z_csgetptn
    procedure, pass(a) :: z_csgetrow
    procedure, pass(a) :: z_csgetblk
    generic, public    :: csget => z_csgetptn, z_csgetrow, z_csgetblk 
    procedure, pass(a) :: csclip
    procedure, pass(a) :: reall => reallocate_nz
    procedure, pass(a) :: get_neigh
    procedure, pass(a) :: z_cscnv
    procedure, pass(a) :: z_cscnv_ip
    generic, public    :: cscnv => z_cscnv, z_cscnv_ip
    procedure, pass(a) :: reinit
    procedure, pass(a) :: print => sparse_print
    procedure, pass(a) :: z_mv_from
    generic, public    :: mv_from => z_mv_from
    procedure, pass(a) :: z_cp_from
    generic, public    :: cp_from => z_cp_from
    
    
    ! Computational routines 
    procedure, pass(a) :: get_diag
    procedure, pass(a) :: csnmi
    procedure, pass(a) :: z_csmv
    procedure, pass(a) :: z_csmm
    generic, public    :: csmm => z_csmm, z_csmv
    procedure, pass(a) :: z_scals
    procedure, pass(a) :: z_scal
    generic, public    :: scal => z_scals, z_scal 
    procedure, pass(a) :: z_cssv
    procedure, pass(a) :: z_cssm
    generic, public    :: cssm => z_cssm, z_cssv

  end type psb_z_sparse_mat

  private :: get_nrows, get_ncols, get_nzeros, get_size, &
       & get_state, get_dupl, is_null, is_bld, is_upd, &
       & is_asb, is_sorted, is_upper, is_lower, is_triangle, &
       & is_unit, get_neigh, csall, csput, z_csgetrow,&
       & z_csgetblk, csclip, z_cscnv, z_cscnv_ip, &
       & reallocate_nz, free, trim, &
       & sparse_print, reinit, &
       & set_nrows, set_ncols, set_dupl, &
       & set_state, set_null, set_bld, &
       & set_upd, set_asb, set_sorted, &
       & set_upper, set_lower, set_triangle, &
       & set_unit, get_diag, get_nz_row, z_csgetptn, &
       & z_mv_from, z_cp_from

  interface psb_sizeof
    module procedure z_sizeof
  end interface

  interface psb_move_alloc 
    module procedure z_sparse_mat_move
  end interface

  interface psb_clone
    module procedure z_sparse_mat_clone
  end interface

  interface psb_csmm
    module procedure z_csmm, z_csmv
  end interface

  interface psb_cssm
    module procedure z_cssm, z_cssv
  end interface

  interface psb_csnmi
    module procedure csnmi
  end interface
  
  interface psb_scal
    module procedure z_scals, z_scal
  end interface

contains 


  !=====================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  !=====================================

  
  function z_sizeof(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    
    res = 0
    if (allocated(a%a)) then 
      res = a%a%sizeof()
    end if
    
  end function z_sizeof



  function sparse_get_fmt(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    character(len=5) :: res

    if (allocated(a%a)) then 
      res = a%a%get_fmt()
    else
      res = 'NULL'
    end if

  end function sparse_get_fmt



  function get_dupl(a) result(res)
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psb_invalid_
    end if
  end function get_dupl


  function get_state(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psb_spmat_null_
    end if
  end function get_state

  function get_nrows(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function get_nrows

  function get_ncols(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function get_ncols

  function is_triangle(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function is_triangle

  function is_unit(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function is_unit

  function is_upper(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function is_upper

  function is_lower(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function is_lower

  function is_null(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function is_null

  function is_bld(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function is_bld

  function is_upd(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function is_upd

  function is_asb(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function is_asb

  function is_sorted(a) result(res)
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function is_sorted



  function get_nzeros(a) result(res)
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    Integer :: err_act, info
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
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res

    Integer :: err_act, info
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


  function get_nz_row(idx,a) result(res)
    use psb_error_mod
    implicit none 
    integer, intent(in)               :: idx
    class(psb_z_sparse_mat), intent(in) :: a
    integer :: res
    
    Integer :: err_act

    res = 0
    
    if (allocated(a%a)) res = a%a%get_nz_row(idx)

  end function get_nz_row



  !=====================================
  !
  !
  !
  ! Setters 
  !
  !
  !
  !
  !
  !
  !=====================================


  subroutine  set_nrows(m,a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    integer, intent(in) :: m
    Integer :: err_act, info
    character(len=20)  :: name='set_nrows'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_nrows(m)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_nrows

  subroutine  set_ncols(n,a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    call a%a%set_ncols(n)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_ncols


  subroutine  set_state(n,a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif
    call a%a%set_state(n)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_state


  subroutine  set_dupl(n,a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_dupl(n)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_dupl

  subroutine  set_null(a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_null()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_null

  subroutine  set_bld(a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_bld()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_bld

  subroutine  set_upd(a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_upd()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if


  end subroutine set_upd

  subroutine  set_asb(a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_asb()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_asb

  subroutine set_sorted(a,val) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_sorted(val)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_sorted

  subroutine set_triangle(a,val) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_triangle(val)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_triangle

  subroutine set_unit(a,val) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_unit(val)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_unit

  subroutine set_lower(a,val) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_lower(val)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_lower

  subroutine set_upper(a,val) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    Integer :: err_act, info
    character(len=20)  :: name='get_nzeros'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%set_upper(val)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine set_upper


  !=====================================
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
  !=====================================  


  subroutine sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_error_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psb_z_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act, info
    character(len=20)  :: name='sparse_print'
    logical, parameter :: debug=.false.

    info = 0
    call psb_get_erraction(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%print(iout,iv,eirs,eics,head,ivr,ivc)

    return

9999 continue

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine sparse_print



  subroutine get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a   
    integer, intent(in)                :: idx 
    integer, intent(out)               :: n   
    integer, allocatable, intent(out)  :: neigh(:)
    integer, intent(out)               :: info
    integer, optional, intent(in)      :: lev 

    Integer :: err_act
    character(len=20)  :: name='get_neigh'
    logical, parameter :: debug=.false.

    info = 0
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


  subroutine csall(nr,nc,a,info,nz) 
    use psb_z_base_mat_mod
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(out) :: a
    integer, intent(in)             :: nr,nc
    integer, intent(out)            :: info
    integer, intent(in), optional   :: nz

    Integer :: err_act 
    character(len=20)  :: name='csall'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)

    info = 0
    allocate(psb_z_coo_sparse_mat :: a%a, stat=info)
    if (info /= 0) then 
      info = 4000 
      call psb_errpush(info, name)
      goto 9999
    end if
    call a%a%allocate(nr,nc,nz)
    call a%set_bld() 

    return

9999 continue

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine csall

  subroutine  reallocate_nz(nz,a) 
    use psb_error_mod
    implicit none 
    integer, intent(in) :: nz
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%reallocate(nz)

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
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='free'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%free()
    deallocate(a%a) 
    return

9999 continue

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine free

  subroutine  trim(a) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='trim'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%trim()

    return

9999 continue

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine trim


  subroutine csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_z_base_mat_mod
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    complex(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)

    Integer :: err_act
    character(len=20)  :: name='csput'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (.not.a%is_bld()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info,gtl) 
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine csput

  subroutine z_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    use psb_z_base_mat_mod
    implicit none
    
    class(psb_z_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: imin,imax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale

    Integer :: err_act
    character(len=20)  :: name='csget'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call a%a%csget(imin,imax,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_csgetptn

  subroutine z_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    use psb_z_base_mat_mod
    implicit none
    
    class(psb_z_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: imin,imax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale

    Integer :: err_act
    character(len=20)  :: name='csget'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call a%a%csget(imin,imax,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_csgetrow



  subroutine z_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    use psb_z_base_mat_mod
    implicit none
    
    class(psb_z_sparse_mat), intent(in) :: a
    class(psb_z_sparse_mat), intent(out) :: b
    integer, intent(in)                  :: imin,imax
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax
    logical, intent(in), optional        :: rscale,cscale

    Integer :: err_act
    character(len=20)  :: name='csget'
    logical, parameter :: debug=.false.
    type(psb_z_coo_sparse_mat), allocatable  :: acoo


    info = 0
    call psb_erractionsave(err_act)
    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    allocate(acoo,stat=info)    
    
    if (info == 0) call a%a%csget(imin,imax,acoo,info,&
       & jmin,jmax,iren,append,rscale,cscale)
    if (info == 0) call move_alloc(acoo,b%a)
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_csgetblk



  subroutine csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    use psb_z_base_mat_mod
    implicit none
    
    class(psb_z_sparse_mat), intent(in) :: a
    class(psb_z_sparse_mat), intent(out) :: b
    integer,intent(out)                  :: info
    integer, intent(in), optional        :: imin,imax,jmin,jmax
    logical, intent(in), optional        :: rscale,cscale

    Integer :: err_act
    character(len=20)  :: name='csclip'
    logical, parameter :: debug=.false.
    type(psb_z_coo_sparse_mat), allocatable  :: acoo

    info = 0
    call psb_erractionsave(err_act)
    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    allocate(acoo,stat=info)    
    if (info == 0) call a%a%csclip(acoo,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
    if (info == 0) call move_alloc(acoo,b%a)
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine csclip



  subroutine z_cscnv(a,b,info,type,mold,upd,dupl)
    use psb_error_mod
    use psb_string_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in)    :: a
    class(psb_z_sparse_mat), intent(out)   :: b
    integer, intent(out)                   :: info
    integer,optional, intent(in)           :: dupl, upd
    character(len=*), optional, intent(in) :: type
    class(psb_z_base_sparse_mat), intent(in), optional :: mold


    class(psb_z_base_sparse_mat), allocatable  :: altmp
    Integer :: err_act
    character(len=20)  :: name='cscnv'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)

    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    if (present(dupl)) then 
      call b%set_dupl(dupl)
    else if (a%is_bld()) then 
      ! Does this make sense at all?? Who knows..
      call b%set_dupl(psb_dupl_def_)
    end if

    if (count( (/present(mold),present(type) /)) > 1) then
      info = 583
      call psb_errpush(info,name,a_err='TYPE, MOLD')
      goto 9999
    end if

    if (present(mold)) then 

      allocate(altmp, source=mold,stat=info) 

    else if (present(type)) then 

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_z_csr_sparse_mat :: altmp, stat=info) 
      case ('COO')
        allocate(psb_z_coo_sparse_mat :: altmp, stat=info) 
      case default
        info = 136 
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    else
      allocate(psb_z_csr_sparse_mat :: altmp, stat=info) 
    end if

    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if

    if (debug) write(0,*) 'Converting from ',&
         & a%get_fmt(),' to ',altmp%get_fmt()

    call altmp%cp_from_fmt(a%a, info)

    if (info /= 0) then
      info = 4010
      call psb_errpush(info,name,a_err="mv_from")
      goto 9999
    end if

    call move_alloc(altmp,b%a)
    call b%set_asb() 
    call b%trim()
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_cscnv


  subroutine z_cscnv_ip(a,info,type,mold,dupl)
    use psb_error_mod
    use psb_string_mod
    implicit none 

    class(psb_z_sparse_mat), intent(inout) :: a
    integer, intent(out)                   :: info
    integer,optional, intent(in)           :: dupl
    character(len=*), optional, intent(in) :: type
    class(psb_z_base_sparse_mat), intent(in), optional :: mold


    class(psb_z_base_sparse_mat), allocatable  :: altmp
    Integer :: err_act
    character(len=20)  :: name='cscnv_ip'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)

    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    if (present(dupl)) then 
      call a%set_dupl(dupl)
    else if (a%is_bld()) then 
      call a%set_dupl(psb_dupl_def_)
    end if

    if (count( (/present(mold),present(type) /)) > 1) then
      info = 583
      call psb_errpush(info,name,a_err='TYPE, MOLD')
      goto 9999
    end if

    if (present(mold)) then 

      allocate(altmp, source=mold,stat=info) 

    else if (present(type)) then 

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_z_csr_sparse_mat :: altmp, stat=info) 
      case ('COO')
        allocate(psb_z_coo_sparse_mat :: altmp, stat=info) 
      case default
        info = 136 
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    else
      allocate(psb_z_csr_sparse_mat :: altmp, stat=info) 
    end if

    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if

    if (debug) write(0,*) 'Converting in-place from ',&
         & a%get_fmt(),' to ',altmp%get_fmt()

    call altmp%mv_from_fmt(a%a, info)

    if (info /= 0) then
      info = 4010
      call psb_errpush(info,name,a_err="mv_from")
      goto 9999
    end if

    call move_alloc(altmp,a%a)
    call a%set_asb() 
    call a%trim()
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_cscnv_ip

  subroutine z_mv_from(a,b)
    use psb_error_mod
    use psb_string_mod
    implicit none 
    class(psb_z_sparse_mat), intent(out) :: a
    class(psb_z_base_sparse_mat), intent(inout) :: b
    integer :: info

    allocate(a%a,source=b, stat=info)
    call a%a%mv_from_fmt(b,info)
    
    return
  end subroutine z_mv_from

  subroutine z_cp_from(a,b)
    use psb_error_mod
    use psb_string_mod
    implicit none 
    class(psb_z_sparse_mat), intent(out) :: a
    class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
    Integer :: err_act, info
    character(len=20)  :: name='clone'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    
    allocate(a%a,source=b,stat=info)
    if (info /= 0) info = 4000
    if (info == 0) call a%a%cp_from_fmt(b, info)    
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine z_cp_from

  subroutine z_sparse_mat_move(a,b,info)
    use psb_error_mod
    use psb_string_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    class(psb_z_sparse_mat), intent(out)   :: b
    integer, intent(out)                   :: info

    Integer :: err_act
    character(len=20)  :: name='move_alloc'
    logical, parameter :: debug=.false.

    info = 0
    call move_alloc(a%a,b%a)
    
    return
  end subroutine z_sparse_mat_move

  subroutine z_sparse_mat_clone(a,b,info)
    use psb_error_mod
    use psb_string_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in)  :: a
    class(psb_z_sparse_mat), intent(out) :: b
    integer, intent(out)                 :: info

    Integer :: err_act
    character(len=20)  :: name='clone'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    
    allocate(b%a,source=a%a,stat=info)
    if (info /= 0) info = 4000
    if (info == 0) call b%a%cp_from_fmt(a%a, info)    
    if (info /= 0) goto 9999 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine z_sparse_mat_clone


  subroutine reinit(a,clear)
    use psb_error_mod
    implicit none 

    class(psb_z_sparse_mat), intent(inout) :: a   
    logical, intent(in), optional :: clear
    Integer :: err_act, info
    character(len=20)  :: name='reinit'

    call psb_erractionsave(err_act)
    if (a%is_null()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%reinit(clear)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if

  end subroutine reinit


  !=====================================
  !
  !
  !
  ! Computational routines
  !
  !
  !
  !
  !
  !
  !=====================================


  subroutine z_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
    complex(psb_dpk_), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psb_csmm'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%csmm(alpha,x,beta,y,info,trans) 
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

  end subroutine z_csmm

  subroutine z_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
    complex(psb_dpk_), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psb_csmv'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%csmm(alpha,x,beta,y,info,trans) 
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

  end subroutine z_csmv

  subroutine z_cssm(alpha,a,x,beta,y,info,trans,side,d) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
    complex(psb_dpk_), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans, side
    complex(psb_dpk_), intent(in), optional :: d(:)
    Integer :: err_act
    character(len=20)  :: name='psb_cssm'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%cssm(alpha,x,beta,y,info,trans,side,d) 
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

  end subroutine z_cssm

  subroutine z_cssv(alpha,a,x,beta,y,info,trans,side,d) 
    use psb_error_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
    complex(psb_dpk_), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans, side
    complex(psb_dpk_), intent(in), optional :: d(:)
    Integer :: err_act
    character(len=20)  :: name='psb_cssv'
    logical, parameter :: debug=.false.

    info = 0 
    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%cssm(alpha,x,beta,y,info,trans,side,d) 

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

  end subroutine z_cssv


  function csnmi(a) result(res)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    real(psb_dpk_)         :: res

    Integer :: err_act, info
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.

    call psb_get_erraction(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    res = a%a%csnmi()


    return

9999 continue

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end function csnmi



  subroutine get_diag(a,d,info)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_z_sparse_mat), intent(in) :: a
    complex(psb_dpk_), intent(out)      :: d(:)
    integer, intent(out)                :: info

    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%get_diag(d,info)
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

  end subroutine get_diag

  subroutine z_scal(d,a,info)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    complex(psb_dpk_), intent(in)          :: d(:)
    integer, intent(out)                   :: info

    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%scal(d,info)
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

  end subroutine z_scal

  subroutine z_scals(d,a,info)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_z_sparse_mat), intent(inout) :: a
    complex(psb_dpk_), intent(in)          :: d
    integer, intent(out)                   :: info

    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%scal(d,info)
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

  end subroutine z_scals


end module psb_z_mat_mod
