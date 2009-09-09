
module psbn_d_mat_mod

  use psbn_d_base_mat_mod
  use psbn_d_csr_mat_mod

  type :: psbn_d_sparse_mat

    class(psbn_d_base_sparse_mat), allocatable  :: a 
    
  contains
    
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
    procedure, pass(a) :: allocate_mnnz
    procedure, pass(a) :: reallocate_nz
    procedure, pass(a) :: free
    procedure, pass(a) :: print => sparse_print
    procedure, pass(a) :: get_fmt => sparse_get_fmt

    generic,   public  :: allocate => allocate_mnnz
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
       & is_unit, get_neigh, allocate_mnnz, &
       & reallocate_nz, free, d_csmv, d_csmm, d_cssv, d_cssm, sparse_print, &
       & set_nrows, set_ncols, set_dupl, set_state, set_null, set_bld, &
       & set_upd, set_asb, set_sorted, set_upper, set_lower, set_triangle, &
       & set_unit



  interface psbn_csall
    subroutine psbn_d_csall(nr,nc,a,info,nz) 
      use psbn_d_base_mat_mod
      import psbn_d_sparse_mat
      type(psbn_d_sparse_mat), intent(out) :: a
      integer, intent(in)             :: nr,nc
      integer, intent(out)            :: info
      integer, intent(in), optional   :: nz
    end subroutine psbn_d_csall
  end interface

  interface psbn_csins
    subroutine psbn_d_csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
      use psbn_d_base_mat_mod
      import psbn_d_sparse_mat
      type(psbn_d_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine psbn_d_csins
  end interface

  interface psbn_cscnv
    subroutine psbn_d_spcnv(a,b,info,type,mold,upd,dupl)
      use psbn_d_base_mat_mod
      import psbn_d_sparse_mat
      type(psbn_d_sparse_mat), intent(in)    :: a
      type(psbn_d_sparse_mat), intent(out)   :: b
      integer, intent(out)                   :: info
      integer, optional, intent(in)           :: dupl, upd
      character(len=*), optional, intent(in) :: type
      class(psbn_d_base_sparse_mat), intent(in), optional :: mold
      
    end subroutine psbn_d_spcnv
    subroutine psbn_d_spcnv_ip(a,info,type,mold,dupl)
      use psbn_d_base_mat_mod
      import psbn_d_sparse_mat
      type(psbn_d_sparse_mat), intent(inout)  :: a
      integer, intent(out)                    :: info
      integer,optional, intent(in)            :: dupl
      character(len=*), optional, intent(in)  :: type
      class(psbn_d_base_sparse_mat), intent(in), optional :: mold
    end subroutine psbn_d_spcnv_ip
  end interface

contains 

 
  function sparse_get_fmt(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
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
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_dupl()
    else
      res = psbn_invalid_
    end if
  end function get_dupl
 
 
  function get_state(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_state()
    else
      res = psbn_spmat_null_
    end if
  end function get_state
 
  function get_nrows(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res
    
    if (allocated(a%a)) then 
      res = a%a%get_nrows()
    else
      res = 0
    end if

  end function get_nrows

  function get_ncols(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    integer :: res

    if (allocated(a%a)) then 
      res = a%a%get_ncols()
    else
      res = 0
    end if

  end function get_ncols

  function is_triangle(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_triangle()
    else
      res = .false.
    end if

  end function is_triangle

  function is_unit(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_unit()
    else
      res = .false.
    end if

  end function is_unit

  function is_upper(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upper()
    else
      res = .false.
    end if

  end function is_upper

  function is_lower(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = .not. a%a%is_upper()
    else
      res = .false.
    end if

  end function is_lower

  function is_null(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res
    
    if (allocated(a%a)) then 
      res = a%a%is_null() 
    else
      res = .true.
    end if

  end function is_null

  function is_bld(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_bld()
    else
      res = .false.
    end if

  end function is_bld

  function is_upd(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_upd()
    else
      res = .false.
    end if

  end function is_upd

  function is_asb(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_asb()
    else
      res = .false.
    end if

  end function is_asb

  function is_sorted(a) result(res)
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    logical :: res

    if (allocated(a%a)) then 
      res = a%a%is_sorted()
    else
      res = .false.
    end if

  end function is_sorted

 
  subroutine  set_nrows(m,a) 
    use psb_error_mod
    implicit none 
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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
    class(psbn_d_sparse_mat), intent(inout) :: a
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




  function get_nzeros(a) result(res)
    use psb_error_mod
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
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
    class(psbn_d_sparse_mat), intent(in) :: a
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

  subroutine sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_error_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psbn_d_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act, info
    character(len=20)  :: name='sparse_print'
    logical, parameter :: debug=.false.

    info = 0
    call psb_erractionsave(err_act)
    if (.not.allocated(a%a)) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    call a%a%print(iout,iv,eirs,eics,head,ivr,ivc)


    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine sparse_print

  subroutine get_neigh(a,idx,neigh,n,info,lev)
    use psb_error_mod
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a   
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


  subroutine  allocate_mnnz(m,n,a,nz,type,mold) 
    use psb_error_mod
    use psb_string_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psbn_d_sparse_mat), intent(inout) :: a
    integer, intent(in), optional :: nz
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
    implicit none 
    integer, intent(in) :: nz
    class(psbn_d_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
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
    implicit none 
    class(psbn_d_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
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
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_csmm'
    logical, parameter :: debug=.false.

    info = 0
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
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_csmv'
    logical, parameter :: debug=.false.

    info = 0
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
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:,:)
    real(kind(1.d0)), intent(inout) :: y(:,:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_cssm'
    logical, parameter :: debug=.false.

    info = 0
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
    implicit none 
    class(psbn_d_sparse_mat), intent(in) :: a
    real(kind(1.d0)), intent(in)    :: alpha, beta, x(:)
    real(kind(1.d0)), intent(inout) :: y(:)
    integer, intent(out)            :: info
    character, optional, intent(in) :: trans
    Integer :: err_act
    character(len=20)  :: name='psbn_cssv'
    logical, parameter :: debug=.false.

    info = 0 
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

