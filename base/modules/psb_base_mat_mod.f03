module psb_base_mat_mod
  
  use psb_const_mod 
  use psi_serial_mod

  type  :: psb_base_sparse_mat
    integer, private     :: m, n
    integer, private     :: state, duplicate 
    logical, private     :: triangle, unitd, upper, sorted
    ! This is a different animal: it's a kitchen sink for
    ! any additional parameters that may be needed
    ! when converting to/from COO. Why here?
    ! Will tell you one day...
    integer, allocatable :: aux(:)
  contains 

    ! == = =================================
    !
    ! Getters 
    !
    !
    ! == = =================================
    procedure, pass(a) :: get_nrows => psb_base_get_nrows
    procedure, pass(a) :: get_ncols => psb_base_get_ncols
    procedure, pass(a) :: get_nzeros => psb_base_get_nzeros
    procedure, pass(a) :: get_nz_row => psb_base_get_nz_row
    procedure, pass(a) :: get_size => psb_base_get_size
    procedure, pass(a) :: get_state => psb_base_get_state
    procedure, pass(a) :: get_dupl => psb_base_get_dupl
    procedure, pass(a) :: get_fmt => psb_base_get_fmt
    procedure, pass(a) :: get_aux => psb_base_get_aux
    procedure, pass(a) :: is_null => psb_base_is_null
    procedure, pass(a) :: is_bld => psb_base_is_bld
    procedure, pass(a) :: is_upd => psb_base_is_upd
    procedure, pass(a) :: is_asb => psb_base_is_asb
    procedure, pass(a) :: is_sorted => psb_base_is_sorted
    procedure, pass(a) :: is_upper => psb_base_is_upper
    procedure, pass(a) :: is_lower => psb_base_is_lower
    procedure, pass(a) :: is_triangle => psb_base_is_triangle
    procedure, pass(a) :: is_unit => psb_base_is_unit
    
    ! == = =================================
    !
    ! Setters 
    !
    ! == = =================================
    procedure, pass(a) :: set_nrows => psb_base_set_nrows
    procedure, pass(a) :: set_ncols => psb_base_set_ncols
    procedure, pass(a) :: set_dupl => psb_base_set_dupl
    procedure, pass(a) :: set_state => psb_base_set_state
    procedure, pass(a) :: set_null => psb_base_set_null
    procedure, pass(a) :: set_bld => psb_base_set_bld
    procedure, pass(a) :: set_upd => psb_base_set_upd
    procedure, pass(a) :: set_asb => psb_base_set_asb
    procedure, pass(a) :: set_sorted => psb_base_set_sorted
    procedure, pass(a) :: set_upper => psb_base_set_upper
    procedure, pass(a) :: set_lower => psb_base_set_lower
    procedure, pass(a) :: set_triangle => psb_base_set_triangle
    procedure, pass(a) :: set_unit => psb_base_set_unit
    procedure, pass(a) :: set_aux => psb_base_set_aux


    ! == = =================================
    !
    ! Data management
    !
    ! == = =================================  
    procedure, pass(a) :: get_neigh => psb_base_get_neigh
    procedure, pass(a) :: free => psb_base_free
    procedure, pass(a) :: trim => psb_base_trim
    procedure, pass(a) :: reinit => psb_base_reinit
    procedure, pass(a) :: allocate_mnnz => psb_base_allocate_mnnz
    procedure, pass(a) :: reallocate_nz => psb_base_reallocate_nz
    generic,   public  :: allocate => allocate_mnnz
    generic,   public  :: reallocate => reallocate_nz
    procedure, pass(a) :: csgetptn => psb_base_csgetptn
    generic, public    :: csget => csgetptn
    procedure, pass(a) :: print => psb_base_sparse_print
    procedure, pass(a) :: sizeof => psb_base_sizeof
    procedure, pass(a) :: psb_base_cp_from
    generic, public    :: cp_from => psb_base_cp_from
    procedure, pass(a) :: psb_base_mv_from
    generic, public    :: mv_from => psb_base_mv_from
    procedure, pass(a) :: transp_1mat => psb_base_transp_1mat
    procedure, pass(a) :: transp_2mat => psb_base_transp_2mat
    generic, public    :: transp => transp_1mat, transp_2mat
    procedure, pass(a) :: transc_1mat => psb_base_transc_1mat
    procedure, pass(a) :: transc_2mat => psb_base_transc_2mat
    generic, public    :: transc => transc_1mat, transc_2mat
 
  end type psb_base_sparse_mat


  interface 
    function psb_base_get_nz_row(idx,a) result(res)
      import psb_base_sparse_mat, psb_long_int_k_
      integer, intent(in)                    :: idx
      class(psb_base_sparse_mat), intent(in) :: a
      integer :: res
    end function psb_base_get_nz_row
  end interface
  
  interface 
    function psb_base_get_nzeros(a) result(res)
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(in) :: a
      integer :: res
    end function psb_base_get_nzeros
  end interface

  interface 
    function psb_base_get_size(a) result(res)
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(in) :: a
      integer :: res
    end function psb_base_get_size
  end interface

  interface 
    subroutine psb_base_reinit(a,clear)
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(inout) :: a   
      logical, intent(in), optional :: clear
    end subroutine psb_base_reinit
  end interface
  
  interface 
    subroutine psb_base_sparse_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      import psb_base_sparse_mat, psb_long_int_k_
      integer, intent(in)               :: iout
      class(psb_base_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
    end subroutine psb_base_sparse_print
  end interface

  interface 
    subroutine psb_base_csgetptn(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine psb_base_csgetptn
  end interface
  
  interface 
    subroutine psb_base_get_neigh(a,idx,neigh,n,info,lev)
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(in) :: a   
      integer, intent(in)                :: idx 
      integer, intent(out)               :: n   
      integer, allocatable, intent(out)  :: neigh(:)
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: lev 
    end subroutine psb_base_get_neigh
  end interface
  
  interface 
    subroutine  psb_base_allocate_mnnz(m,n,a,nz) 
      import psb_base_sparse_mat, psb_long_int_k_
      integer, intent(in) :: m,n
      class(psb_base_sparse_mat), intent(inout) :: a
      integer, intent(in), optional  :: nz
    end subroutine psb_base_allocate_mnnz
  end interface

  interface 
    subroutine psb_base_reallocate_nz(nz,a) 
      import psb_base_sparse_mat, psb_long_int_k_
      integer, intent(in) :: nz
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_reallocate_nz
  end interface

  interface 
    subroutine psb_base_free(a) 
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_free
  end interface
  
  interface 
    subroutine psb_base_trim(a) 
      import psb_base_sparse_mat, psb_long_int_k_
      class(psb_base_sparse_mat), intent(inout) :: a
    end subroutine psb_base_trim
  end interface
  
 
contains

  function psb_base_sizeof(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8
  end function psb_base_sizeof
 
  function psb_base_get_fmt(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'NULL'
  end function psb_base_get_fmt
  
  function psb_base_get_dupl(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%duplicate
  end function psb_base_get_dupl
 
 
  function psb_base_get_state(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%state
  end function psb_base_get_state
 
  function psb_base_get_nrows(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%m
  end function psb_base_get_nrows

  function psb_base_get_ncols(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer :: res
    res = a%n
  end function psb_base_get_ncols

 
  subroutine  psb_base_set_aux(v,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: v(:)
    ! TBD
    write(0,*) 'SET_AUX is empty right now '
  end subroutine psb_base_set_aux

  subroutine  psb_base_get_aux(v,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    integer, intent(out), allocatable  :: v(:)
    ! TBD
    write(0,*) 'GET_AUX is empty right now '
  end subroutine psb_base_get_aux
 
  subroutine  psb_base_set_nrows(m,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: m
    a%m = m
  end subroutine psb_base_set_nrows

  subroutine  psb_base_set_ncols(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%n = n
  end subroutine psb_base_set_ncols
  

  subroutine  psb_base_set_state(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%state = n
  end subroutine psb_base_set_state


  subroutine  psb_base_set_dupl(n,a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    integer, intent(in) :: n
    a%duplicate = n
  end subroutine psb_base_set_dupl

  subroutine  psb_base_set_null(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_null_
  end subroutine psb_base_set_null

  subroutine  psb_base_set_bld(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_bld_
  end subroutine psb_base_set_bld

  subroutine  psb_base_set_upd(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_upd_
  end subroutine psb_base_set_upd

  subroutine  psb_base_set_asb(a) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a

    a%state = psb_spmat_asb_
  end subroutine psb_base_set_asb

  subroutine psb_base_set_sorted(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%sorted = val
    else
      a%sorted = .true.
    end if
  end subroutine psb_base_set_sorted

  subroutine psb_base_set_triangle(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%triangle = val
    else
      a%triangle = .true.
    end if
  end subroutine psb_base_set_triangle

  subroutine psb_base_set_unit(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%unitd = val
    else
      a%unitd = .true.
    end if
  end subroutine psb_base_set_unit

  subroutine psb_base_set_lower(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = .not.val
    else
      a%upper = .false.
    end if
  end subroutine psb_base_set_lower

  subroutine psb_base_set_upper(a,val) 
    implicit none 
    class(psb_base_sparse_mat), intent(inout) :: a
    logical, intent(in), optional :: val
    
    if (present(val)) then 
      a%upper = val
    else
      a%upper = .true.
    end if
  end subroutine psb_base_set_upper

  function psb_base_is_triangle(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%triangle
  end function psb_base_is_triangle

  function psb_base_is_unit(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%unitd
  end function psb_base_is_unit

  function psb_base_is_upper(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%upper
  end function psb_base_is_upper

  function psb_base_is_lower(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = .not.a%upper
  end function psb_base_is_lower

  function psb_base_is_null(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_null_)
  end function psb_base_is_null

  function psb_base_is_bld(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_bld_)
  end function psb_base_is_bld

  function psb_base_is_upd(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_upd_)
  end function psb_base_is_upd

  function psb_base_is_asb(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = (a%state == psb_spmat_asb_)
  end function psb_base_is_asb

  function psb_base_is_sorted(a) result(res)
    implicit none 
    class(psb_base_sparse_mat), intent(in) :: a
    logical :: res
    res = a%sorted
  end function psb_base_is_sorted

  subroutine psb_base_mv_from(a,b)
    implicit none 

    class(psb_base_sparse_mat), intent(out)   :: a
    type(psb_base_sparse_mat), intent(inout) :: b

    a%m         = b%m
    a%n         = b%n
    a%state     = b%state
    a%duplicate = b%duplicate
    a%triangle  = b%triangle
    a%unitd     = b%unitd
    a%upper     = b%upper
    a%sorted    = b%sorted
    call move_alloc(b%aux,a%aux)

  end subroutine psb_base_mv_from
  
  subroutine psb_base_cp_from(a,b)
    implicit none 

    class(psb_base_sparse_mat), intent(out) :: a
    type(psb_base_sparse_mat), intent(in)  :: b

    a%m         = b%m
    a%n         = b%n
    a%state     = b%state
    a%duplicate = b%duplicate
    a%triangle  = b%triangle
    a%unitd     = b%unitd
    a%upper     = b%upper
    a%sorted    = b%sorted
    if (allocated(b%aux)) then 
      allocate(a%aux(size(b%aux)))
      a%aux(:)       = b%aux(:)
    end if

  end subroutine psb_base_cp_from

  subroutine psb_base_transp_2mat(a,b)
    implicit none 
    
    class(psb_base_sparse_mat), intent(out) :: a
    class(psb_base_sparse_mat), intent(in)  :: b
    
    a%m         = b%n
    a%n         = b%m
    a%state     = b%state
    a%duplicate = b%duplicate
    a%triangle  = b%triangle
    a%unitd     = b%unitd
    a%upper     = .not.b%upper
    a%sorted    = .false.
    if (allocated(b%aux)) then 
      allocate(a%aux(size(b%aux)))
      a%aux(:)       = b%aux(:)
    end if
    
  end subroutine psb_base_transp_2mat

  subroutine psb_base_transc_2mat(a,b)
    implicit none 
    
    class(psb_base_sparse_mat), intent(out) :: a
    class(psb_base_sparse_mat), intent(in)   :: b

    call a%transp(b) 
  end subroutine psb_base_transc_2mat

  subroutine psb_base_transp_1mat(a)
    implicit none 
    
    class(psb_base_sparse_mat), intent(inout) :: a
    integer :: itmp

    itmp        = a%m
    a%m         = a%n
    a%n         = itmp
    a%state     = a%state
    a%duplicate = a%duplicate
    a%triangle  = a%triangle
    a%unitd     = a%unitd
    a%upper     = .not.a%upper
    a%sorted    = .false.
    
  end subroutine psb_base_transp_1mat

  subroutine psb_base_transc_1mat(a)
    implicit none 
    
    class(psb_base_sparse_mat), intent(inout) :: a
    
    call a%transp() 
  end subroutine psb_base_transc_1mat


end module psb_base_mat_mod

