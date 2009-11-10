module psb_s_csr_mat_mod

  use psb_s_base_mat_mod

  type, extends(psb_s_base_sparse_mat) :: psb_s_csr_sparse_mat

    integer, allocatable :: irp(:), ja(:)
    real(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_nzeros => s_csr_get_nzeros
    procedure, pass(a) :: get_fmt  => s_csr_get_fmt
    procedure, pass(a) :: get_diag => s_csr_get_diag
    procedure, pass(a) :: s_base_csmm => s_csr_csmm
    procedure, pass(a) :: s_base_csmv => s_csr_csmv
    procedure, pass(a) :: s_base_cssm => s_csr_cssm
    procedure, pass(a) :: s_base_cssv => s_csr_cssv
    procedure, pass(a) :: s_scals => s_csr_scals
    procedure, pass(a) :: s_scal => s_csr_scal
    procedure, pass(a) :: csnmi => s_csr_csnmi
    procedure, pass(a) :: reallocate_nz => s_csr_reallocate_nz
    procedure, pass(a) :: csput => s_csr_csput
    procedure, pass(a) :: allocate_mnnz => s_csr_allocate_mnnz
    procedure, pass(a) :: cp_to_coo => s_cp_csr_to_coo
    procedure, pass(a) :: cp_from_coo => s_cp_csr_from_coo
    procedure, pass(a) :: cp_to_fmt => s_cp_csr_to_fmt
    procedure, pass(a) :: cp_from_fmt => s_cp_csr_from_fmt
    procedure, pass(a) :: mv_to_coo => s_mv_csr_to_coo
    procedure, pass(a) :: mv_from_coo => s_mv_csr_from_coo
    procedure, pass(a) :: mv_to_fmt => s_mv_csr_to_fmt
    procedure, pass(a) :: mv_from_fmt => s_mv_csr_from_fmt
    procedure, pass(a) :: csgetptn => s_csr_csgetptn
    procedure, pass(a) :: s_csgetrow => s_csr_csgetrow
    procedure, pass(a) :: get_nz_row => s_csr_get_nz_row
    procedure, pass(a) :: get_size => s_csr_get_size
    procedure, pass(a) :: free => s_csr_free
    procedure, pass(a) :: trim => s_csr_trim
    procedure, pass(a) :: print => s_csr_print
    procedure, pass(a) :: sizeof => s_csr_sizeof
    procedure, pass(a) :: reinit => s_csr_reinit
    procedure, pass(a) :: s_csr_cp_from
    generic, public    :: cp_from => s_csr_cp_from
    procedure, pass(a) :: s_csr_mv_from
    generic, public    :: mv_from => s_csr_mv_from

  end type psb_s_csr_sparse_mat

  private :: s_csr_get_nzeros, s_csr_csmm, s_csr_csmv, s_csr_cssm, s_csr_cssv, &
       & s_csr_csput, s_csr_reallocate_nz, s_csr_allocate_mnnz, &
       & s_csr_free,  s_csr_print, s_csr_get_fmt, s_csr_csnmi, get_diag, &
       & s_cp_csr_to_coo, s_cp_csr_from_coo, &
       & s_mv_csr_to_coo, s_mv_csr_from_coo, &
       & s_cp_csr_to_fmt, s_cp_csr_from_fmt, &
       & s_mv_csr_to_fmt, s_mv_csr_from_fmt, &
       & s_csr_scals, s_csr_scal, s_csr_trim, s_csr_csgetrow, s_csr_get_size, &
       & s_csr_sizeof, s_csr_csgetptn, s_csr_get_nz_row, s_csr_reinit
!!$, &
!!$       & s_csr_mv_from, s_csr_mv_from


  interface 
    subroutine s_cp_csr_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      class(psb_s_base_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine s_cp_csr_to_fmt_impl
  end interface

  interface 
    subroutine s_cp_csr_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(in) :: b
      integer, intent(out)                        :: info
    end subroutine s_cp_csr_from_fmt_impl
  end interface


  interface 
    subroutine s_cp_csr_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      class(psb_s_coo_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine s_cp_csr_to_coo_impl
  end interface

  interface 
    subroutine s_cp_csr_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine s_cp_csr_from_coo_impl
  end interface

  interface 
    subroutine s_mv_csr_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      class(psb_s_base_sparse_mat), intent(out)  :: b
      integer, intent(out)            :: info
    end subroutine s_mv_csr_to_fmt_impl
  end interface

  interface 
    subroutine s_mv_csr_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout)  :: a
      class(psb_s_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine s_mv_csr_from_fmt_impl
  end interface


  interface 
    subroutine s_mv_csr_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(out)   :: b
      integer, intent(out)            :: info
    end subroutine s_mv_csr_to_coo_impl
  end interface

  interface 
    subroutine s_mv_csr_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_s_base_mat_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      class(psb_s_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine s_mv_csr_from_coo_impl
  end interface

  interface 
    subroutine s_csr_csput_impl(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(inout) :: a
      real(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine s_csr_csput_impl
  end interface

  interface 
    subroutine s_csr_csgetptn_impl(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      use psb_const_mod
      import psb_s_csr_sparse_mat
      implicit none
      
      class(psb_s_csr_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine s_csr_csgetptn_impl
  end interface

  interface 
    subroutine s_csr_csgetrow_impl(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      use psb_const_mod
      import psb_s_csr_sparse_mat
      implicit none
      
      class(psb_s_csr_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      real(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine s_csr_csgetrow_impl
  end interface

  interface s_csr_cssm_impl
    subroutine s_csr_cssv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine s_csr_cssv_impl
    subroutine s_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine s_csr_cssm_impl
  end interface

  interface s_csr_csmm_impl
    subroutine s_csr_csmv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:)
      real(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine s_csr_csmv_impl
    subroutine s_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine s_csr_csmm_impl
  end interface

  interface s_csr_csnmi_impl
    function s_csr_csnmi_impl(a) result(res)
      use psb_const_mod
      import psb_s_csr_sparse_mat
      class(psb_s_csr_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function s_csr_csnmi_impl
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

  
  function s_csr_sizeof(a) result(res)
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + psb_sizeof_sp  * size(a%val)
    res = res + psb_sizeof_int * size(a%irp)
    res = res + psb_sizeof_int * size(a%ja)
      
  end function s_csr_sizeof

  function s_csr_get_fmt(a) result(res)
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'CSR'
  end function s_csr_get_fmt
  
  function s_csr_get_nzeros(a) result(res)
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    integer :: res
    res = a%irp(a%get_nrows()+1)-1
  end function s_csr_get_nzeros

  function s_csr_get_size(a) result(res)
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    integer :: res

    res = -1
    
    if (allocated(a%ja)) then 
      if (res >= 0) then 
        res = min(res,size(a%ja))
      else 
        res = size(a%ja)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function s_csr_get_size



  function  s_csr_get_nz_row(idx,a) result(res)
    use psb_const_mod
    implicit none
    
    class(psb_s_csr_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_nrows())) then 
      res = a%irp(idx+1)-a%irp(idx)
    end if
    
  end function s_csr_get_nz_row



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


  subroutine  s_csr_reallocate_nz(nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: nz
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='s_csr_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    call psb_realloc(nz,a%ja,info)
    if (info == 0) call psb_realloc(nz,a%val,info)
    if (info == 0) call psb_realloc(&
         & max(nz,a%get_nrows()+1,a%get_ncols()+1),a%irp,info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_reallocate_nz

  subroutine s_csr_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_const_mod
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    real(psb_spk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)


    Integer            :: err_act
    character(len=20)  :: name='s_csr_csput'
    logical, parameter :: debug=.false.
    integer            :: nza, i,j,k, nzl, isza, int_err(5)

    call psb_erractionsave(err_act)
    info = 0

    if (nz <= 0) then 
      info = 10
      int_err(1)=1
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (size(ia) < nz) then 
      info = 35
      int_err(1)=2
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    if (size(ja) < nz) then 
      info = 35
      int_err(1)=3
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
    if (size(val) < nz) then 
      info = 35
      int_err(1)=4
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    if (nz == 0) return

    call s_csr_csput_impl(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
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
  end subroutine s_csr_csput

  subroutine s_csr_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_s_csr_sparse_mat), intent(in) :: a
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

    call psb_erractionsave(err_act)
    info = 0

    call s_csr_csgetptn_impl(imin,imax,a,nz,ia,ja,info,&
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
    return

  end subroutine s_csr_csgetptn


  subroutine s_csr_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_s_csr_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: imin,imax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    real(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax, nzin
    logical, intent(in), optional        :: rscale,cscale
    Integer :: err_act
    character(len=20)  :: name='csget'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0

    call s_csr_csgetrow_impl(imin,imax,a,nz,ia,ja,val,info,&
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
    return

  end subroutine s_csr_csgetrow


  subroutine s_csr_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_s_csr_sparse_mat), intent(in) :: a
    class(psb_s_coo_sparse_mat), intent(inout) :: b
    integer, intent(in)                  :: imin,imax
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax
    logical, intent(in), optional        :: rscale,cscale
    Integer :: err_act, nzin, nzout
    character(len=20)  :: name='csget'
    logical :: appens_
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0

    if (present(append)) then 
      appens_ = append
    else
      appens_ = .false.
    endif
    if (appens_) then 
      nzin = a%get_nzeros()
    else
      nzin = 0
    endif

    call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
         & jmin=jmin, jmax=jmax, iren=iren, append=appens_, &
         & nzin=nzin, rscale=rscale, cscale=cscale)

    if (info /= 0) goto 9999

    call b%set_nzeros(nzin+nzout)
    call b%fix(info)
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

  end subroutine s_csr_csgetblk


  subroutine s_csr_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_s_csr_sparse_mat), intent(in) :: a
    class(psb_s_coo_sparse_mat), intent(out) :: b
    integer,intent(out)                  :: info
    integer, intent(in), optional        :: imin,imax,jmin,jmax
    logical, intent(in), optional        :: rscale,cscale

    Integer :: err_act, nzin, nzout, imin_, imax_, jmin_, jmax_, mb,nb
    character(len=20)  :: name='csget'
    logical :: rscale_, cscale_
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0

    nzin = 0
    if (present(imin)) then 
      imin_ = imin
    else
      imin_ = 1
    end if
    if (present(imax)) then 
      imax_ = imax
    else
      imax_ = a%get_nrows()
    end if
    if (present(jmin)) then 
      jmin_ = jmin
    else
      jmin_ = 1
    end if
    if (present(jmax)) then 
      jmax_ = jmax
    else
      jmax_ = a%get_ncols()
    end if
    if (present(rscale)) then 
      rscale_ = rscale
    else
      rscale_ = .true.
    end if
    if (present(cscale)) then 
      cscale_ = cscale
    else
      cscale_ = .true.
    end if

    if (rscale_) then 
      mb = imax_ - imin_ +1
    else 
      mb = a%get_nrows() ! Should this be imax_ ?? 
    endif
    if (cscale_) then 
      nb = jmax_ - jmin_ +1
    else 
      nb = a%get_ncols()  ! Should this be jmax_ ?? 
    endif
    call b%allocate(mb,nb)

    call a%csget(imin_,imax_,nzout,b%ia,b%ja,b%val,info,&
         & jmin=jmin_, jmax=jmax_, append=.false., &
         & nzin=nzin, rscale=rscale_, cscale=cscale_)

    if (info /= 0) goto 9999

    call b%set_nzeros(nzin+nzout)
    call b%fix(info)

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
    
  end subroutine s_csr_csclip


  subroutine  s_csr_free(a) 
    implicit none 

    class(psb_s_csr_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine s_csr_free

  subroutine s_csr_reinit(a,clear)
    use psb_error_mod
    implicit none 

    class(psb_s_csr_sparse_mat), intent(inout) :: a   
    logical, intent(in), optional :: clear

    Integer :: err_act, info
    character(len=20)  :: name='reinit'
    logical  :: clear_
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0


    if (present(clear)) then 
      clear_ = clear
    else
      clear_ = .true.
    end if

    if (a%is_bld() .or. a%is_upd()) then 
      ! do nothing
      return
    else if (a%is_asb()) then 
      if (clear_) a%val(:) = szero
      call a%set_upd()
    else
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_reinit


  subroutine  s_csr_trim(a)
    use psb_realloc_mod
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    Integer :: err_act, info, nz, m 
    character(len=20)  :: name='trim'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    m   = a%get_nrows()
    nz  = a%get_nzeros()
    if (info == 0) call psb_realloc(m+1,a%irp,info)
    if (info == 0) call psb_realloc(nz,a%ja,info)
    if (info == 0) call psb_realloc(nz,a%val,info)

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

  end subroutine s_csr_trim


  subroutine s_cp_csr_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    class(psb_s_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_cp_csr_to_coo_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_cp_csr_to_coo
  
  subroutine s_cp_csr_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    class(psb_s_coo_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_cp_csr_from_coo_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_cp_csr_from_coo


  subroutine s_cp_csr_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    class(psb_s_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_cp_csr_to_fmt_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_cp_csr_to_fmt
  
  subroutine s_cp_csr_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    class(psb_s_base_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_cp_csr_from_fmt_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_cp_csr_from_fmt


  subroutine s_mv_csr_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    class(psb_s_coo_sparse_mat), intent(out)   :: b
    integer, intent(out)            :: info 

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_mv_csr_to_coo_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_mv_csr_to_coo
  
  subroutine s_mv_csr_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    class(psb_s_coo_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_mv_csr_from_coo_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_mv_csr_from_coo


  subroutine s_mv_csr_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    class(psb_s_base_sparse_mat), intent(out)  :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_mv_csr_to_fmt_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_mv_csr_to_fmt
  
  subroutine s_mv_csr_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout)  :: a
    class(psb_s_base_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call s_mv_csr_from_fmt_impl(a,b,info)
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)
          
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_mv_csr_from_fmt


  subroutine  s_csr_allocate_mnnz(m,n,a,nz) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    integer, intent(in), optional :: nz
    Integer :: err_act, info, nz_
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    if (m < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
      goto 9999
    endif
    if (n < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
      goto 9999
    endif
        if (present(nz)) then 
      nz_ = nz
    else
      nz_ = max(7*m,7*n,1)
    end if
    if (nz_ < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
      goto 9999
    endif
      
    if (info == 0) call psb_realloc(m+1,a%irp,info)
    if (info == 0) call psb_realloc(nz_,a%ja,info)
    if (info == 0) call psb_realloc(nz_,a%val,info)
    if (info == 0) then 
      a%irp=0
      call a%set_nrows(m)
      call a%set_ncols(n)
      call a%set_bld()
      call a%set_triangle(.false.)
      call a%set_unit(.false.)
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_allocate_mnnz


  subroutine s_csr_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_string_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psb_s_csr_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act
    character(len=20)  :: name='s_csr_print'
    logical, parameter :: debug=.false.

    character(len=80)                 :: frmtv 
    integer  :: irs,ics,i,j, nmx, ni, nr, nc, nz
    
    if (present(eirs)) then 
      irs = eirs
    else
      irs = 0
    endif
    if (present(eics)) then 
      ics = eics
    else
      ics = 0
    endif
    
    if (present(head)) then 
      write(iout,'(a)') '%%MatrixMarket matrix coordinate real general'
      write(iout,'(a,a)') '% ',head 
      write(iout,'(a)') '%'    
      write(iout,'(a,a)') '% COO'
    endif
    
    nr = a%get_nrows()
    nc = a%get_ncols()
    nz = a%get_nzeros()
    nmx = max(nr,nc,1)
    ni  = floor(log10(1.0*nmx)) + 1
    
    write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
    write(iout,*) nr, nc, nz 
    if(present(iv)) then 
      do i=1, nr
        do j=a%irp(i),a%irp(i+1)-1 
          write(iout,frmtv) iv(i),iv(a%ja(j)),a%val(j)
        end do
      enddo
    else      
      if (present(ivr).and..not.present(ivc)) then 
        do i=1, nr
          do j=a%irp(i),a%irp(i+1)-1 
            write(iout,frmtv) ivr(i),(a%ja(j)),a%val(j)
          end do
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do i=1, nr
          do j=a%irp(i),a%irp(i+1)-1 
            write(iout,frmtv) ivr(i),ivc(a%ja(j)),a%val(j)
          end do
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do i=1, nr
          do j=a%irp(i),a%irp(i+1)-1 
            write(iout,frmtv) (i),ivc(a%ja(j)),a%val(j)
          end do
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do i=1, nr
          do j=a%irp(i),a%irp(i+1)-1 
            write(iout,frmtv) (i),(a%ja(j)),a%val(j)
          end do
        enddo
      endif
    endif

  end subroutine s_csr_print


  subroutine s_csr_cp_from(a,b)
    use psb_error_mod
    implicit none 

    class(psb_s_csr_sparse_mat), intent(out) :: a
    type(psb_s_csr_sparse_mat), intent(in)   :: b


    Integer :: err_act, info
    character(len=20)  :: name='cp_from'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    info = 0

    call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
    call a%psb_s_base_sparse_mat%cp_from(b%psb_s_base_sparse_mat)
    a%irp = b%irp 
    a%ja  = b%ja
    a%val = b%val 

    if (info /= 0) goto 9999
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_csr_cp_from

  subroutine s_csr_mv_from(a,b)
    use psb_error_mod
    implicit none 

    class(psb_s_csr_sparse_mat), intent(out)  :: a
    type(psb_s_csr_sparse_mat), intent(inout) :: b
    

    Integer :: err_act, info
    character(len=20)  :: name='mv_from'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call a%psb_s_base_sparse_mat%mv_from(b%psb_s_base_sparse_mat)
    call move_alloc(b%irp, a%irp)
    call move_alloc(b%ja,  a%ja)
    call move_alloc(b%val, a%val)
    call b%free()

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    call psb_errpush(info,name)

    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return

  end subroutine s_csr_mv_from



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


  subroutine s_csr_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(in)          :: alpha, beta, x(:)
    real(psb_spk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_spk_) :: acc
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='s_csr_csmv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call s_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine s_csr_csmv

  subroutine s_csr_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_spk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_spk_), allocatable  :: acc(:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='s_csr_csmm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)



    call s_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine s_csr_csmm


  subroutine s_csr_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(in)          :: alpha, beta, x(:)
    real(psb_spk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_spk_) :: acc
    real(psb_spk_), allocatable :: tmp(:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='s_csr_cssv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    
    if (.not. (a%is_triangle())) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call s_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine s_csr_cssv



  subroutine s_csr_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_spk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_spk_) :: acc
    real(psb_spk_), allocatable :: tmp(:,:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='s_csr_csmm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif

    
    if (.not. (a%is_triangle())) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call s_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 
    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_cssm
 
  function s_csr_csnmi(a) result(res)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_)         :: res
    
    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.
    
   
    res = s_csr_csnmi_impl(a)
    
    return

  end function s_csr_csnmi

  subroutine s_csr_get_diag(a,d,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(out)     :: d(:)
    integer, intent(out)            :: info

    Integer :: err_act, mnm, i, j, k
    character(len=20)  :: name='get_diag'
    logical, parameter :: debug=.false.

    info  = 0
    call psb_erractionsave(err_act)

    mnm = min(a%get_nrows(),a%get_ncols())
    if (size(d) < mnm) then 
      info=35
      call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
      goto 9999
    end if


    do i=1, mnm
      do k=a%irp(i),a%irp(i+1)-1
        j=a%ja(k)
        if ((j==i) .and.(j <= mnm )) then 
          d(i) = a%val(k)
        endif
      enddo
    end do
    do i=mnm+1,size(d) 
      d(i) = szero 
    end do
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_get_diag


  subroutine s_csr_scal(d,a,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    real(psb_spk_), intent(in)      :: d(:)
    integer, intent(out)            :: info

    Integer :: err_act,mnm, i, j, m
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.

    info  = 0
    call psb_erractionsave(err_act)

    m = a%get_nrows()
    if (size(d) < m) then 
      info=35
      call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
      goto 9999
    end if

    do i=1, m 
      do j = a%irp(i), a%irp(i+1) -1 
        a%val(j) = a%val(j) * d(i)
      end do
    enddo

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_scal


  subroutine s_csr_scals(d,a,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_s_csr_sparse_mat), intent(inout) :: a
    real(psb_spk_), intent(in)      :: d
    integer, intent(out)            :: info

    Integer :: err_act,mnm, i, j, m
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.

    info  = 0
    call psb_erractionsave(err_act)


    do i=1,a%get_nzeros()
      a%val(i) = a%val(i) * d
    enddo

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_csr_scals



end module psb_s_csr_mat_mod