module psb_c_csc_mat_mod

  use psb_c_base_mat_mod

  type, extends(psb_c_base_sparse_mat) :: psb_c_csc_sparse_mat

    integer, allocatable :: icp(:), ia(:)
    complex(psb_spk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_nzeros => c_csc_get_nzeros
    procedure, pass(a) :: get_fmt  => c_csc_get_fmt
    procedure, pass(a) :: get_diag => c_csc_get_diag
    procedure, pass(a) :: c_base_csmm => c_csc_csmm
    procedure, pass(a) :: c_base_csmv => c_csc_csmv
    procedure, pass(a) :: c_base_cssm => c_csc_cssm
    procedure, pass(a) :: c_base_cssv => c_csc_cssv
    procedure, pass(a) :: c_scals => c_csc_scals
    procedure, pass(a) :: c_scal => c_csc_scal
    procedure, pass(a) :: csnmi => c_csc_csnmi
    procedure, pass(a) :: reallocate_nz => c_csc_reallocate_nz
    procedure, pass(a) :: csput => c_csc_csput
    procedure, pass(a) :: allocate_mnnz => c_csc_allocate_mnnz
    procedure, pass(a) :: cp_to_coo => c_cp_csc_to_coo
    procedure, pass(a) :: cp_from_coo => c_cp_csc_from_coo
    procedure, pass(a) :: cp_to_fmt => c_cp_csc_to_fmt
    procedure, pass(a) :: cp_from_fmt => c_cp_csc_from_fmt
    procedure, pass(a) :: mv_to_coo => c_mv_csc_to_coo
    procedure, pass(a) :: mv_from_coo => c_mv_csc_from_coo
    procedure, pass(a) :: mv_to_fmt => c_mv_csc_to_fmt
    procedure, pass(a) :: mv_from_fmt => c_mv_csc_from_fmt
    procedure, pass(a) :: csgetptn => c_csc_csgetptn
    procedure, pass(a) :: c_csgetrow => c_csc_csgetrow
    procedure, pass(a) :: get_size => c_csc_get_size
    procedure, pass(a) :: free => c_csc_free
    procedure, pass(a) :: trim => c_csc_trim
    procedure, pass(a) :: print => c_csc_print
    procedure, pass(a) :: sizeof => c_csc_sizeof
    procedure, pass(a) :: reinit => c_csc_reinit
    procedure, pass(a) :: c_csc_cp_from
    generic, public    :: cp_from => c_csc_cp_from
    procedure, pass(a) :: c_csc_mv_from
    generic, public    :: mv_from => c_csc_mv_from
  end type psb_c_csc_sparse_mat

  private :: c_csc_get_nzeros, c_csc_csmm, c_csc_csmv, c_csc_cssm, c_csc_cssv, &
       & c_csc_csput, c_csc_reallocate_nz, c_csc_allocate_mnnz, &
       & c_csc_free,  c_csc_print, c_csc_get_fmt, c_csc_csnmi, get_diag, &
       & c_cp_csc_to_coo, c_cp_csc_from_coo, &
       & c_mv_csc_to_coo, c_mv_csc_from_coo, &
       & c_cp_csc_to_fmt, c_cp_csc_from_fmt, &
       & c_mv_csc_to_fmt, c_mv_csc_from_fmt, &
       & c_csc_scals, c_csc_scal, c_csc_trim, c_csc_csgetrow, c_csc_get_size, &
       & c_csc_sizeof, c_csc_csgetptn, c_csc_get_nz_row, c_csc_reinit


  interface 
    subroutine c_cp_csc_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      class(psb_c_base_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine c_cp_csc_to_fmt_impl
  end interface

  interface 
    subroutine c_cp_csc_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(in) :: b
      integer, intent(out)                        :: info
    end subroutine c_cp_csc_from_fmt_impl
  end interface


  interface 
    subroutine c_cp_csc_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      class(psb_c_coo_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine c_cp_csc_to_coo_impl
  end interface

  interface 
    subroutine c_cp_csc_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine c_cp_csc_from_coo_impl
  end interface

  interface 
    subroutine c_mv_csc_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_base_sparse_mat), intent(out)  :: b
      integer, intent(out)            :: info
    end subroutine c_mv_csc_to_fmt_impl
  end interface

  interface 
    subroutine c_mv_csc_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout)  :: a
      class(psb_c_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine c_mv_csc_from_fmt_impl
  end interface


  interface 
    subroutine c_mv_csc_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(out)   :: b
      integer, intent(out)            :: info
    end subroutine c_mv_csc_to_coo_impl
  end interface

  interface 
    subroutine c_mv_csc_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psb_c_base_mat_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      class(psb_c_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine c_mv_csc_from_coo_impl
  end interface

  interface 
    subroutine c_csc_csput_impl(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(inout) :: a
      complex(psb_spk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine c_csc_csput_impl
  end interface

  interface 
    subroutine c_csc_csgetptn_impl(imin,imax,a,nz,ia,ja,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      use psb_const_mod
      import psb_c_csc_sparse_mat
      implicit none
      
      class(psb_c_csc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine c_csc_csgetptn_impl
  end interface

  interface 
    subroutine c_csc_csgetrow_impl(imin,imax,a,nz,ia,ja,val,info,&
         & jmin,jmax,iren,append,nzin,rscale,cscale)
      use psb_const_mod
      import psb_c_csc_sparse_mat
      implicit none
      
      class(psb_c_csc_sparse_mat), intent(in) :: a
      integer, intent(in)                  :: imin,imax
      integer, intent(out)                 :: nz
      integer, allocatable, intent(inout)  :: ia(:), ja(:)
      complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
      integer,intent(out)                  :: info
      logical, intent(in), optional        :: append
      integer, intent(in), optional        :: iren(:)
      integer, intent(in), optional        :: jmin,jmax, nzin
      logical, intent(in), optional        :: rscale,cscale
    end subroutine c_csc_csgetrow_impl
  end interface

  interface c_csc_cssm_impl
    subroutine c_csc_cssv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine c_csc_cssv_impl
    subroutine c_csc_cssm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine c_csc_cssm_impl
  end interface

  interface c_csc_csmm_impl
    subroutine c_csc_csmv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
      complex(psb_spk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine c_csc_csmv_impl
    subroutine c_csc_csmm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
      complex(psb_spk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine c_csc_csmm_impl
  end interface

  interface c_csc_csnmi_impl
    function c_csc_csnmi_impl(a) result(res)
      use psb_const_mod
      import psb_c_csc_sparse_mat
      class(psb_c_csc_sparse_mat), intent(in) :: a
      real(psb_spk_)         :: res
    end function c_csc_csnmi_impl
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

  
  function c_csc_sizeof(a) result(res)
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res = 8 
    res = res + 2*psb_sizeof_sp  * size(a%val)
    res = res + psb_sizeof_int * size(a%icp)
    res = res + psb_sizeof_int * size(a%ia)
      
  end function c_csc_sizeof

  function c_csc_get_fmt(a) result(res)
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'CSC'
  end function c_csc_get_fmt
  
  function c_csc_get_nzeros(a) result(res)
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer :: res
    res = a%icp(a%get_ncols()+1)-1
  end function c_csc_get_nzeros

  function c_csc_get_size(a) result(res)
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer :: res

    res = -1
    
    if (allocated(a%ia)) then 
      if (res >= 0) then 
        res = min(res,size(a%ia))
      else 
        res = size(a%ia)
      end if
    end if
    if (allocated(a%val)) then 
      if (res >= 0) then 
        res = min(res,size(a%val))
      else 
        res = size(a%val)
      end if
    end if

  end function c_csc_get_size



  function  c_csc_get_nz_col(idx,a) result(res)
    use psb_const_mod
    implicit none
    
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    
    res = 0 
 
    if ((1<=idx).and.(idx<=a%get_ncols())) then 
      res = a%icp(idx+1)-a%icp(idx)
    end if
    
  end function c_csc_get_nz_col



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


  subroutine  c_csc_reallocate_nz(nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: nz
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='c_csc_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    call psb_realloc(nz,a%ia,info)
    if (info == 0) call psb_realloc(nz,a%val,info)
    if (info == 0) call psb_realloc(max(nz,a%get_nrows()+1,a%get_ncols()+1),a%icp,info)
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

  end subroutine c_csc_reallocate_nz

  subroutine c_csc_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_const_mod
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    complex(psb_spk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)


    Integer            :: err_act
    character(len=20)  :: name='c_csc_csput'
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

    call c_csc_csput_impl(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
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
  end subroutine c_csc_csput

  subroutine c_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_c_csc_sparse_mat), intent(in) :: a
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

    call c_csc_csgetptn_impl(imin,imax,a,nz,ia,ja,info,&
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

  end subroutine c_csc_csgetptn


  subroutine c_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_c_csc_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: imin,imax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
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

    call c_csc_csgetrow_impl(imin,imax,a,nz,ia,ja,val,info,&
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

  end subroutine c_csc_csgetrow


  subroutine c_csc_csgetblk(imin,imax,a,b,info,&
       & jmin,jmax,iren,append,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_c_csc_sparse_mat), intent(in) :: a
    class(psb_c_coo_sparse_mat), intent(inout) :: b
    integer, intent(in)                  :: imin,imax
    integer,intent(out)                  :: info
    logical, intent(in), optional        :: append
    integer, intent(in), optional        :: iren(:)
    integer, intent(in), optional        :: jmin,jmax
    logical, intent(in), optional        :: rscale,cscale
    Integer :: err_act, nzin, nzout
    character(len=20)  :: name='csget'
    logical :: append_
    logical, parameter :: debug=.false.
    
    call psb_erractionsave(err_act)
    info = 0

    if (present(append)) then 
      append_ = append
    else
      append_ = .false.
    endif
    if (append_) then 
      nzin = a%get_nzeros()
    else
      nzin = 0
    endif

    call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
         & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
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

  end subroutine c_csc_csgetblk


  subroutine c_csc_csclip(a,b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
    ! Output is always in  COO format 
    use psb_error_mod
    use psb_const_mod
    implicit none
    
    class(psb_c_csc_sparse_mat), intent(in) :: a
    class(psb_c_coo_sparse_mat), intent(out) :: b
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
    
  end subroutine c_csc_csclip


  subroutine  c_csc_free(a) 
    implicit none 

    class(psb_c_csc_sparse_mat), intent(inout) :: a

    if (allocated(a%icp)) deallocate(a%icp)
    if (allocated(a%ia)) deallocate(a%ia)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine c_csc_free

  subroutine c_csc_reinit(a,clear)
    use psb_error_mod
    implicit none 

    class(psb_c_csc_sparse_mat), intent(inout) :: a   
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

  end subroutine c_csc_reinit

  subroutine  c_csc_trim(a)
    use psb_realloc_mod
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    Integer :: err_act, info, nz, n
    character(len=20)  :: name='trim'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    n   = a%get_ncols()
    nz  = a%get_nzeros()
    if (info == 0) call psb_realloc(n+1,a%icp,info)
    if (info == 0) call psb_realloc(nz,a%ia,info)
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

  end subroutine c_csc_trim

  subroutine c_cp_csc_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    class(psb_c_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_cp_csc_to_coo_impl(a,b,info)
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

  end subroutine c_cp_csc_to_coo
  
  subroutine c_cp_csc_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    class(psb_c_coo_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_cp_csc_from_coo_impl(a,b,info)
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

  end subroutine c_cp_csc_from_coo


  subroutine c_cp_csc_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    class(psb_c_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_cp_csc_to_fmt_impl(a,b,info)
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

  end subroutine c_cp_csc_to_fmt
  
  subroutine c_cp_csc_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    class(psb_c_base_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_cp_csc_from_fmt_impl(a,b,info)
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

  end subroutine c_cp_csc_from_fmt


  subroutine c_mv_csc_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    class(psb_c_coo_sparse_mat), intent(out)   :: b
    integer, intent(out)            :: info 

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_mv_csc_to_coo_impl(a,b,info)
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

  end subroutine c_mv_csc_to_coo
  
  subroutine c_mv_csc_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    class(psb_c_coo_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_mv_csc_from_coo_impl(a,b,info)
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

  end subroutine c_mv_csc_from_coo


  subroutine c_mv_csc_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    class(psb_c_base_sparse_mat), intent(out)  :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_mv_csc_to_fmt_impl(a,b,info)
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

  end subroutine c_mv_csc_to_fmt
  
  subroutine c_mv_csc_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout)  :: a
    class(psb_c_base_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call c_mv_csc_from_fmt_impl(a,b,info)
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

  end subroutine c_mv_csc_from_fmt

  subroutine  c_csc_allocate_mnnz(m,n,a,nz) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psb_c_csc_sparse_mat), intent(inout) :: a
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
      
    if (info == 0) call psb_realloc(n+1,a%icp,info)
    if (info == 0) call psb_realloc(nz_,a%ia,info)
    if (info == 0) call psb_realloc(nz_,a%val,info)
    if (info == 0) then 
      a%icp=0
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

  end subroutine c_csc_allocate_mnnz


  subroutine c_csc_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_string_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psb_c_csc_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act
    character(len=20)  :: name='c_csc_print'
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
        do j=a%icp(i),a%icp(i+1)-1 
          write(iout,frmtv) iv(a%ia(j)),iv(i),a%val(j)
        end do
      enddo
    else      
      if (present(ivr).and..not.present(ivc)) then 
        do i=1, nr
          do j=a%icp(i),a%icp(i+1)-1 
            write(iout,frmtv) ivr(a%ia(j)),i,a%val(j)
          end do
        enddo
      else if (present(ivr).and.present(ivc)) then 
        do i=1, nr
          do j=a%icp(i),a%icp(i+1)-1 
            write(iout,frmtv) ivr(a%ia(j)),ivc(i),a%val(j)
          end do
        enddo
      else if (.not.present(ivr).and.present(ivc)) then 
        do i=1, nr
          do j=a%icp(i),a%icp(i+1)-1 
            write(iout,frmtv) (a%ia(j)),ivc(i),a%val(j)
          end do
        enddo
      else if (.not.present(ivr).and..not.present(ivc)) then 
        do i=1, nr
          do j=a%icp(i),a%icp(i+1)-1 
            write(iout,frmtv) (a%ia(j)),(i),a%val(j)
          end do
        enddo
      endif
    endif

  end subroutine c_csc_print


  subroutine c_csc_cp_from(a,b)
    use psb_error_mod
    implicit none 

    class(psb_c_csc_sparse_mat), intent(out) :: a
    type(psb_c_csc_sparse_mat), intent(in)   :: b


    Integer :: err_act, info
    character(len=20)  :: name='cp_from'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    info = 0

    call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
    call a%psb_c_base_sparse_mat%cp_from(b%psb_c_base_sparse_mat)
    a%icp = b%icp 
    a%ia  = b%ia
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

  end subroutine c_csc_cp_from

  subroutine c_csc_mv_from(a,b)
    use psb_error_mod
    implicit none 

    class(psb_c_csc_sparse_mat), intent(out)  :: a
    type(psb_c_csc_sparse_mat), intent(inout) :: b
    

    Integer :: err_act, info
    character(len=20)  :: name='mv_from'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call a%psb_c_base_sparse_mat%mv_from(b%psb_c_base_sparse_mat)
    call move_alloc(b%icp, a%icp)
    call move_alloc(b%ia,  a%ia)
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

  end subroutine c_csc_mv_from



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


  subroutine c_csc_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
    complex(psb_spk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    Integer :: err_act
    character(len=20)  :: name='c_csc_csmv'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call c_csc_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine c_csc_csmv

  subroutine c_csc_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
    complex(psb_spk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='c_csc_csmm'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)



    call c_csc_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine c_csc_csmm


  subroutine c_csc_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
    complex(psb_spk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    Integer :: err_act
    character(len=20)  :: name='c_csc_cssv'
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
    
    call c_csc_cssm_impl(alpha,a,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine c_csc_cssv



  subroutine c_csc_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
    complex(psb_spk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    Integer :: err_act
    character(len=20)  :: name='c_csc_csmm'
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
    
    call c_csc_cssm_impl(alpha,a,x,beta,y,info,trans) 
    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine c_csc_cssm
 
  function c_csc_csnmi(a) result(res)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    real(psb_spk_)         :: res
    
    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.
    
   
    res = c_csc_csnmi_impl(a)
    
    return

  end function c_csc_csnmi

  subroutine c_csc_get_diag(a,d,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(out)     :: d(:)
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
      do k=a%icp(i),a%icp(i+1)-1
        j=a%ia(k)
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

  end subroutine c_csc_get_diag


  subroutine c_csc_scal(d,a,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    complex(psb_spk_), intent(in)      :: d(:)
    integer, intent(out)            :: info

    Integer :: err_act,mnm, i, j, n
    character(len=20)  :: name='scal'
    logical, parameter :: debug=.false.

    info  = 0
    call psb_erractionsave(err_act)

    n = a%get_ncols()
    if (size(d) < n) then 
      info=35
      call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
      goto 9999
    end if

    do i=1, n
      do j = a%icp(i), a%icp(i+1) -1 
        a%val(j) = a%val(j) * d(a%ia(j))
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

  end subroutine c_csc_scal


  subroutine c_csc_scals(d,a,info) 
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psb_c_csc_sparse_mat), intent(inout) :: a
    complex(psb_spk_), intent(in)      :: d
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

  end subroutine c_csc_scals



end module psb_c_csc_mat_mod
