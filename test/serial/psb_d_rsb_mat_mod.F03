!
!
! FIXME:
! * some RSB constants are used in their value form, and with no explanation
! * error handling
! * PSBLAS interface adherence
! * should test and fix all the problems that for sure will occur
! * ..
! 
module psb_d_rsb_mat_mod
  use psb_d_base_mat_mod
  use rsb_mod
#ifdef HAVE_LIBRSB
    use iso_c_binding  
#endif
  integer :: c_f_order=2        ! FIXME: here should use RSB_FLAG_WANT_COLUMN_MAJOR_ORDER
  integer :: c_f_index=256*16 ! 0x001000 ! FIXME: here should use RSB_FLAG_FORTRAN_INDICES_INTERFACE
  integer :: c_d_typecode=0 ! FIXME: here should use ..
  integer :: c_def_flags =0 ! FIXME: here should use ..
  integer :: c_upd_flags =0 ! FIXME: here should use ..
  type, extends(psb_d_base_sparse_mat) :: psb_d_rsb_sparse_mat
#ifdef HAVE_LIBRSB
    type(c_ptr) :: rsbmptr
    contains 
    procedure, pass(a) :: get_size     => d_rsb_get_size
    procedure, pass(a) :: get_nzeros   => d_rsb_get_nzeros
    procedure, pass(a) :: get_fmt      => d_rsb_get_fmt
    procedure, pass(a) :: sizeof       => d_rsb_sizeof
    procedure, pass(a) :: d_csmm       => psb_d_rsb_csmm
    procedure, pass(a) :: d_csmv       => psb_d_rsb_csmv
    procedure, pass(a) :: d_inner_cssm => psb_d_rsb_cssm
    procedure, pass(a) :: d_inner_cssv => psb_d_rsb_cssv
    procedure, pass(a) :: d_scals      => psb_d_rsb_scals
    procedure, pass(a) :: d_scal       => psb_d_rsb_scal
    procedure, pass(a) :: csnmi        => psb_d_rsb_csnmi
    procedure, pass(a) :: csnm1        => psb_d_rsb_csnm1
    procedure, pass(a) :: rowsum       => psb_d_rsb_rowsum
    procedure, pass(a) :: arwsum       => psb_d_rsb_arwsum
    procedure, pass(a) :: colsum       => psb_d_rsb_colsum
    procedure, pass(a) :: aclsum       => psb_d_rsb_aclsum
!    procedure, pass(a) :: reallocate_nz => psb_d_rsb_reallocate_nz ! FIXME
!    procedure, pass(a) :: allocate_mnnz => psb_d_rsb_allocate_mnnz ! FIXME
    procedure, pass(a) :: cp_to_coo    => psb_d_cp_rsb_to_coo
    procedure, pass(a) :: cp_from_coo  => psb_d_cp_rsb_from_coo
    procedure, pass(a) :: cp_to_fmt    => psb_d_cp_rsb_to_fmt
    procedure, pass(a) :: cp_from_fmt  => psb_d_cp_rsb_from_fmt
    procedure, pass(a) :: mv_to_coo    => psb_d_mv_rsb_to_coo
    procedure, pass(a) :: mv_from_coo  => psb_d_mv_rsb_from_coo
    procedure, pass(a) :: mv_to_fmt    => psb_d_mv_rsb_to_fmt
    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_rsb_from_fmt
    procedure, pass(a) :: csput        => psb_d_rsb_csput
    procedure, pass(a) :: get_diag     => psb_d_rsb_get_diag
    procedure, pass(a) :: csgetptn     => psb_d_rsb_csgetptn
    procedure, pass(a) :: d_csgetrow   => psb_d_rsb_csgetrow
    procedure, pass(a) :: get_nz_row   => d_rsb_get_nz_row
    procedure, pass(a) :: reinit       => psb_d_rsb_reinit
    procedure, pass(a) :: trim         => psb_d_rsb_trim
    procedure, pass(a) :: print        => psb_d_rsb_print
    procedure, pass(a) :: free         => d_rsb_free
    procedure, pass(a) :: mold         => psb_d_rsb_mold
    procedure, pass(a) :: psb_d_rsb_cp_from
    generic, public    :: cp_from => psb_d_rsb_cp_from
    procedure, pass(a) :: psb_d_rsb_mv_from
    generic, public    :: mv_from => psb_d_rsb_mv_from

#endif
  end type
  ! FIXME: complete the following
  !private :: d_rsb_get_nzeros, d_rsb_get_fmt
  private :: d_rsb_to_psb_info
#ifdef HAVE_LIBRSB
  contains 

  function d_rsb_to_psb_info(info) result(res)
    implicit none 
    integer :: res,info
    res=info
  end function d_rsb_to_psb_info

  function d_rsb_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    res=rsb_get_matrix_nnz(a%rsbmptr)
  end function d_rsb_get_nzeros

  function d_rsb_get_fmt(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'RSB'
  end function d_rsb_get_fmt
  
  function d_rsb_get_size(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    res = d_rsb_get_nzeros(a)
  end function d_rsb_get_size

  function d_rsb_sizeof(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    res=rsb_sizeof(a%rsbmptr)
  end function d_rsb_sizeof

subroutine psb_d_rsb_csmv(alpha,a,x,beta,y,info,trans) 
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in) :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  character :: trans_
  info = psb_success_

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  info=d_rsb_to_psb_info(rsb_spmv(a%rsbmptr,x,y,alpha,beta,1,1,rsb_psblas_trans_to_rsb_trans(trans_)))
end subroutine psb_d_rsb_csmv

subroutine psb_d_rsb_cssv(alpha,a,x,beta,y,info,trans) 
  ! FIXME: and what when x is an alias of y ?
  ! FIXME: ignoring beta
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in) :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  character :: trans_
  info = psb_success_

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  info=d_rsb_to_psb_info(rsb_spsv(a%rsbmptr,x,y,alpha,1,1,rsb_psblas_trans_to_rsb_trans(trans_)))
end subroutine psb_d_rsb_cssv

subroutine psb_d_rsb_scals(d,a,info) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer, intent(out)            :: info
  info=d_rsb_to_psb_info(rsb_elemental_scale(a%rsbmptr,d))
end subroutine psb_d_rsb_scals

subroutine psb_d_rsb_scal(d,a,info) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer, intent(out)            :: info
  info=d_rsb_to_psb_info(rsb_scale_rows(a%rsbmptr,d))
end subroutine psb_d_rsb_scal

  subroutine  d_rsb_free(a) 
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(inout) :: a
    type(c_ptr) :: dummy
    dummy=rsb_free_sparse_matrix(a%rsbmptr)
  end subroutine d_rsb_free

subroutine  psb_d_rsb_trim(a)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  ! FIXME: this is supposed to remain empty for RSB
end subroutine  psb_d_rsb_trim

    subroutine psb_d_rsb_print(iout,a,iv,eirs,eics,head,ivr,ivc)
      integer, intent(in)               :: iout
      class(psb_d_rsb_sparse_mat), intent(in) :: a   
      integer, intent(in), optional     :: iv(:)
      integer, intent(in), optional     :: eirs,eics
      character(len=*), optional        :: head
      integer, intent(in), optional     :: ivr(:), ivc(:)
      integer             :: info
      ! FIXME: UNFINISHED
      info=rsb_print_matrix_t(a%rsbmptr)
    end subroutine psb_d_rsb_print

    subroutine psb_d_rsb_get_diag(a,d,info) 
      class(psb_d_rsb_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
      info=rsb_getdiag(a%rsbmptr,d)
    end subroutine psb_d_rsb_get_diag
  
function psb_d_rsb_csnmi(a) result(res)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res
  integer             :: info
  info=rsb_infinity_norm(a%rsbmptr,res,rsb_psblas_trans_to_rsb_trans('N'))
end function psb_d_rsb_csnmi

function psb_d_rsb_csnm1(a) result(res)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res
  integer             :: info
  info=rsb_one_norm(a%rsbmptr,res,rsb_psblas_trans_to_rsb_trans('N'))
end function psb_d_rsb_csnm1

subroutine psb_d_rsb_aclsum(d,a) 
  use psb_sparse_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)
  info=rsb_absolute_columns_sums(a%rsbmptr,d)
end subroutine psb_d_rsb_aclsum

subroutine psb_d_rsb_arwsum(d,a) 
  use psb_sparse_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)
  info=rsb_absolute_rows_sums(a%rsbmptr,d)
end subroutine psb_d_rsb_arwsum

subroutine psb_d_rsb_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer :: ldy,ldx,nc

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=d_rsb_to_psb_info(rsb_spmm(a%rsbmptr,x,y,ldx,ldy,nc,rsb_psblas_trans_to_rsb_trans(trans_),alpha,beta,c_f_order))
end subroutine psb_d_rsb_csmm

subroutine psb_d_rsb_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  integer :: ldy,ldx,nc
  character :: trans_
  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=d_rsb_to_psb_info(rsb_spsm(a%rsbmptr,y,ldy,nc,rsb_psblas_trans_to_rsb_trans(trans_),alpha,beta,c_f_order))
end subroutine

subroutine psb_d_rsb_rowsum(d,a) 
  use psb_sparse_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)
  integer :: info
  info=d_rsb_to_psb_info(rsb_rows_sums(a%rsbmptr,d))
end subroutine psb_d_rsb_rowsum

subroutine psb_d_rsb_colsum(d,a) 
  use psb_sparse_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)
  integer :: info
  info=d_rsb_to_psb_info(rsb_columns_sums(a%rsbmptr,d))
end subroutine psb_d_rsb_colsum

subroutine psb_d_rsb_mold(a,b,info) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer, intent(out)                    :: info
  Integer :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_d_rsb_sparse_mat :: b, stat=info)

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return
end subroutine psb_d_rsb_mold

subroutine psb_d_rsb_reinit(a,clear)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear
  Integer :: info
  info=d_rsb_to_psb_info(rsb_reinit(a%rsbmptr))
end subroutine psb_d_rsb_reinit


  function  d_rsb_get_nz_row(idx,a) result(res)
    implicit none
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    integer                              :: info
    res=0
    res=rsb_get_rows_nnz(a%rsbmptr,idx,idx,c_f_index,info)
    info=d_rsb_to_psb_info(info)
    if(info.ne.0)res=0
  end function d_rsb_get_nz_row

subroutine psb_d_cp_rsb_to_coo(a,b,info) 
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                      :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, idl,err_act
  integer             :: debug_level, debug_unit
  character(len=20)   :: name
  info = psb_success_
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()
  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)
  info=d_rsb_to_psb_info(rsb_get_coo(a%rsbmptr,b%val,b%ia,b%ja,c_f_index))
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)
end subroutine psb_d_cp_rsb_to_coo

subroutine psb_d_cp_rsb_to_fmt(a,b,info) 
  use psb_sparse_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(in)   :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                       :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_d_rsb_sparse_mat) 
    call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)! FIXME: ?
    b%rsbmptr=rsb_clone(a%rsbmptr) ! FIXME is thi enough ?
    ! FIXME: error handling needed here

  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_rsb_to_fmt

subroutine psb_d_cp_rsb_from_coo(a,b,info) 
  use psb_sparse_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

   info = psb_success_
   ! This is to have fix_coo called behind the scenes
   call b%cp_to_coo(tmp,info)
   if (info == psb_success_) call a%mv_from_coo(tmp,info)
   a%rsbmptr=rsb_allocate_rsb_sparse_matrix_const&
  &(tmp%val,tmp%ia,tmp%ja,tmp%get_nzeros(),c_d_typecode,tmp%get_nrows(),tmp%get_ncols(),1,1,c_def_flags,info)
    info=d_rsb_to_psb_info(info)
   ! FIXME: should destroy tmp ?
end subroutine psb_d_cp_rsb_from_coo

subroutine psb_d_cp_rsb_from_fmt(a,b,info) 
  use psb_sparse_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nz, nr, i,j,irw, idl,err_act, nc
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_rsb_sparse_mat) 
    call b%cp_to_fmt(a,info)  ! FIXME
    ! FIXME: missing error handling 

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_rsb_from_fmt


subroutine psb_d_rsb_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_sparse_mod
  implicit none

  class(psb_d_rsb_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  ! FIXME: MISSING THE HANDLING OF OPTIONS, HERE

  call psb_erractionsave(err_act)
  info = psb_success_
  
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif

  if ((imax<imin).or.(jmax_<jmin_)) then 
    nz = 0
    return
  end if

  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif
  if ((append_).and.(present(nzin))) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .false.
  endif
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .false.
  endif
  if ((rscale_.or.cscale_).and.(present(iren))) then 
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='iren (rscale.or.cscale)')
    goto 9999
  end if

    nzin_=rsb_get_block_nnz(a%rsbmptr,imin,imax,imin,imax,c_f_index,info)
    ! FIXME: unfinished; missing error handling ..

    call psb_ensure_size(nzin_,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_,ja,info)
    if (info == psb_success_) call psb_ensure_size(nzin_,val,info)
    if (info /= psb_success_) return
    
    info=d_rsb_to_psb_info(rsb_get_block_sparse&
    &(a%rsbmptr,val,imin,imax,jmin,jmax,ia,ja,c_null_ptr,c_null_ptr,nzin_,c_f_index))
    ! FIXME: unfinished; missing error handling ..

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if

  end subroutine psb_d_rsb_csgetrow

subroutine psb_d_rsb_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_sparse_mod
  implicit none

  class(psb_d_rsb_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    !nzt = ..
    nz = 0 

    call psb_ensure_size(nzin_,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_,ja,info)

    if (info /= psb_success_) return
    nzin_=rsb_get_block_nnz(a%rsbmptr,imin,imax,imin,imax,c_f_index,info)
    ! FIXME: unfinished; missing error handling ..

    call psb_ensure_size(nzin_,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_,ja,info)
    if (info /= psb_success_) return
    
    info=d_rsb_to_psb_info(rsb_get_block_sparse_pattern&
    &(a%rsbmptr,imin,imax,jmin,jmax,ia,ja,c_null_ptr,c_null_ptr,nzin_,c_f_index))
    ! FIXME: unfinished; missing error handling ..
end subroutine psb_d_rsb_csgetptn

subroutine psb_d_rsb_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_sparse_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_rsb_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)
  info=d_rsb_to_psb_info(rsb_update_elements(a%rsbmptr,val,ia,ja,nz,c_upd_flags))
end subroutine psb_d_rsb_csput

subroutine psb_d_mv_rsb_to_coo(a,b,info) 
  use psb_sparse_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer, intent(out)                        :: info
  call psb_d_cp_rsb_to_coo(a,b,info)
  call d_rsb_free(a) 
end subroutine psb_d_mv_rsb_to_coo

    subroutine psb_d_mv_rsb_to_fmt(a,b,info) 
      class(psb_d_rsb_sparse_mat), intent(inout) :: a
      class(psb_d_base_sparse_mat), intent(inout)  :: b
      integer, intent(out)                        :: info
  call psb_d_cp_rsb_to_fmt(a,b,info)
  call d_rsb_free(a)
    end subroutine psb_d_mv_rsb_to_fmt
  
subroutine psb_d_mv_rsb_from_fmt(a,b,info) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info
  type(psb_d_coo_sparse_mat) :: tmp
  info = psb_success_
  select type (b)
  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_mv_rsb_from_fmt

subroutine psb_d_mv_rsb_from_coo(a,b,info) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info
  call a%cp_from_coo(b,info)
  call b%free()
end subroutine psb_d_mv_rsb_from_coo

subroutine psb_d_rsb_cp_from(a,b)
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  type(psb_d_rsb_sparse_mat), intent(in)   :: b
  Integer :: info
  type(psb_d_coo_sparse_mat) :: tmp
  call b%cp_to_coo(tmp,info)
  call a%mv_from_coo(tmp,info)
  call tmp%free()
end subroutine psb_d_rsb_cp_from

subroutine psb_d_rsb_mv_from(a,b)
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  type(psb_d_rsb_sparse_mat), intent(inout)   :: b
  Integer :: info
  type(psb_d_coo_sparse_mat) :: tmp
  call b%mv_to_coo(tmp,info)
  call a%mv_from_coo(tmp,info)
end subroutine psb_d_rsb_mv_from


#endif
end module psb_d_rsb_mat_mod
