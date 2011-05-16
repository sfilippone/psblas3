!
!
! FIXME/TODO:
! * some RSB constants are used in their value form, and with no explanation
! * error handling
! * PSBLAS interface adherence
! * should test and fix all the problems that for sure will occur
! * duplicate handling is not defined
! * the printing function is not complete
! * should substitute -1 with another valid PSBLAS error code
! * ..
! 
module psb_d_rsb_mat_mod
  use psb_d_base_mat_mod
  use rsb_d_mod
#ifdef HAVE_LIBRSB
  use iso_c_binding  
#endif
#if 0
#define PSBRSB_DEBUG(MSG) write(*,*) __FILE__,':',__LINE__,':',MSG
#define PSBRSB_ERROR(MSG) write(*,*) __FILE__,':',__LINE__,':'," ERROR: ",MSG
#define PSBRSB_WARNING(MSG) write(*,*) __FILE__,':',__LINE__,':'," WARNING: ",MSG
#else
#define PSBRSB_DEBUG(MSG)
#define PSBRSB_ERROR(MSG)
#define PSBRSB_WARNING(MSG)
#endif
  integer, parameter :: c_typecode=68 ! this is module specific 
  integer, parameter :: c_for_flags=1 ! : here should use RSB_FLAG_FORTRAN_INDICES_INTERFACE
  integer, parameter :: c_srt_flags =4 ! flags if rsb input is row major sorted ..
  !integer, parameter :: c_own_flags =-1 ! flags if rsb input shall not be freed by rsb
  integer, parameter :: c_tri_flags =8 ! flags for specifying a triangle
  integer, parameter :: c_low_flags =16 ! flags for specifying a lower triangle/symmetry
  integer, parameter :: c_upp_flags =32 ! flags for specifying a lower triangle/symmetry
  integer, parameter :: c_idi_flags =64 ! flags for specifying diagonal implicit
  integer, parameter :: c_def_flags =c_for_flags ! FIXME: here should use ..
  integer :: c_f_order=c_for_flags        ! FIXME: here should use RSB_FLAG_WANT_COLUMN_MAJOR_ORDER
  integer, parameter :: c_upd_flags =c_for_flags ! flags for when updating the assembled rsb matrix
  integer, parameter :: c_psbrsb_err_ =psb_err_internal_error_
  type, extends(psb_d_base_sparse_mat) :: psb_d_rsb_sparse_mat
#ifdef HAVE_LIBRSB
    type(c_ptr) :: rsbmptr=c_null_ptr
  contains 
    procedure, pass(a) :: get_size     => d_rsb_get_size
    procedure, pass(a) :: get_nzeros   => d_rsb_get_nzeros
    procedure, pass(a) :: get_ncols   => d_rsb_get_ncols
    procedure, pass(a) :: get_nrows   => d_rsb_get_nrows
    procedure, nopass  :: get_fmt      => d_rsb_get_fmt
    procedure, pass(a) :: sizeof       => d_rsb_sizeof
    procedure, pass(a) :: d_csmm       => psb_d_rsb_csmm
    !procedure, pass(a) :: d_csmv_nt       => psb_d_rsb_csmv_nt ! FIXME: a placeholder for future memory
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
    procedure, pass(a) :: trim         => psb_d_rsb_trim ! evil
    procedure, pass(a) :: print        => psb_d_rsb_print
    procedure, pass(a) :: free         => d_rsb_free
    procedure, pass(a) :: mold         => psb_d_rsb_mold
    procedure, pass(a) :: psb_d_rsb_cp_from
    generic, public    :: cp_from => psb_d_rsb_cp_from
    procedure, pass(a) :: psb_d_rsb_mv_from
    generic, public    :: mv_from => psb_d_rsb_mv_from

#endif
  end type psb_d_rsb_sparse_mat
  ! FIXME: complete the following
  !private :: d_rsb_get_nzeros, d_rsb_get_fmt
  private :: d_rsb_to_psb_info
#ifdef HAVE_LIBRSB
  contains 

  function psb_rsb_matmod_init() result(res)
    implicit none 
    integer :: res
    !PSBRSB_DEBUG('')
    res=-1 ! FIXME
#ifdef HAVE_LIBRSB
    res=d_rsb_to_psb_info(rsb_init(c_null_ptr))
#endif
  end function psb_rsb_matmod_init

  function psb_rsb_matmod_exit() result(res)
    implicit none 
    integer :: res
    !PSBRSB_DEBUG('')
    res=-1 ! FIXME
#ifdef HAVE_LIBRSB
    res=d_rsb_to_psb_info(rsb_exit())
#endif
  end function psb_rsb_matmod_exit

 function d_rsb_to_psb_info(info) result(res)
    implicit none 
    integer , intent(in) :: info
    integer :: res
    !PSBRSB_DEBUG('')
    if(info.ne.0)then
      res=-1
    else
      res=psb_success_
    end if
  end function d_rsb_to_psb_info

 function d_rsb_get_flags(a) result(flags)
    implicit none 
    integer :: flags
    class(psb_d_base_sparse_mat), intent(in) :: a
    !PSBRSB_DEBUG('')
    flags=c_def_flags 
    if(a%is_sorted()) flags=flags+c_srt_flags 
    if(a%is_triangle()) flags=flags+c_tri_flags 
    if(a%is_upper()) flags=flags+c_upp_flags 
    if(a%is_unit()) flags=flags+c_idi_flags 
    if(a%is_lower()) flags=flags+c_low_flags 
  end function d_rsb_get_flags

  function d_rsb_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    !PSBRSB_DEBUG('')
    res=rsb_get_matrix_nnz(a%rsbmptr)
  end function d_rsb_get_nzeros

  function d_rsb_get_nrows(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    !PSBRSB_DEBUG('')
    res=rsb_get_matrix_n_rows(a%rsbmptr)
  end function d_rsb_get_nrows

  function d_rsb_get_ncols(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    !PSBRSB_DEBUG('')
    res=rsb_get_matrix_n_columns(a%rsbmptr)
  end function d_rsb_get_ncols

  function d_rsb_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    !the following printout is harmful, here, if happening during a write :) (causes a deadlock)
    !PSBRSB_DEBUG('')
    res = 'RSB'
  end function d_rsb_get_fmt
  
  function d_rsb_get_size(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    !PSBRSB_DEBUG('')
    res = d_rsb_get_nzeros(a)
  end function d_rsb_get_size

  function d_rsb_sizeof(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer(psb_long_int_k_) :: res
    !PSBRSB_DEBUG('')
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
!    PSBRSB_DEBUG('')
  info = psb_success_

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  info=d_rsb_to_psb_info(rsb_spmv(rsb_psblas_trans_to_rsb_trans(trans_),alpha,a%rsbmptr,x,1,beta,y,1))
end subroutine psb_d_rsb_csmv

subroutine psb_d_rsb_csmv_nt(alpha,a,x1,x2,beta,y1,y2,info) 
  ! FIXME: this routine is here as a placeholder for a specialized implementation of 
  ! joint spmv and spmv transposed.
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in) :: alpha, beta, x1(:), x2(:)
  real(psb_dpk_), intent(inout)       :: y1(:), y2(:)
  integer, intent(out)                :: info
!    PSBRSB_DEBUG('')
  info = psb_success_
  info=d_rsb_to_psb_info(rsb_spmv_nt(alpha,a%rsbmptr,x1,x2,1,beta,y1,y2,1))
  return
end subroutine psb_d_rsb_csmv_nt

subroutine psb_d_rsb_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  ! FIXME: and what when x is an alias of y ?
  ! FIXME: ignoring beta
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in) :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  character :: trans_
  Integer :: err_act, i
  character(len=20)  :: name='rsb_cssv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

!    PSBRSB_DEBUG('')

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  info=d_rsb_to_psb_info(rsb_spsv(rsb_psblas_trans_to_rsb_trans(trans_),alpha,a%rsbmptr,x,1,y,1))
  if (info /= 0) then 
    i = info
    info = psb_err_from_subroutine_ai_
    call psb_errpush(info,name,&
         & i_err=(/i,0,0,0,0/),a_err="rsb_spsv")
    goto 9999
  end if
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    PSBRSB_ERROR("!")
    call psb_error()
    return
  end if
  return
    
end subroutine psb_d_rsb_cssv

subroutine psb_d_rsb_scals(d,a,info) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer, intent(out)            :: info
    PSBRSB_DEBUG('')
  info=d_rsb_to_psb_info(rsb_elemental_scale(a%rsbmptr,d))
end subroutine psb_d_rsb_scals

subroutine psb_d_rsb_scal(d,a,info) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer, intent(out)            :: info
    PSBRSB_DEBUG('')
  info=d_rsb_to_psb_info(rsb_scale_rows(a%rsbmptr,d))
end subroutine psb_d_rsb_scal

  subroutine  d_rsb_free(a) 
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(inout) :: a
    type(c_ptr) :: dummy
    !PSBRSB_DEBUG('freeing RSB matrix')
    dummy=rsb_free_sparse_matrix(a%rsbmptr)
  end subroutine d_rsb_free

subroutine  psb_d_rsb_trim(a)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
    !PSBRSB_DEBUG('')
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
    PSBRSB_DEBUG('')
      ! FIXME: UNFINISHED
      info=rsb_print_matrix_t(a%rsbmptr)
    end subroutine psb_d_rsb_print

    subroutine psb_d_rsb_get_diag(a,d,info) 
      class(psb_d_rsb_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(out)     :: d(:)
      integer, intent(out)            :: info
      !PSBRSB_DEBUG('')
      info=rsb_getdiag(a%rsbmptr,d)
    end subroutine psb_d_rsb_get_diag
  
function psb_d_rsb_csnmi(a) result(csnmi_res)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_),target         :: csnmi_res ! please DO NOT rename this variable (see the Makefile)
  real(psb_dpk_)         :: resa(1)
  integer             :: info
    !PSBRSB_DEBUG('')
  info=rsb_infinity_norm(a%rsbmptr,resa,rsb_psblas_trans_to_rsb_trans('N'))
  !info=rsb_infinity_norm(a%rsbmptr,c_loc(res),rsb_psblas_trans_to_rsb_trans('N'))
  csnmi_res=resa(1)
end function psb_d_rsb_csnmi

function psb_d_rsb_csnm1(a) result(res)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res
  real(psb_dpk_)         :: resa(1)
  integer             :: info
    PSBRSB_DEBUG('')
  info=rsb_one_norm(a%rsbmptr,resa,rsb_psblas_trans_to_rsb_trans('N'))
  !info=rsb_one_norm(a%rsbmptr,res,rsb_psblas_trans_to_rsb_trans('N'))
end function psb_d_rsb_csnm1

subroutine psb_d_rsb_aclsum(d,a) 
  use psb_base_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)
    PSBRSB_DEBUG('')
  info=rsb_absolute_columns_sums(a%rsbmptr,d)
end subroutine psb_d_rsb_aclsum

subroutine psb_d_rsb_arwsum(d,a) 
  use psb_base_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)
    PSBRSB_DEBUG('')
  info=rsb_absolute_rows_sums(a%rsbmptr,d)
end subroutine psb_d_rsb_arwsum

subroutine psb_d_rsb_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer :: ldy,ldx,nc
    PSBRSB_DEBUG('')
    PSBRSB_DEBUG('ERROR: UNIMPLEMENTED')

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=-1
  info=d_rsb_to_psb_info(rsb_spmm(rsb_psblas_trans_to_rsb_trans(trans_),alpha,a%rsbmptr,nc,c_f_order,x,ldx,beta,y,ldy))
end subroutine psb_d_rsb_csmm

subroutine psb_d_rsb_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  integer :: ldy,ldx,nc
  character :: trans_
    PSBRSB_DEBUG('')
    PSBRSB_DEBUG('ERROR: UNIMPLEMENTED')
  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=-1
  info=d_rsb_to_psb_info(rsb_spsm(rsb_psblas_trans_to_rsb_trans(trans_),alpha,a%rsbmptr,nc,c_f_order,beta,x,ldx,y,ldy))
end subroutine

subroutine psb_d_rsb_rowsum(d,a) 
  use psb_base_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)
  integer :: info
    PSBRSB_DEBUG('')
  info=d_rsb_to_psb_info(rsb_rows_sums(a%rsbmptr,d))
end subroutine psb_d_rsb_rowsum

subroutine psb_d_rsb_colsum(d,a) 
  use psb_base_mod
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)
  integer :: info
    PSBRSB_DEBUG('')
  info=d_rsb_to_psb_info(rsb_columns_sums(a%rsbmptr,d))
end subroutine psb_d_rsb_colsum

subroutine psb_d_rsb_mold(a,b,info) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer, intent(out)                    :: info
  Integer :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.
    PSBRSB_DEBUG('')

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
    PSBRSB_ERROR("!")
    call psb_error()
  end if
  return
end subroutine psb_d_rsb_mold

subroutine psb_d_rsb_reinit(a,clear)
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear
  Integer :: info
    PSBRSB_DEBUG('')
  info=d_rsb_to_psb_info(rsb_reinit_matrix(a%rsbmptr))
end subroutine psb_d_rsb_reinit


  function  d_rsb_get_nz_row(idx,a) result(res)
    implicit none
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer, intent(in)                  :: idx
    integer                              :: res
    integer                              :: info
    PSBRSB_DEBUG('')
    res=0
    res=rsb_get_rows_nnz(a%rsbmptr,idx,idx,c_for_flags,info)
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
   ! PSBRSB_DEBUG('')
  info = psb_success_
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()
  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)
  info=d_rsb_to_psb_info(rsb_get_coo(a%rsbmptr,b%val,b%ia,b%ja,c_for_flags))
  call b%set_nzeros(a%get_nzeros())
  call b%set_nrows(a%get_nrows())
  call b%set_ncols(a%get_ncols())
  call b%fix(info)
  !write(*,*)b%val
  !write(*,*)b%ia
  !write(*,*)b%ja
  !write(*,*)b%get_nrows()
  !write(*,*)b%get_ncols()
  !write(*,*)b%get_nzeros()
  !write(*,*)a%get_nrows()
  !write(*,*)a%get_ncols()
  !write(*,*)a%get_nzeros()
end subroutine psb_d_cp_rsb_to_coo

subroutine psb_d_cp_rsb_to_fmt(a,b,info) 
  use psb_base_mod
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
    PSBRSB_DEBUG('')

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
  use psb_base_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  integer             :: debug_level, debug_unit
  integer             :: flags
  character(len=20)   :: name
   ! PSBRSB_DEBUG('')

   flags=d_rsb_get_flags(b)
  
   info = psb_success_
   call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)

   !write (*,*) b%val
   ! FIXME: and if sorted ? the process could be speeded up !
   a%rsbmptr=rsb_allocate_rsb_sparse_matrix_const&
  &(b%val,b%ia,b%ja,b%get_nzeros(),c_typecode,b%get_nrows(),b%get_ncols(),1,1,flags,info)
    info=d_rsb_to_psb_info(info)
   ! FIXME: should destroy tmp ?
end subroutine psb_d_cp_rsb_from_coo

subroutine psb_d_cp_rsb_from_fmt(a,b,info) 
  use psb_base_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nz, nr, i,j,irw, idl,err_act, nc
  integer              :: debug_level, debug_unit
  integer             :: flags
  character(len=20)   :: name
    PSBRSB_DEBUG('')

  info = psb_success_
  flags=d_rsb_get_flags(b)

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_csr_sparse_mat) 
    call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
    a%rsbmptr=rsb_allocate_rsb_sparse_matrix_from_csr_const&
        &(b%val,b%irp,b%ja,b%get_nzeros(),c_typecode,b%get_nrows(),b%get_ncols(),1,1,flags,info)
    info=d_rsb_to_psb_info(info)

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
  use psb_base_mod
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
  integer :: nzin_, jmin_, jmax_, err_act, i, nzrsb 
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  ! FIXME: MISSING THE HANDLING OF OPTIONS, HERE
  PSBRSB_DEBUG('')

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
    !info=c_psbrsb_err_ 
    PSBRSB_WARNING("imax < imin ? or jmax < jmin ? !")
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
!  if ((rscale_.or.cscale_).and.(present(iren))) then 
!    PSBRSB_ERROR("!")
!    info = psb_err_many_optional_arg_
!    call psb_errpush(info,name,a_err='iren (rscale.or.cscale)')
!    goto 9999
!  end if

  nzrsb = rsb_get_block_nnz(a%rsbmptr,imin,imax,jmin_,jmax_,c_for_flags,info)
  ! FIXME: unfinished; missing error handling ..
  
  call psb_ensure_size(nzin_+nzrsb,ia,info)
  if (info == psb_success_) call psb_ensure_size(nzin_+nzrsb,ja,info)
  if (info == psb_success_) call psb_ensure_size(nzin_+nzrsb,val,info)
  if (info /= psb_success_)then
    PSBRSB_ERROR("psb_ensure_size failed !")
    return
  endif
  

  info=d_rsb_to_psb_info(rsb_get_block_sparse(a%rsbmptr,&
       & val(nzin_+1:),imin,imax,jmin_,jmax_,&
       & ia(nzin_+1:),ja(nzin_+1:),&
       & c_null_ptr,c_null_ptr,nz,c_for_flags))
    ! FIXME: unfinished; missing error handling ..
  if (nz /= nzrsb) then 
    info=c_psbrsb_err_ 
    PSBRSB_ERROR("Mismatch in output from rsb_getblk")
    !write(*,*) 'Mismatch in output from rsb_getblk: ',nz,nzrsb
 end if

  if (rscale_) then 
    do i=nzin_+1, nzin_+nz
      ia(i) = ia(i) - imin + 1
    end do
  end if
  if (cscale_) then 
    do i=nzin_+1, nzin_+nz
      ja(i) = ja(i) - jmin_ + 1
    end do
  end if



9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    PSBRSB_ERROR("!")
    call psb_error()
    return
  end if

end subroutine psb_d_rsb_csgetrow

subroutine psb_d_rsb_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_base_mod
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
  PSBRSB_DEBUG('')

  if (present(iren).or.present(rscale).or.present(cscale)) then 
    ! FIXME: error condition
    PSBRSB_ERROR("unsupported optional arguments!")
    call psb_error()
  endif

  if (present(append)) then 
    append_ = append
  else
    append_ = .false.
  endif
  if (present(append).and.append.and.present(nzin)) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif

  if (present(jmin)) then 
    jmin_ = jmin
  else
    jmin_ = 1
  endif

  if (present(jmax)) then 
    jmax_ = jmax
  else
    jmax_ = a%get_nrows()
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
  if (present(iren)) then 
    info = c_psbrsb_err_ 
    PSBRSB_ERROR("ERROR: the RSB pattern get needs iren support !!")
    goto 9999
  end if


  !nzt = ..
  nz = 0 

  call psb_ensure_size(nzin_,ia,info)
  if (info == psb_success_) call psb_ensure_size(nzin_,ja,info)

  if (info /= psb_success_) return
  nz=rsb_get_block_nnz(a%rsbmptr,imin,imax,jmin_,jmax_,c_for_flags,info)
  !write(*,*) 'debug:',nzin_,nz,imin,imax,jmin_,jmax_
  ! FIXME: unfinished; missing error handling ..

  call psb_ensure_size(nzin_+nz,ia,info)
  if (info == psb_success_) call psb_ensure_size(nzin_+nz,ja,info)
  if (info /= psb_success_)then
    PSBRSB_ERROR("!")
    return
  endif

  info=d_rsb_to_psb_info(rsb_get_block_sparse_pattern&
       &(a%rsbmptr,imin,imax,jmin_,jmax_,ia,ja,c_null_ptr,c_null_ptr,nzin_,c_for_flags))
  ! FIXME: unfinished; missing error handling ..

  !write(*,*) 'debug:',nzin_,nz,imin,imax,jmin_,jmax_
  if (rscale_) then 
    do i=nzin_+1, nzin_+nz
      ia(i) = ia(i) - imin + 1
    end do
  end if
  if (cscale_) then 
    do i=nzin_+1, nzin_+nz
      ja(i) = ja(i) - jmin_ + 1
    end do
  end if

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    PSBRSB_ERROR("!")
    call psb_error()
    return
  endif

end subroutine psb_d_rsb_csgetptn

subroutine psb_d_rsb_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_base_mod
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
  PSBRSB_DEBUG('')
  if(present(gtl))then
    PSBRSB_ERROR("!")
  endif
  info=d_rsb_to_psb_info(rsb_update_elements(a%rsbmptr,val,ia,ja,nz,c_upd_flags))
end subroutine psb_d_rsb_csput

subroutine psb_d_mv_rsb_to_coo(a,b,info) 
  use psb_base_mod
  implicit none 

  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer, intent(out)                        :: info
    PSBRSB_DEBUG('')
  ! FIXME: use rsb_switch_rsb_matrix_to_coo_sorted !
  call psb_d_cp_rsb_to_coo(a,b,info)
  call a%free()
end subroutine psb_d_mv_rsb_to_coo

subroutine psb_d_mv_rsb_to_fmt(a,b,info) 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout)  :: b
  integer, intent(out)                        :: info
  PSBRSB_DEBUG('')
  ! FIXME: could use here rsb_switch_rsb_matrix_to_csr_sorted
  call psb_d_cp_rsb_to_fmt(a,b,info)
  call d_rsb_free(a)
  a%rsbmptr=c_null_ptr
end subroutine psb_d_mv_rsb_to_fmt
  
subroutine psb_d_mv_rsb_from_fmt(a,b,info) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info
  ! FIXME: could use here rsb_allocate_rsb_sparse_matrix_from_csr_inplace
  !if(b%is_sorted()) flags=flags+c_srt_flags 
  type(psb_d_coo_sparse_mat) :: tmp
  !  PSBRSB_DEBUG('')
  info = psb_success_
  select type (b)
  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_mv_rsb_from_fmt

subroutine psb_d_mv_rsb_from_coo(a,b,info) 
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info
  !  PSBRSB_DEBUG('')
  ! FIXME: should use rsb_allocate_rsb_sparse_matrix_inplace
  !if(b%is_sorted()) flags=flags+c_srt_flags 
  !if(b%is_triangle()) flags=flags+c_tri_flags 
  call a%cp_from_coo(b,info)
  call b%free()
end subroutine psb_d_mv_rsb_from_coo

subroutine psb_d_rsb_cp_from(a,b)
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  type(psb_d_rsb_sparse_mat), intent(in)   :: b
  Integer :: info
  type(psb_d_coo_sparse_mat) :: tmp
    PSBRSB_DEBUG('')
  call b%cp_to_coo(tmp,info)
  call a%mv_from_coo(tmp,info)
  call tmp%free()
end subroutine psb_d_rsb_cp_from

subroutine psb_d_rsb_mv_from(a,b)
  use psb_base_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(inout) :: a
  type(psb_d_rsb_sparse_mat), intent(inout)   :: b
  Integer :: info
  type(psb_d_coo_sparse_mat) :: tmp
    PSBRSB_DEBUG('')
  call b%mv_to_coo(tmp,info)
  call a%mv_from_coo(tmp,info)
end subroutine psb_d_rsb_mv_from


#endif
end module psb_d_rsb_mat_mod
