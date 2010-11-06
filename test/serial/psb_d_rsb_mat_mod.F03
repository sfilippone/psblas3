module psb_d_rsb_mat_mod
  use psb_d_base_mat_mod
  use rsb_mod
#ifdef HAVE_LIBRSB
    use iso_c_binding  
#endif
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
!    procedure, pass(a) :: cp_to_coo    => psb_d_cp_rsb_to_coo
!    procedure, pass(a) :: cp_from_coo  => psb_d_cp_rsb_from_coo
!    procedure, pass(a) :: cp_to_fmt    => psb_d_cp_rsb_to_fmt
!    procedure, pass(a) :: cp_from_fmt  => psb_d_cp_rsb_from_fmt
!    procedure, pass(a) :: mv_to_coo    => psb_d_mv_rsb_to_coo
!    procedure, pass(a) :: mv_from_coo  => psb_d_mv_rsb_from_coo
!    procedure, pass(a) :: mv_to_fmt    => psb_d_mv_rsb_to_fmt
!    procedure, pass(a) :: mv_from_fmt  => psb_d_mv_rsb_from_fmt
!    procedure, pass(a) :: csput        => psb_d_rsb_csput
    procedure, pass(a) :: get_diag     => psb_d_rsb_get_diag
!    procedure, pass(a) :: csgetptn     => psb_d_rsb_csgetptn
!    procedure, pass(a) :: d_csgetrow   => psb_d_rsb_csgetrow
    procedure, pass(a) :: get_nz_row   => d_rsb_get_nz_row
    procedure, pass(a) :: reinit       => psb_d_rsb_reinit
    procedure, pass(a) :: trim         => psb_d_rsb_trim
    procedure, pass(a) :: print        => psb_d_rsb_print
    procedure, pass(a) :: free         => d_rsb_free
    procedure, pass(a) :: mold         => psb_d_rsb_mold
!    procedure, pass(a) :: psb_d_rsb_cp_from
!    generic, public    :: cp_from => psb_d_rsb_cp_from
!    procedure, pass(a) :: psb_d_rsb_mv_from
!    generic, public    :: mv_from => psb_d_rsb_mv_from

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
  integer :: order=2 ! FIXME

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=d_rsb_to_psb_info(rsb_spmm(a%rsbmptr,x,y,ldx,ldy,nc,rsb_psblas_trans_to_rsb_trans(trans_),alpha,beta,order))
end subroutine psb_d_rsb_csmm

subroutine psb_d_rsb_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_sparse_mod
  implicit none 
  class(psb_d_rsb_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans
  integer :: order=2 ! FIXME
  integer :: ldy,ldx,nc
  character :: trans_
  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  ldx=size(x,1); ldy=size(y,1)
  nc=min(size(x,2),size(y,2) )
  info=d_rsb_to_psb_info(rsb_spsm(a%rsbmptr,y,ldy,nc,rsb_psblas_trans_to_rsb_trans(trans_),alpha,beta,order))
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
    res=rsb_get_rows_nnz(a%rsbmptr,idx-1,idx-1,info)
    info=d_rsb_to_psb_info(info)
    if(info.ne.0.0)res=0
  end function d_rsb_get_nz_row


#endif
end module psb_d_rsb_mat_mod
