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
!    procedure, pass(a) :: d_csmm       => psb_d_rsb_csmm
    procedure, pass(a) :: d_csmv       => psb_d_rsb_csmv
!    procedure, pass(a) :: d_inner_cssm => psb_d_rsb_cssm
!    procedure, pass(a) :: d_inner_cssv => psb_d_rsb_cssv
!    procedure, pass(a) :: d_scals      => psb_d_rsb_scals
!    procedure, pass(a) :: d_scal       => psb_d_rsb_scal
!    procedure, pass(a) :: csnmi        => psb_d_rsb_csnmi
!    procedure, pass(a) :: csnm1        => psb_d_rsb_csnm1
!    procedure, pass(a) :: rowsum       => psb_d_rsb_rowsum
!    procedure, pass(a) :: arwsum       => psb_d_rsb_arwsum
!    procedure, pass(a) :: colsum       => psb_d_rsb_colsum
!    procedure, pass(a) :: aclsum       => psb_d_rsb_aclsum
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
!    procedure, pass(a) :: get_diag     => psb_d_rsb_get_diag
!    procedure, pass(a) :: csgetptn     => psb_d_rsb_csgetptn
!    procedure, pass(a) :: d_csgetrow   => psb_d_rsb_csgetrow
!    procedure, pass(a) :: get_nz_row   => d_rsb_get_nz_row
!    procedure, pass(a) :: reinit       => psb_d_rsb_reinit
!    procedure, pass(a) :: trim         => psb_d_rsb_trim
!    procedure, pass(a) :: print        => psb_d_rsb_print
!    procedure, pass(a) :: free         => d_rsb_free
!    procedure, pass(a) :: mold         => psb_d_rsb_mold
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

  function d_psb_to_rsb_trans(trans) result(res)
    implicit none 
    character :: trans
    integer :: res
    res=0 !FIXME
  end function d_psb_to_rsb_trans

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
  info=d_rsb_to_psb_info(rsb_spmv(a%rsbmptr,x,y,alpha,beta,1,1,d_psb_to_rsb_trans(trans_)))
end subroutine psb_d_rsb_csmv

#endif
end module psb_d_rsb_mat_mod
