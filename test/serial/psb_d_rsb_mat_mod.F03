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
    procedure, pass(a) :: get_nzeros   => d_rsb_get_nzeros
#endif
  end type
  private :: d_rsb_get_nzeros
#ifdef HAVE_LIBRSB
  contains 

  function d_rsb_get_nzeros(a) result(res)
    implicit none 
    class(psb_d_rsb_sparse_mat), intent(in) :: a
    integer :: res
    res = rsb_get_matrix_nnz(a%rsbmptr)
  end function d_rsb_get_nzeros

#endif
end module psb_d_rsb_mat_mod
