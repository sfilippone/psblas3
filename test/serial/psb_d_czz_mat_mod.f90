module psb_d_czz_mat_mod

  use psb_d_base_mat_mod

  type, extends(psb_d_base_sparse_mat) :: psb_d_czz_sparse_mat

    integer, allocatable :: irp(:), ja(:)
    real(psb_dpk_), allocatable :: val(:)

  contains
  end type psb_d_czz_sparse_mat
contains 
end module psb_d_czz_mat_mod
