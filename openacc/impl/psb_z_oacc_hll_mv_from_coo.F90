submodule (psb_z_oacc_hll_mat_mod) psb_z_oacc_hll_mv_from_coo_impl
  use psb_base_mod
contains
  module subroutine psb_z_oacc_hll_mv_from_coo(a, b, info)
    implicit none 

    class(psb_z_oacc_hll_sparse_mat), intent(inout) :: a
    class(psb_z_coo_sparse_mat), intent(inout)      :: b
    integer(psb_ipk_), intent(out)                  :: info

    info = psb_success_

    call a%free_dev_space()
    call a%psb_z_hll_sparse_mat%mv_from_coo(b, info)
    if (info /= 0) goto 9999
    call a%sync_dev_space()
    call a%set_host()
    call a%sync()

    return

9999 continue
    info = psb_err_alloc_dealloc_
    return

  end subroutine psb_z_oacc_hll_mv_from_coo
end submodule psb_z_oacc_hll_mv_from_coo_impl
