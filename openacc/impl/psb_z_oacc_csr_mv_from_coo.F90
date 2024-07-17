submodule (psb_z_oacc_csr_mat_mod) psb_z_oacc_csr_mv_from_coo_impl
  use psb_base_mod
contains
  module subroutine psb_z_oacc_csr_mv_from_coo(a, b, info)
    implicit none 

    class(psb_z_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_z_coo_sparse_mat), intent(inout)  :: b
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_

    call a%psb_z_csr_sparse_mat%mv_from_coo(b, info)
    if (info /= 0) goto 9999

    !$acc update device(a%val, a%ja, a%irp)

    return

9999 continue
    info = psb_err_alloc_dealloc_
    return

  end subroutine psb_z_oacc_csr_mv_from_coo
end submodule psb_z_oacc_csr_mv_from_coo_impl
