submodule (psb_s_oacc_csr_mat_mod) psb_s_oacc_csr_cp_from_coo_impl
  use psb_base_mod
contains
  module subroutine psb_s_oacc_csr_cp_from_coo(a, b, info)
    implicit none

    class(psb_s_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_s_coo_sparse_mat), intent(in)         :: b
    integer(psb_ipk_), intent(out)                  :: info

    info = psb_success_

    call a%free_space()
    call a%psb_s_csr_sparse_mat%cp_from_coo(b, info)
    if (info /= 0) goto 9999
    call a%sync_space()
    call a%set_host()
    call a%sync()

    return

9999 continue
    info = psb_err_alloc_dealloc_
    return

  end subroutine psb_s_oacc_csr_cp_from_coo
end submodule psb_s_oacc_csr_cp_from_coo_impl
