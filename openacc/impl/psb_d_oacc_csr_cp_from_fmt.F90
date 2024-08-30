submodule (psb_d_oacc_csr_mat_mod) psb_d_oacc_csr_cp_from_fmt_impl
  use psb_base_mod
contains
  module subroutine psb_d_oacc_csr_cp_from_fmt(a, b, info)
    implicit none 

    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_d_base_sparse_mat), intent(in) :: b
    integer(psb_ipk_), intent(out)            :: info

    info = psb_success_

    select type(b)
    type is (psb_d_coo_sparse_mat)
      call a%cp_from_coo(b, info)
    class default
      call a%free_dev_space()
      call a%psb_d_csr_sparse_mat%cp_from_fmt(b, info)
      if (info /= 0) return
      call a%sync_dev_space()
      call a%set_host()
      call a%sync()
    end select

  end subroutine psb_d_oacc_csr_cp_from_fmt
end submodule psb_d_oacc_csr_cp_from_fmt_impl
