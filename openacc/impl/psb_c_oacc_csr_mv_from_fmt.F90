submodule (psb_c_oacc_csr_mat_mod) psb_c_oacc_csr_mv_from_fmt_impl
  use psb_base_mod
contains
  module subroutine psb_c_oacc_csr_mv_from_fmt(a, b, info)
    implicit none 

    class(psb_c_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_c_base_sparse_mat), intent(inout) :: b
    integer(psb_ipk_), intent(out)              :: info

    info = psb_success_

    select type(b)
    type is (psb_c_coo_sparse_mat)
      call a%mv_from_coo(b, info)
    class default
      call a%psb_c_csr_sparse_mat%mv_from_fmt(b, info)
      if (info /= 0) return

      !$acc update device(a%val, a%ja, a%irp)
    end select

  end subroutine psb_c_oacc_csr_mv_from_fmt
end submodule psb_c_oacc_csr_mv_from_fmt_impl
