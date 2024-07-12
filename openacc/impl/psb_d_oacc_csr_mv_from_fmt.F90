subroutine psb_d_oacc_csr_mv_from_fmt(a, b, info)
    use psb_base_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_mv_from_fmt
    implicit none 
  
    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_d_base_sparse_mat), intent(inout) :: b
    integer(psb_ipk_), intent(out)              :: info
  
    info = psb_success_
  
    select type(b)
    type is (psb_d_coo_sparse_mat)
      call a%mv_from_coo(b, info)
    class default
      call a%psb_d_csr_sparse_mat%mv_from_fmt(b, info)
      if (info /= 0) return
  
      !$acc update device(a%val, a%ja, a%irp)
    end select
  
end subroutine psb_d_oacc_csr_mv_from_fmt
  