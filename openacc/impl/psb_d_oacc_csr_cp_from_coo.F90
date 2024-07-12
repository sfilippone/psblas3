subroutine psb_d_oacc_csr_cp_from_coo(a, b, info)
    use psb_base_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_cp_from_coo
    implicit none
  
    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    class(psb_d_coo_sparse_mat), intent(in)         :: b
    integer(psb_ipk_), intent(out)                  :: info
  
    info = psb_success_
  
    call a%psb_d_csr_sparse_mat%cp_from_coo(b, info)
    if (info /= 0) goto 9999
  
    call a%set_dev()
    if (info /= 0) goto 9999
  
    return
  
  9999 continue
    info = psb_err_alloc_dealloc_
    return
  
  end subroutine psb_d_oacc_csr_cp_from_coo
  