subroutine psb_d_oacc_csr_allocate_mnnz(m, n, a, nz)
    use psb_base_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_allocate_mnnz
    implicit none 
    integer(psb_ipk_), intent(in) :: m, n
    class(psb_d_oacc_csr_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in), optional :: nz
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: err_act, nz_
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.
  
    call psb_erractionsave(err_act)
    info = psb_success_
  
    call a%psb_d_csr_sparse_mat%allocate(m, n, nz)
    
    if (.not.allocated(a%val)) then
      allocate(a%val(nz))
      allocate(a%ja(nz))
      allocate(a%irp(m+1))
    end if
  
    call a%set_dev()
    if (info /= 0) goto 9999
  
    call psb_erractionrestore(err_act)
    return
  
  9999 call psb_error_handler(err_act)
    return
  
end subroutine psb_d_oacc_csr_allocate_mnnz
  