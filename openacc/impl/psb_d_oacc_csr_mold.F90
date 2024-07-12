subroutine psb_d_oacc_csr_mold(a, b, info)
    use psb_base_mod
    use psb_d_oacc_csr_mat_mod, psb_protect_name => psb_d_oacc_csr_mold
    implicit none 
    class(psb_d_oacc_csr_sparse_mat), intent(in)                 :: a
    class(psb_d_base_sparse_mat), intent(inout), allocatable :: b
    integer(psb_ipk_), intent(out)                    :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='csr_mold'
    logical, parameter :: debug=.false.
  
    call psb_get_erraction(err_act)
    
    info = 0 
    if (allocated(b)) then 
      call b%free()
      deallocate(b, stat=info)
    end if
    if (info == 0) allocate(psb_d_oacc_csr_sparse_mat :: b, stat=info)
  
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_ 
      call psb_errpush(info, name)
      goto 9999
    end if
    return
  
  9999 call psb_error_handler(err_act)
  
    return
  
  end subroutine psb_d_oacc_csr_mold
  