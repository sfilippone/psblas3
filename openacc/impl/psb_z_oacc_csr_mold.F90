submodule (psb_z_oacc_csr_mat_mod) psb_z_oacc_csr_mold_impl
  use psb_base_mod
contains  
  module subroutine psb_z_oacc_csr_mold(a, b, info)    
    implicit none 
    class(psb_z_oacc_csr_sparse_mat), intent(in)                 :: a
    class(psb_z_base_sparse_mat), intent(inout), allocatable :: b
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
    if (info == 0) allocate(psb_z_oacc_csr_sparse_mat :: b, stat=info)

    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_ 
      call psb_errpush(info, name)
      goto 9999
    end if
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_oacc_csr_mold
end submodule psb_z_oacc_csr_mold_impl

