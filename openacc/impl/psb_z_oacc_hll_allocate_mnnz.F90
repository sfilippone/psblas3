submodule (psb_z_oacc_hll_mat_mod) psb_z_oacc_hll_allocate_mnnz_impl
  use psb_base_mod
contains
  module subroutine psb_z_oacc_hll_allocate_mnnz(m, n, a, nz)
    implicit none 
    integer(psb_ipk_), intent(in) :: m, n
    class(psb_z_oacc_hll_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in), optional :: nz
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: err_act, nz_
    character(len=20)  :: name='allocate_mnnz'
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: hksz, nhacks

    call psb_erractionsave(err_act)
    info = psb_success_

    if (present(nz)) then
      nz_ = nz
    else
      nz_ = 10  
    end if

    call a%psb_z_hll_sparse_mat%allocate(m, n, nz_)

    hksz = a%hksz
    nhacks = (m + hksz - 1) / hksz

    if (.not.allocated(a%val)) then
      allocate(a%val(nz_ * m))
      allocate(a%ja(nz_ * m))
      allocate(a%irn(m))
      allocate(a%idiag(m))
      allocate(a%hkoffs(nhacks))
    end if

    a%val = zzero
    a%ja = -1
    a%irn = 0
    a%idiag = 0
    a%hkoffs = 0

    call a%set_dev()
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_z_oacc_hll_allocate_mnnz
end submodule psb_z_oacc_hll_allocate_mnnz_impl
