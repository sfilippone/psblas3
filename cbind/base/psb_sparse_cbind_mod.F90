! Shared validation for the distributed sparse C interfaces.
module psb_sparse_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  implicit none
  private
  public :: psb_c_sparse_get_desc, psb_c_sparse_check_pair
  public :: psb_c_sparse_same_owned, psb_c_sparse_agree

contains

  subroutine psb_c_sparse_get_desc(cdh, descp, info)
    type(psb_c_descriptor), intent(in) :: cdh
    type(psb_desc_type), pointer :: descp
    integer(psb_ipk_), intent(out) :: info

    nullify(descp)
    info = psb_err_invalid_cd_state_
    if (.not.c_associated(cdh%item)) return
    call c_f_pointer(cdh%item, descp)
    if (.not.descp%is_ok()) return
    info = psb_success_
  end subroutine psb_c_sparse_get_desc

  function psb_c_sparse_check_pair(desc_r, desc_c, assembled) result(info)
    type(psb_desc_type), intent(in) :: desc_r, desc_c
    logical, intent(in) :: assembled
    integer(psb_ipk_) :: info

    info = psb_err_context_error_
    if (.not.psb_cmp_ctxt(desc_r%get_context(), desc_c%get_context())) return
    info = psb_err_invalid_cd_state_
    if (.not.(desc_r%is_ok().and.desc_c%is_ok())) return
    if (assembled) then
      if (.not.(desc_r%is_asb().and.desc_c%is_asb())) return
    end if
    info = psb_success_
  end function psb_c_sparse_check_pair

  function psb_c_sparse_same_owned(desc_a, desc_b, info) result(same)
    type(psb_desc_type), intent(in) :: desc_a, desc_b
    integer(psb_ipk_), intent(out) :: info
    logical :: same
    integer(psb_lpk_), allocatable :: ia(:), ib(:)
    integer(psb_ipk_) :: i, nr, stat

    same = .false.
    info = psb_success_
    if (desc_a%get_global_rows() /= desc_b%get_global_rows()) return
    if (desc_a%get_local_rows() /= desc_b%get_local_rows()) return
    nr = desc_a%get_local_rows()
    allocate(ia(nr), ib(nr), stat=stat)
    if (stat /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    do i = 1, nr
      ia(i) = i
      ib(i) = i
    end do
    call desc_a%l2gip(ia, info, owned=.true.)
    if (info == psb_success_) call desc_b%l2gip(ib, info, owned=.true.)
    if (info /= psb_success_) return
    same = all(ia == ib)
  end function psb_c_sparse_same_owned

  subroutine psb_c_sparse_agree(desc, info)
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(inout) :: info

    ! Every participating rank must stop before entering a native collective
    ! if any rank has rejected the arguments or failed to allocate a workspace.
    if (info /= psb_success_) info = abs(info)
    call psb_max(desc%get_context(), info)
  end subroutine psb_c_sparse_agree

end module psb_sparse_cbind_mod
