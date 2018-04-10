module psb_mat_mod
  use psb_s_mat_mod
  use psb_d_mat_mod
  use psb_c_mat_mod
  use psb_z_mat_mod

contains

  subroutine psb_init_mat_defaults()
    implicit none
    !
    ! Defaults for matrices
    !
    type(psb_s_csr_sparse_mat) :: smatdef
    type(psb_d_csr_sparse_mat) :: dmatdef
    type(psb_c_csr_sparse_mat) :: cmatdef
    type(psb_z_csr_sparse_mat) :: zmatdef

    call psb_set_mat_default(smatdef)
    call psb_set_mat_default(dmatdef)
    call psb_set_mat_default(cmatdef)
    call psb_set_mat_default(zmatdef)
    
  end subroutine psb_init_mat_defaults

end module psb_mat_mod
