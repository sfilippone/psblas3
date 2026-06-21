module psb_vect_mod
  use psb_i2_vect_mod
  use psb_i_vect_mod
  use psb_l_vect_mod  
  use psb_s_vect_mod
  use psb_d_vect_mod
  use psb_c_vect_mod
  use psb_z_vect_mod
  use psb_i_multivect_mod
  use psb_l_multivect_mod
  use psb_s_multivect_mod
  use psb_d_multivect_mod
  use psb_c_multivect_mod
  use psb_z_multivect_mod

contains

  subroutine psb_init_vect_defaults()
    implicit none
    !
    ! Defaults for vectors 
    !
    type(psb_i_base_vect_type)  :: ivetdef
    type(psb_l_base_vect_type)  :: lvetdef
    type(psb_s_base_vect_type)  :: svetdef
    type(psb_d_base_vect_type)  :: dvetdef
    type(psb_c_base_vect_type)  :: cvetdef
    type(psb_z_base_vect_type)  :: zvetdef
    type(psb_i_base_multivect_type)  :: imvetdef
    type(psb_l_base_multivect_type)  :: lmvetdef
    type(psb_s_base_multivect_type)  :: smvetdef
    type(psb_d_base_multivect_type)  :: dmvetdef
    type(psb_c_base_multivect_type)  :: cmvetdef
    type(psb_z_base_multivect_type)  :: zmvetdef

    call psb_i_set_vect_default(ivetdef)
    call psb_l_set_vect_default(lvetdef)
    call psb_s_set_vect_default(svetdef)
    call psb_d_set_vect_default(dvetdef)
    call psb_c_set_vect_default(cvetdef)
    call psb_z_set_vect_default(zvetdef)

    call psb_i_set_multivect_default(imvetdef)
    call psb_l_set_multivect_default(lmvetdef)
    call psb_s_set_multivect_default(smvetdef)
    call psb_d_set_multivect_default(dmvetdef)
    call psb_c_set_multivect_default(cmvetdef)
    call psb_z_set_multivect_default(zmvetdef)

  end subroutine psb_init_vect_defaults

  subroutine psb_clear_vect_defaults()
    implicit none

    call psb_i_clear_vect_default()
    call psb_l_clear_vect_default()
    call psb_s_clear_vect_default()
    call psb_d_clear_vect_default()
    call psb_c_clear_vect_default()
    call psb_z_clear_vect_default()

    call psb_i_clear_multivect_default()
    call psb_l_clear_multivect_default()
    call psb_s_clear_multivect_default()
    call psb_d_clear_multivect_default()
    call psb_c_clear_multivect_default()
    call psb_z_clear_multivect_default()

  end subroutine psb_clear_vect_defaults

end module psb_vect_mod
