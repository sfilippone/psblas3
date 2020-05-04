module psb_base_util_cbind_mod
  use iso_c_binding
  use psb_util_mod
  use psb_cutil_cbind_mod
  use psb_dutil_cbind_mod
  use psb_sutil_cbind_mod
  use psb_zutil_cbind_mod

contains

  ! Routines for managing indexes are type independent
  ! so we have them defined only in the common module
  ! for all the index lengths:

  function psb_c_i_ijk2idx(ijk,sizes,modes,base) bind(c) result(idx)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    use psb_util_mod
    implicit none

    integer(psb_c_ipk_)        :: idx
    integer(psb_c_ipk_), value :: modes, base
    integer(psb_c_ipk_)        :: ijk(modes)
    integer(psb_c_ipk_)        :: sizes(modes)

    integer(psb_ipk_)          :: fijk(modes), fsizes(modes)

    fijk(1:modes) = ijk(1:modes)
    fsizes(1:modes) = sizes(1:modes)

    call ijk2idx(idx,fijk,fsizes,base)

  end function psb_c_i_ijk2idx

  function psb_c_l_ijk2idx(ijk,sizes,modes,base) bind(c) result(idx)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    use psb_util_mod
    implicit none

    integer(psb_c_lpk_)        :: idx
    integer(psb_c_ipk_), value :: modes, base
    integer(psb_c_ipk_)        :: ijk(modes)
    integer(psb_c_ipk_)        :: sizes(modes)

    integer(psb_ipk_)          :: fijk(modes), fsizes(modes)

    fijk(1:modes) = ijk(1:modes)
    fsizes(1:modes) = sizes(1:modes)

    call ijk2idx(idx,fijk,fsizes,base)

  end function psb_c_l_ijk2idx

  function psb_c_i_idx2ijk(ijk,idx,sizes,modes,base) bind(c) result(res)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none

    integer(psb_c_ipk_)        :: res
    integer(psb_c_ipk_), value :: idx
    integer(psb_c_ipk_), value :: modes, base
    integer(psb_c_ipk_)        :: ijk(modes)
    integer(psb_c_ipk_)        :: sizes(modes)

    integer(psb_ipk_)          :: fijk(modes), fsizes(modes)

    res = -1

    fsizes(1:modes) = sizes(1:modes)
    call idx2ijk(fijk,idx,fsizes,base=base)

    ijk(1:modes) = fijk(1:modes)

    res = 0

  end function

  function psb_c_l_idx2ijk(ijk,idx,sizes,modes,base) bind(c) result(res)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none

    integer(psb_c_ipk_)        :: res
    integer(psb_c_lpk_), value :: idx
    integer(psb_c_ipk_), value :: modes, base
    integer(psb_c_ipk_)        :: ijk(modes)
    integer(psb_c_ipk_)        :: sizes(modes)

    integer(psb_ipk_)  :: fijk(modes), fsizes(modes)

    res = -1

    fsizes(1:modes) = sizes(1:modes)
    call idx2ijk(fijk,idx,fsizes,base=base)

    ijk(1:modes) = fijk(1:modes)

    res = 0

  end function

end module psb_base_util_cbind_mod
