module psb_base_util_cbind_mod
  use iso_c_binding
  use psb_cutil_cbind_mod
  use psb_dutil_cbind_mod
  use psb_sutil_cbind_mod
  use psb_zutil_cbind_mod

contains

  ! Routines for managing indexes are type independent
  ! so we have them defined only in the common module
  ! for all the types

  function psb_c_idx2ijk(i,j,idx,nx,ny,base) bind(c) result(res)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none

    integer(psb_c_ipk_)         :: res
    integer(psb_c_ipk_), value  :: idx,nx,ny,base
    integer(psb_c_ipk_)         :: i,j

    res = -1

    call idx2ijk(i,j,idx,nx,ny,base=base)

    res = 0

  end function psb_c_idx2ijk

  function psb_c_lidx2ijk(i,j,idx,nx,ny,base) bind(c) result(res)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none

    integer(psb_c_ipk_)         :: res
    integer(psb_c_lpk_), value  :: idx 
    integer(psb_c_ipk_), value  :: nx,ny,base
    integer(psb_c_ipk_)         :: i,j

    res = -1

    call idx2ijk(i,j,idx,nx,ny,base=base)

    res = 0

  end function psb_c_lidx2ijk

end module psb_base_util_cbind_mod
