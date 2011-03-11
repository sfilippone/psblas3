function psb_s_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1
  integer, intent(in) :: iv1(*)
  real(psb_spk_), intent(in) :: v1(*),v2(*)
  real(psb_spk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = szero 
  ip1 = 1

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_s_spge_dot

function psb_d_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1
  integer, intent(in) :: iv1(*)
  real(psb_dpk_), intent(in) :: v1(*),v2(*)
  real(psb_dpk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = dzero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_d_spge_dot

function psb_c_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1
  integer, intent(in) :: iv1(*)
  complex(psb_spk_), intent(in) :: v1(*),v2(*)
  complex(psb_spk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = czero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_c_spge_dot

function psb_z_spge_dot(nv1,iv1,v1,v2) result(dot) 
  use psb_const_mod
  integer, intent(in) :: nv1
  integer, intent(in) :: iv1(*)
  complex(psb_dpk_), intent(in) :: v1(*),v2(*)
  complex(psb_dpk_)      :: dot

  integer :: i,j,k, ip1, ip2

  dot = zzero 

  do i=1, nv1
    dot = dot + v1(i)*v2(iv1(i))
  end do

end function psb_z_spge_dot
