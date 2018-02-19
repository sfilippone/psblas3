program tryidxijk
  use psb_base_mod
  use psb_util_mod

  integer(psb_lpk_) :: idx,idxm
  integer(psb_ipk_) :: nx,ny,nz
  integer(psb_ipk_) :: i,j,k, sidx

  idxm = 1000
  idxm = idxm*2000*1000
  nx = 2000
  ny = 2000
  nz = 2000
  do idx = idxm+300*1000*1000, idxm+300*1000*1000+50000
    call idx2ijk(i,j,k,idx,nx,ny,nz)
    sidx = idx
    write(*,*) 'idx2ijk: ',idx,i,j,k, sidx
  end do
end program tryidxijk
