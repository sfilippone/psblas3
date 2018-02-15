program tryyidxijk
  use psb_base_mod, only : psb_ipk_
  implicit none 

  integer(psb_ipk_) :: nx,ny,nz, base, i, j, k, idx
  integer(psb_ipk_) :: ic, jc, kc, idxc
  write(*,*) 'nx,ny,nz,base? '
  read(*,*) nx,ny,nz, base

  idx = base
  do i=base,nx-(1-base)
    do j=base,ny-(1-base)
      do k=base,nz-(1-base)
        call idx2ijk(ic,jc,kc,idx,nx,ny,nz,base)
        call ijk2idx(idxc,i,j,k,nx,ny,nz,base)
        write(*,*) i,j,k,idx,':',ic,jc,kc,idxc
        idx = idx + 1
      end do
    end do
  end do
  
contains
  !
  ! Given  a global index IDX and the domain size (NX,NY,NZ)
  ! compute the point coordinates (I,J,K) 
  ! Optional argument: base 0 or 1, default 1
  !
  ! This mapping is equivalent to a loop nesting:
  !  idx = base
  !  do i=1,nx
  !    do j=1,ny
  !      do k=1,nz
  !         ijk2idx(i,j,k) = idx
  !         idx = idx + 1
  subroutine  idx2ijk(i,j,k,idx,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_
    implicit none 
    integer(psb_ipk_), intent(out) :: i,j,k
    integer(psb_ipk_), intent(in)  :: idx,nx,ny,nz
    integer(psb_ipk_), intent(in), optional :: base
    
    integer(psb_ipk_) :: base_, idx_
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if

    idx_ = idx - base_
    
    i = idx_/(nx*ny)
    j = (idx_ - i*nx*ny)/ny
    k = (idx_ - i*nx*ny -j*ny) 

    i = i + base_
    j = j + base_
    k = k + base_
  end subroutine idx2ijk
  
  !
  ! Given  a triple (I,J,K) and  the domain size (NX,NY,NZ)
  ! compute the global index IDX 
  ! Optional argument: base 0 or 1, default 1
  !
  ! This mapping is equivalent to a loop nesting:
  !  idx = base
  !  do i=1,nx
  !    do j=1,ny
  !      do k=1,nz
  !         ijk2idx(i,j,k) = idx
  !         idx = idx + 1
  subroutine  ijk2idx(idx,i,j,k,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_
    implicit none 
    integer(psb_ipk_), intent(out) :: idx
    integer(psb_ipk_), intent(in)  :: i,j,k,nx,ny,nz
    integer(psb_ipk_), intent(in), optional :: base
    
    integer(psb_ipk_) :: base_
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if


    idx = ((i-base_)*nx*ny + (j-base_)*nx + k -base_) + base_
    
  end subroutine ijk2idx
end program tryyidxijk
