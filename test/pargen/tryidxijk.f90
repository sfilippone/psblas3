program tryyidxijk
  use psb_base_mod, only : psb_ipk_
  implicit none 

  integer(psb_ipk_) :: nx,ny,nz, base, i, j, k, idx,npx,npy,npz
  integer(psb_ipk_) :: ic, jc, kc, idxc
  integer(psb_ipk_), allocatable :: v(:)
  write(*,*) 'nx,ny,nz,base? '
  read(*,*) nx,ny,nz, base

  idx = base
  do i=base,nx-(1-base)
    do j=base,ny-(1-base)
      do k=base,nz-(1-base)
        call idx2ijk(ic,jc,kc,idx,nx,ny,nz,base)
        call ijk2idx(idxc,i,j,k,nx,ny,nz,base)
        ! write(*,*) i,j,k,idx,':',ic,jc,kc,idxc
        if ((i/=ic).or.&
             & (j/=jc).or.&
             & (k/=kc).or.&
             & (idxc/=idx))then 
          write(*,*) 'Error:',i,j,k,idx,':',ic,jc,kc,idxc
          stop
        end if
        idx = idx + 1
      end do
    end do
  end do
  write(*,*) 'Ok '
  write(*,*) 'npx,npy,npz? '
  read(*,*) npx,npy,npz
  
  call dist1Didx(v,nx,npx)
  write(*,*) '  X:',v
  write(*,*) 'SZX:',v(2:npx+1)-v(1:npx)
  call dist1Didx(v,ny,npy)
  write(*,*) '  Y:',v
  write(*,*) 'SZY:',v(2:npy+1)-v(1:npy)
  call dist1Didx(v,nz,npz)
  write(*,*) '  Z:',v
  write(*,*) 'SZZ:',v(2:npz+1)-v(1:npz)
  
  
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
    
    k = mod(idx_,nz)         + base_
    j = mod(idx_/nz,ny)      + base_
    i = mod(idx_/(ny*nz),nx) + base_

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

    idx = ((i-base_)*nz*ny + (j-base_)*nz + k - base_) + base_
    
  end subroutine ijk2idx

  !
  ! dist1Didx
  !   Given an index space [base : N-(1-base)] and
  !   a set of NP processes, split the index base as
  !   evenly as possible, then return the boundaries
  !   in a vector such that
  !     V(P)   : first index owned by process P
  !     V(P+1) : first index owned by process P+1
  !
  subroutine dist1Didx(v,n,np,base)
    use psb_base_mod, only : psb_ipk_
    implicit none 
    integer(psb_ipk_), allocatable, intent(out) :: v(:)
    integer(psb_ipk_), intent(in)  :: n, np
    integer(psb_ipk_), intent(in), optional :: base
    !
    integer(psb_ipk_) :: base_, nb, i

    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if
    allocate(v(np+1))
    
    nb = n/np
    do i=1,mod(n,np)
      v(i) = nb + 1
    end do
    do i=mod(n,np)+1,np
      v(i) = nb 
    end do
    v(2:np+1)  = v(1:np)
    v(1) = base_
    do i=2,np+1
      v(i) = v(i) + v(i-1)
    end do
  end subroutine dist1Didx
  
end program tryyidxijk
