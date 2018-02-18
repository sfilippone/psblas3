program tryyidxijk
  use psb_util_mod
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
  allocate(v(0:max(npx,npy,npz)))
  call dist1Didx(v,nx,npx)
  write(*,*) '  X:',v(0:npx)
  write(*,*) 'SZX:',v(1:npx)-v(0:npx-1)
  call dist1Didx(v,ny,npy)
  write(*,*) '  Y:',v(0:npy)
  write(*,*) 'SZY:',v(1:npy)-v(0:npy-1)
  call dist1Didx(v,nz,npz)
  write(*,*) '  Z:',v(0:npz)
  write(*,*) 'SZZ:',v(1:npz)-v(0:npz-1)
  
  
end program tryyidxijk
