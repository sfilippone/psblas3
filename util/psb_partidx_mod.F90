!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
! Purpose: 
!  Provide  functions to handle a distribution of a general
!  rectangular 2/3/n-dimensional domain onto
!  a rectangular 2/3/n-dimensional grid of processes
!
!  See test/pargen/psb_X_pdeNd for examples of usage
!
module psb_partidx_mod
  use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_

  interface idx2ijk
    module procedure idx2ijk3d, idx2ijkv, idx2ijk2d
  end interface idx2ijk

  interface ijk2idx
    module procedure ijk2idx3d, ijk2idxv, ijk2idx2d
  end interface ijk2idx
  interface idx2ijk
    module procedure lidx2ijk3d, lidx2ijkv, lidx2ijk2d,&
         & lidx2lijk3d, lidx2lijkv, lidx2lijk2d
  end interface idx2ijk

  interface ijk2idx
    module procedure ijk2lidx3d, ijk2lidxv, ijk2lidx2d,&
         & lijk2lidx3d, lijk2lidxv, lijk2lidx2d
  end interface ijk2idx

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
  subroutine  idx2ijk3d(i,j,k,idx,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: i,j,k
    integer(psb_mpk_), intent(in)  :: idx,nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: coords(3)

    call idx2ijk(coords,idx,[nx,ny,nz],base)
    
    k = coords(3)
    j = coords(2)
    i = coords(1)

  end subroutine idx2ijk3d
  
  subroutine  idx2ijk2d(i,j,idx,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: i,j
    integer(psb_mpk_), intent(in)  :: idx,nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: coords(2)

    call idx2ijk(coords,idx,[nx,ny],base)
    
    j = coords(2)
    i = coords(1)
    
  end subroutine idx2ijk2d
  
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
  subroutine  idx2ijkv(coords,idx,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: coords(:)
    integer(psb_mpk_), intent(in)  :: idx,dims(:)
    integer(psb_mpk_), intent(in), optional :: base

    integer(psb_mpk_) :: base_, idx_, i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if

    idx_ = idx - base_

    if (size(coords) < size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      coords = 0
      return
    end if

    !
    ! This code is equivalent to (3D case)
    ! k = mod(idx_,nz) + base_
    ! j = mod(idx_/nz,ny) + base_
    ! i = mod(idx_/(nx*ny),nx) + base_
    !
    do i=size(dims),1,-1
      coords(i) = mod(idx_,dims(i)) + base_
      idx_ = idx_ / dims(i)
    end do

  end subroutine idx2ijkv

  subroutine  lidx2ijk3d(i,j,k,idx,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: i,j,k
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_mpk_), intent(in)  :: nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: coords(3)

    call idx2ijk(coords,idx,[nx,ny,nz],base)
    
    k = coords(3)
    j = coords(2)
    i = coords(1)

  end subroutine lidx2ijk3d
  
  subroutine  lidx2ijk2d(i,j,idx,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: i,j
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_mpk_), intent(in)  :: nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: coords(2)

    call idx2ijk(coords,idx,[nx,ny],base)
    
    j = coords(2)
    i = coords(1)
    
  end subroutine lidx2ijk2d
  
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
  subroutine  lidx2ijkv(coords,idx,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: coords(:)
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_mpk_), intent(in)  :: dims(:)
    integer(psb_mpk_), intent(in), optional :: base

    integer(psb_epk_) :: base_, idx_
    integer(psb_mpk_) :: i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if

    idx_ = idx - base_

    if (size(coords) < size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      coords = 0
      return
    end if

    !
    ! This code is equivalent to (3D case)
    ! k = mod(idx_,nz) + base_
    ! j = mod(idx_/nz,ny) + base_
    ! i = mod(idx_/(nx*ny),nx) + base_
    !
    do i=size(dims),1,-1
      coords(i) = mod(idx_,dims(i)) + base_
      idx_ = idx_ / dims(i)
    end do

  end subroutine lidx2ijkv

  subroutine  lidx2lijk3d(i,j,k,idx,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: i,j,k
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_epk_), intent(in)  :: nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_epk_) :: coords(3)

    call idx2ijk(coords,idx,[nx,ny,nz],base)
    
    k = coords(3)
    j = coords(2)
    i = coords(1)

  end subroutine lidx2lijk3d
  
  subroutine  lidx2lijk2d(i,j,idx,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: i,j
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_epk_), intent(in)  :: nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_epk_) :: coords(2)

    call idx2ijk(coords,idx,[nx,ny],base)
    
    j = coords(2)
    i = coords(1)
    
  end subroutine lidx2lijk2d
  
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
  subroutine  lidx2lijkv(coords,idx,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: coords(:)
    integer(psb_epk_), intent(in)  :: idx
    integer(psb_epk_), intent(in)  :: dims(:)
    integer(psb_mpk_), intent(in), optional :: base

    integer(psb_epk_) :: base_, idx_
    integer(psb_epk_) :: i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if

    idx_ = idx - base_

    if (size(coords) < size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      coords = 0
      return
    end if

    !
    ! This code is equivalent to (3D case)
    ! k = mod(idx_,nz) + base_
    ! j = mod(idx_/nz,ny) + base_
    ! i = mod(idx_/(nx*ny),nx) + base_
    !
    do i=size(dims),1,-1
      coords(i) = mod(idx_,dims(i)) + base_
      idx_ = idx_ / dims(i)
    end do

  end subroutine lidx2lijkv
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
  subroutine  ijk2idxv(idx,coords,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(in)  :: coords(:),dims(:)
    integer(psb_mpk_), intent(out) :: idx
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: base_, i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if
    sz = size(coords)
    if (sz /= size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      idx = 0
      return
    end if

    idx = coords(1) - base_
    do i=2, sz
      idx = (idx * dims(i)) + coords(i) - base_
    end do
    idx = idx + base_
    
  end subroutine ijk2idxv
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
  subroutine  ijk2idx3d(idx,i,j,k,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: idx
    integer(psb_mpk_), intent(in)  :: i,j,k,nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*nz*ny + (j-base_)*nz + k - base_) + base_
    call ijk2idx(idx,[i,j,k],[nx,ny,nz],base)
  end subroutine ijk2idx3d

  subroutine  ijk2idx2d(idx,i,j,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(out) :: idx
    integer(psb_mpk_), intent(in)  :: i,j,nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*ny + (j-base_) + base_
    call ijk2idx(idx,[i,j],[nx,ny],base)
  end subroutine ijk2idx2d

  subroutine  ijk2lidxv(idx,coords,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_mpk_), intent(in)  :: coords(:),dims(:)
    integer(psb_epk_), intent(out) :: idx
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_mpk_) :: base_, i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if
    sz = size(coords)
    if (sz /= size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      idx = 0
      return
    end if

    idx = coords(1) - base_
    do i=2, sz
      idx = (idx * dims(i)) + coords(i) - base_
    end do
    idx = idx + base_
    
  end subroutine ijk2lidxv
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
  subroutine  ijk2lidx3d(idx,i,j,k,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: idx
    integer(psb_mpk_), intent(in)  :: i,j,k,nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*nz*ny + (j-base_)*nz + k - base_) + base_
    call ijk2idx(idx,[i,j,k],[nx,ny,nz],base)
  end subroutine ijk2lidx3d

  subroutine  ijk2lidx2d(idx,i,j,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: idx
    integer(psb_mpk_), intent(in)  :: i,j,nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*ny + (j-base_) + base_
    call ijk2idx(idx,[i,j],[nx,ny],base)
  end subroutine ijk2lidx2d


  subroutine  lijk2lidxv(idx,coords,dims,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(in)  :: coords(:),dims(:)
    integer(psb_epk_), intent(out) :: idx
    integer(psb_mpk_), intent(in), optional :: base
    
    integer(psb_epk_) :: base_, i, sz
    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if
    sz = size(coords)
    if (sz /= size(dims)) then
      write(0,*) 'Error: size mismatch ',size(coords),size(dims)
      idx = 0
      return
    end if

    idx = coords(1) - base_
    do i=2, sz
      idx = (idx * dims(i)) + coords(i) - base_
    end do
    idx = idx + base_
    
  end subroutine lijk2lidxv
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
  subroutine  lijk2lidx3d(idx,i,j,k,nx,ny,nz,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: idx
    integer(psb_epk_), intent(in)  :: i,j,k,nx,ny,nz
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*nz*ny + (j-base_)*nz + k - base_) + base_
    call ijk2idx(idx,[i,j,k],[nx,ny,nz],base)
  end subroutine lijk2lidx3d

  subroutine  lijk2lidx2d(idx,i,j,nx,ny,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_epk_), intent(out) :: idx
    integer(psb_epk_), intent(in)  :: i,j,nx,ny
    integer(psb_mpk_), intent(in), optional :: base
    
    !    idx = ((i-base_)*ny + (j-base_) + base_
    call ijk2idx(idx,[i,j],[nx,ny],base)
  end subroutine lijk2lidx2d
  
  !
  ! dist1Didx
  !   Given an index space [base : N-(1-base)] and
  !   a set of NP processes, split the index base as
  !   evenly as possible, i.e. difference in size 
  !   between any two processes is either 0 or 1,
  !   then return the boundaries in a vector
  !   such that
  !     V(P)   : first index owned by process P
  !     V(P+1) : first index owned by process P+1
  !
  subroutine dist1Didx(v,n,np,base)
    use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_
    implicit none 
    integer(psb_ipk_), intent(out) :: v(:)
    integer(psb_ipk_), intent(in)  :: n, np
    integer(psb_ipk_), intent(in), optional :: base
    !
    integer(psb_ipk_) :: base_, nb, i

    if (present(base)) then
      base_  = base
    else
      base_ = 1
    end if
    
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
  
end module psb_partidx_mod

