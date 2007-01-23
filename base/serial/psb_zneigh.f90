!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
! File:  psb_zneigh.f90 
! Subroutine: 
! Parameters:

subroutine psb_zneigh(a,idx,neigh,n,info,lev)

  use psb_realloc_mod
  use psb_const_mod
  use psb_spmat_type
  implicit none


  type(psb_zspmat_type), intent(in) :: a   ! the sparse matrix
  integer, intent(in)               :: idx ! the index whose neighbours we want to find
  integer, intent(out)              :: n, info   ! the number of neighbours and the info
  integer, allocatable                  :: neigh(:) ! the neighbours
  integer, optional                 :: lev ! level of neighbours to find

  integer, allocatable   :: tmpn(:)
  integer :: level, dim, i, j, k, r, c, brow,&
       & elem_pt, ii, n1, col_idx, ne, err_act, nn, nidx
  character(len=20)                 :: name, ch_err

  name='psb_zneigh'
  info  = 0
  call psb_erractionsave(err_act)

  n    = 0
  info = 0
  if(present(lev)) then 
    if(lev.le.2) then
      level=lev
    else
      write(0,'("Too many levels!!!")')
      return
    endif
  else
    level=1
  end if

  call psb_zneigh1l(a,idx,neigh,n)
  
  if(level.eq.2) then
     n1=n
     allocate(tmpn(max(10,2*n)))
     if (size(neigh)<n*n1) call psb_realloc(n*n1,neigh,info)
     do i=1,n1
        nidx=neigh(i)
        if((nidx.ne.idx).and.(nidx.gt.0).and.(nidx.le.a%m)) then
           call psb_zneigh1l(a,nidx,tmpn,nn)
           if((n+nn).gt.size(neigh)) call psb_realloc(n+nn,neigh,info)
           neigh(n+1:n+nn)=tmpn(1:nn)
           n=n+nn
        end if
     end do
     deallocate(tmpn)
  end if
  
  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

     
contains

  subroutine csr_zneigh1l(a,idx,neigh,n)
    
    use psb_realloc_mod
    use psb_const_mod
    use psb_spmat_type
    implicit none
    
    
    type(psb_zspmat_type), intent(in) :: a   ! the sparse matrix
    integer, intent(in)               :: idx ! the index whose neighbours we want to find
    integer, intent(out)              :: n   ! the number of neighbours and the info
    integer, allocatable                  :: neigh(:) ! the neighbours
    
    integer :: dim, i, iidx

    if(a%pl(1).ne.0) then
       iidx=a%pl(idx)
    else
       iidx=idx
    end if

    dim=a%ia2(iidx+1)-a%ia2(iidx)
    if(dim.gt.size(neigh)) call psb_realloc(dim,neigh,info)

    n=0
    do i=a%ia2(iidx), a%ia2(iidx+1)-1
       n=n+1
       neigh(n)=a%ia1(i)
    end do
    
  end subroutine csr_zneigh1l


  subroutine coo_zneigh1l(a,idx,neigh,n)

    use psb_realloc_mod
    use psb_const_mod
    use psb_spmat_type
    implicit none


    type(psb_zspmat_type), intent(in) :: a   ! the sparse matrix
    integer, intent(in)               :: idx ! the index whose neighbours we want to find
    integer, intent(out)              :: n   ! the number of neighbours and the info
    integer, allocatable                  :: neigh(:) ! the neighbours

    integer :: dim, i, iidx, ip, nza

    if(a%pl(1).ne.0) then
      iidx=a%pl(idx)
    else
      iidx=idx
    end if

    nza=a%infoa(psb_nnz_)

    if (a%infoa(psb_srtd_) == psb_isrtdcoo_) then 
      call ibsrch(ip,iidx,nza,a%ia1)
      if (ip /= -1) then 
        ! bring ip backward to the beginning of the row
        do 
          if (ip < 2) exit
          if (a%ia1(ip-1) == iidx) then  
            ip = ip -1 
          else 
            exit
          end if
        end do
      end if

      dim=0
      do
        if(a%ia1(ip).eq.iidx) then
          dim=dim+1
          if(dim.gt.size(neigh)) call psb_realloc(dim*3/2,neigh,info)
          neigh(dim)=a%ia2(ip)
          ip=ip+1
        else
          exit
        end if
      end do

    else

      dim=0
      do i=1,nza
        if(a%ia1(i).eq.iidx) then
          dim=dim+1
          if(dim.gt.size(neigh)) call psb_realloc(dim*3/2,neigh,info)
          neigh(dim)=a%ia2(ip)
        end if
      end do
    end if

    n=dim

  end subroutine coo_zneigh1l
    


  subroutine jad_zneigh1l(a,idx,neigh,n)
    
    use psb_realloc_mod
    use psb_const_mod
    use psb_spmat_type
    implicit none
    
    
    type(psb_zspmat_type), intent(in),target :: a   ! the sparse matrix
    integer, intent(in)               :: idx ! the index whose neighbours we want to find
    integer, intent(out)              :: n   ! the number of neighbours and the info
    integer, allocatable              :: neigh(:) ! the neighbours
    
    integer :: dim, i, iidx, ip, nza
    integer, pointer                      :: ia1(:), ia2(:), ia3(:),&
         & ja(:), ka(:)
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ng

    if(a%pl(1).ne.0) then
       iidx=a%pl(idx)
    else
       iidx=idx
    end if

    nza=a%infoa(psb_nnz_)

    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)
    
    ng  =  a%ia2(png)              ! the number of blocks
    ja  => a%ia2(pja:)             ! the array containing the pointers to ka and aspk
    ka  => a%ia1(:)                ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. in ja to the first csr column


    i=0
    dim=0
    blkfnd: do
       i=i+1
       if(ia1(i).eq.iidx) then
          blk=i
          dim=dim+ia3(i)-ia2(i)
          ipx = ia1(i)         ! the first row index of the block
          rb  = iidx-ipx   ! the row offset within the block
          row = ia3(i)+rb
          dim  = dim+ja(row+1)-ja(row)
          exit blkfnd
       else if(ia1(i).gt.iidx) then
          blk=i-1
          dim=dim+ia3(i-1)-ia2(i-1)
          ipx = ia1(i-1)         ! the first row index of the block
          rb  = iidx-ipx   ! the row offset within the block
          row = ia3(i-1)+rb
          dim  = dim+ja(row+1)-ja(row)
          exit blkfnd
       end if
    end do blkfnd
    
    if(dim.gt.size(neigh)) call psb_realloc(dim,neigh,info)
    
    ipx = ia1(blk)             ! the first row index of the block
    k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
    rb  = iidx-ipx       ! the row offset within the block
    npg = ja(k_pt+1)-ja(k_pt)  ! the number of rows in the block
    
    k=0
    do  col = ia2(blk), ia3(blk)-1 
       k=k+1
       neigh(k)  = ka(ja(col)+rb)
    end do
    
    ! extract second part of the row from the csr tail
    row=ia3(blk)+rb
    do j=ja(row), ja(row+1)-1
       k=k+1
       neigh(k)  = ka(j)
    end do
    
    n=k

  end subroutine jad_zneigh1l

  subroutine psb_zneigh1l(a,idx,neigh,n)
    
    use psb_realloc_mod
    use psb_const_mod
    use psb_spmat_type
    use psb_string_mod
    implicit none
    
    
    type(psb_zspmat_type), intent(in) :: a   ! the sparse matrix
    integer, intent(in)               :: idx ! the index whose neighbours we want to find
    integer, intent(out)              :: n   ! the number of neighbours and the info
    integer, allocatable                  :: neigh(:) ! the neighbours


    select case(toupper(a%fida(1:3)))
    case('CSR')
       call csr_zneigh1l(a,idx,neigh,n)
    case('COO')
       call coo_zneigh1l(a,idx,neigh,n)
    case('JAD')
       call jad_zneigh1l(a,idx,neigh,n)
    end select

  end subroutine psb_zneigh1l

end subroutine psb_zneigh
  
