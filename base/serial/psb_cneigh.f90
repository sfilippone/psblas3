!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
! File:  psb_cneigh.f90 
! Subroutine: 
! Arguments:

subroutine psb_cneigh(a,idx,neigh,n,info,lev)

  use psb_realloc_mod
  use psb_const_mod
  use psb_spmat_type
  use psb_serial_mod, psb_protect_name => psb_cneigh
  implicit none


  type(psb_cspmat_type), intent(in) :: a   ! the sparse matrix
  integer, intent(in)               :: idx ! the index whose neighbours we want to find
  integer, intent(out)              :: n, info   ! the number of neighbours and the info
  integer, allocatable              :: neigh(:) ! the neighbours
  integer, optional, intent(in)     :: lev ! level of neighbours to find

  integer :: lev_, i, nl, ifl,ill,&
       &  n1, err_act, nn, nidx,ntl
  character(len=20)                 :: name
  integer, allocatable :: ia(:), ja(:)
  complex(psb_spk_), allocatable :: val(:)

  name='psb_cneigh'
  info  = 0
  call psb_erractionsave(err_act)

  n    = 0
  info = 0
  if(present(lev)) then 
    lev_ = lev
  else
    lev_=1
  end if

  call psb_sp_getrow(idx,a,n,ia,ja,val,info)
  if (info == 0) call psb_realloc(n,neigh,info)
  if (info /= 0) then 
    call psb_errpush(4000,name)
    goto 9999
  end if
  neigh(1:n) = ja(1:n)
  ifl = 1
  ill = n

  
  do nl = 2, lev_ 
    n1 = ill - ifl + 1
    call psb_ensure_size(ill+n1*n1,neigh,info)
    if (info /= 0) then 
      call psb_errpush(4000,name)
      goto 9999
    end if
    ntl = 0
    do i=ifl,ill
      nidx=neigh(i)
      if ((nidx /= idx).and.(nidx > 0).and.(nidx <= a%m)) then
        call psb_sp_getrow(nidx,a,nn,ia,ja,val,info)
        if (info==0) call psb_ensure_size(ill+ntl+nn,neigh,info)
        if (info /= 0) then 
          call psb_errpush(4000,name)
          goto 9999
        end if        
        neigh(ill+ntl+1:ill+ntl+nn)=ja(1:nn)
        ntl = ntl+nn
      end if
    end do
    call psb_msort_unique(neigh(ill+1:ill+ntl),nn)
    ifl = ill + 1
    ill = ill + nn
  end do
  call imsru(ill,neigh,psb_sort_up_,nn)
  n = nn

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error()
     return
  end if
  return


end subroutine psb_cneigh
  
