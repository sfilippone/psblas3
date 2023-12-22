!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

subroutine psb_c_cp_hdia_to_coo(a,b,info) 
  
  use psb_base_mod
  use psb_c_hdia_mat_mod, psb_protect_name => psb_c_cp_hdia_to_coo
  use psi_ext_util_mod
  implicit none 

  class(psb_c_hdia_sparse_mat), intent(in)    :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  integer(psb_ipk_)   :: k,i,j,nc,nr,nza, nzd,nd,hacksize,nhacks,iszd,&
       &  ib, ir, kfirst, klast1, hackfirst, hacknext

  info = psb_success_
  if (a%is_dev()) call a%sync()

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat
  call b%set_nzeros(nza)
  call b%set_sort_status(psb_unsorted_)
  nhacks   = a%nhacks
  hacksize = a%hacksize
  j = 0 
  do k=1, nhacks
    i = (k-1)*hacksize + 1
    ib = min(hacksize,nr-i+1) 
    hackfirst = a%hackoffsets(k)
    hacknext  = a%hackoffsets(k+1)
    call psi_c_xtr_coo_from_dia(nr,nc,&
           & b%ia(j+1:), b%ja(j+1:), b%val(j+1:), nzd, &
           & hacksize,(hacknext-hackfirst),&
           & a%val((hacksize*hackfirst)+1:hacksize*hacknext),&
           & a%diaOffsets(hackfirst+1:hacknext),info,rdisp=(i-1))
!!$    write(*,*) 'diaoffsets',ib,' : ',ib - abs(a%diaOffsets(hackfirst+1:hacknext))
!!$    write(*,*) 'sum',ib,j,' : ',sum(ib - abs(a%diaOffsets(hackfirst+1:hacknext)))
    j = j + nzd
  end do
  if (nza /= j) then 
    write(*,*) 'Wrong counts in hdia_to_coo',j,nza
    info = -8 
    return 
  end if
  call b%set_host()
  call b%fix(info)

end subroutine psb_c_cp_hdia_to_coo
