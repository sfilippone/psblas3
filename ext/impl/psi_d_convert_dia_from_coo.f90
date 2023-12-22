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
subroutine psi_d_convert_dia_from_coo(a,tmp,info)
  use psb_base_mod
  use psb_d_dia_mat_mod, psb_protect_name => psi_d_convert_dia_from_coo
  use psi_ext_util_mod
  implicit none 
  class(psb_d_dia_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: tmp
  integer(psb_ipk_), intent(out)             :: info
  
  !locals
  integer(psb_ipk_)              :: ndiag,nd
  integer(psb_ipk_),allocatable  :: d(:)
  integer(psb_ipk_)              :: k,i,j,nc,nr,nza,ir,ic
  
  info = psb_success_
  nr  = tmp%get_nrows()
  nc  = tmp%get_ncols()
  nza = tmp%get_nzeros()
  ! If it is sorted then we can lessen memory impact 
  a%psb_d_base_sparse_mat = tmp%psb_d_base_sparse_mat
  
  ndiag = nr+nc-1
  allocate(d(ndiag),stat=info)
  if (info /= 0) return
  call psb_realloc(ndiag,a%offset,info)
  if (info /= 0) return
  
  call psi_dia_offset_from_coo(nr,nc,nza,tmp%ia,tmp%ja, &
       & nd,d,a%offset,info,initd=.true.,cleard=.false.) 
  
  call psb_realloc(nd,a%offset,info)
  if (info /= 0) return
  call psb_realloc(nr,nd,a%data,info) 
  if (info /= 0) return
  a%nzeros = nza
  
  call psi_xtr_dia_from_coo(nr,nc,nza,tmp%ia,tmp%ja,tmp%val,&
       & d,nr,nd,a%data,info,initdata=.true.)
  
  deallocate(d,stat=info)
  if (info /= 0) return
  
end subroutine psi_d_convert_dia_from_coo
