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

subroutine psb_z_hdia_csmv(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_z_hdia_mat_mod, psb_protect_name => psb_z_hdia_csmv
  implicit none 
  class(psb_z_hdia_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)        :: alpha, beta, x(:)
  complex(psb_dpk_), intent(inout)     :: y(:)
  integer(psb_ipk_), intent(out)     :: info
  character, optional, intent(in)    :: trans

  character :: trans_
  integer(psb_ipk_)  :: i,j,k,m,n, nnz, ir, jc,nr,nc
  integer(psb_ipk_)  :: irs,ics, nmx, ni
  integer(psb_ipk_)  :: nhacks, hacksize,maxnzhack, ncd,ib, nzhack, &
       & hackfirst, hacknext
  logical            :: tra, ctra
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_hdia_csmv'
  logical, parameter :: debug=.false.
  real :: start, finish
  call psb_erractionsave(err_act)
  info = psb_success_

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (a%is_dev()) call a%sync()
  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')
  if (tra.or.ctra) then 
    m = a%get_ncols()
    n = a%get_nrows()
    info = psb_err_transpose_not_n_unsupported_
    call psb_errpush(info,name)
    goto 9999
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if
  
  if (size(x,1)<n) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3*ione,n/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5*ione,m/))
    goto 9999
  end if

  nhacks   = a%nhacks
  hacksize = a%hacksize
  
  do k=1, nhacks
    i = (k-1)*hacksize + 1
    ib = min(hacksize,m-i+1) 
    hackfirst = a%hackoffsets(k)
    hacknext  = a%hackoffsets(k+1)
    ncd = hacknext-hackfirst
    
    call psi_z_inner_dia_csmv(m,n,&
         & alpha,hacksize,ncd,&
         & a%val((hacksize*hackfirst)+1:hacksize*hacknext),&
         & a%diaOffsets(hackfirst+1:hacknext),x,beta,y,info,rdisp=(i-1))

  end do
    
  
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psi_z_inner_dia_csmv(nr,nc,alpha,nrd,ncd,data,offsets,&
       & x,beta,y,info,rdisp) 
    implicit none 
    integer(psb_ipk_), intent(in)  :: nr,nc,nrd,ncd,offsets(*)
    integer(psb_ipk_)              :: rdisp, info
    complex(psb_dpk_), intent(in)     :: alpha, beta, x(*),data(nrd,ncd)
    complex(psb_dpk_), intent(inout)  :: y(*)
        

    integer(psb_ipk_) :: i,j,k, ir, jc, m4, ir1, ir2, nrcmdisp, rdisp1
    
    info = 0
    nrcmdisp = min(nr-rdisp,nc-rdisp) 
    rdisp1   = 1-rdisp
    if (beta == dzero) then
      do i = 1, min(nrd,nr-rdisp)
        y(rdisp+i) = dzero
      enddo
    else
      do  i = 1, min(nrd,nr-rdisp)
        y(rdisp+i) = beta*y(i)
      end do
    endif
    do j=1, ncd
      if (offsets(j)>=0) then 
        ir1 = 1
        !  min(nrd,nr - offsets(j) - rdisp_,nc-offsets(j)-rdisp_)
        ir2 = min(nrd, nrcmdisp - offsets(j))
      else
        !  max(1,1-offsets(j)-rdisp_) 
        ir1 = max(1, rdisp1 - offsets(j))
        ir2 = min(nrd, nrcmdisp) 
      end if
      jc = ir1 + rdisp + offsets(j)
      do i=ir1,ir2 
        y(rdisp+i) = y(rdisp+i) + alpha*data(i,j)*x(jc)
        jc = jc + 1 
      enddo
    end do
  end subroutine psi_z_inner_dia_csmv

end subroutine psb_z_hdia_csmv
