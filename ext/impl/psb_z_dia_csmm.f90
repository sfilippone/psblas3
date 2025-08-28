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
subroutine psb_z_dia_csmm(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_z_dia_mat_mod, psb_protect_name => psb_z_dia_csmm
  implicit none 
  class(psb_z_dia_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  complex(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)       :: info
  character, optional, intent(in)      :: trans

  character :: trans_
  integer(psb_ipk_)   :: i,j,k,m,n, nnz, ir, jc, nxy
  logical   :: tra, ctra
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='z_dia_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (a%is_dev()) call a%sync()

  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')
  if (tra.or.ctra) then 
    m = a%get_ncols()
    n = a%get_nrows()
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

  nxy = min(size(x,2) , size(y,2) )

  call  psb_z_dia_csmm_inner(m,n,nxy,alpha,&
       & a%data,size(a%data,1,kind=psb_ipk_), size(a%data,2,kind=psb_ipk_), a%offset,&
       & x,size(x,1,kind=psb_ipk_), beta, y,size(y,1,kind=psb_ipk_))

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_z_dia_csmm_inner(m,n,nxy,alpha,data,nr,nc,off,&
       &x,ldx,beta,y,ldy) 
    integer(psb_ipk_), intent(in)   :: m,n,nr,nc,off(*), ldx,ldy,nxy
    complex(psb_dpk_), intent(in)     :: alpha, beta, x(ldx,*),data(nr,*)
    complex(psb_dpk_), intent(inout)  :: y(ldy,*)

    integer(psb_ipk_) :: i,j,k, ir, jc, m4, ir1, ir2

    if (beta == zzero) then
      do i = 1, m
        y(i,1:nxy) = zzero
      enddo
    else
      do  i = 1, m
        y(i,1:nxy) = beta*y(i,1:nxy)
      end do
    endif

    do j=1,nc
      if (off(j) > 0) then 
        ir1 = 1
        ir2 = nr - off(j)
      else
        ir1 = 1 - off(j)
        ir2 = nr
      end if
      do i=ir1, ir2
        y(i,1:nxy) = y(i,1:nxy) + alpha*data(i,j)*x(i+off(j),1:nxy)
      enddo
    enddo

  end subroutine psb_z_dia_csmm_inner

end subroutine psb_z_dia_csmm
