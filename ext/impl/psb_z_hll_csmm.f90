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
  

subroutine psb_z_hll_csmm(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_z_hll_mat_mod, psb_protect_name => psb_z_hll_csmm
  implicit none 
  class(psb_z_hll_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)              :: alpha, beta, x(:,:)
  complex(psb_dpk_), intent(inout)           :: y(:,:)
  integer(psb_ipk_), intent(out)          :: info
  character, optional, intent(in)         :: trans

  character :: trans_
  integer(psb_ipk_)   :: i,j,k,m,n, nnz, ir, jc, nxy,ldx,ldy,hksz,mxrwl
  complex(psb_dpk_), allocatable  :: acc(:)
  logical            :: tra, ctra
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='z_hll_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  nxy = min(size(x,2) , size(y,2) )


  ldx = size(x,1)
  ldy = size(y,1)
  if (a%is_dev()) call a%sync()
  
  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')


  if (tra.or.ctra) then 

    m = a%get_ncols()
    n = a%get_nrows()
    if (ldx<n) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3*ione,n/))
      goto 9999
    end if

    if (ldy<m) then 
      info = 36
      call psb_errpush(info,name,i_err=(/5*ione,m/))
      goto 9999
    end if

    if (beta == zzero) then
      do i = 1, m
        y(i,1:nxy) = zzero
      enddo
    else
      do  i = 1, m
        y(i,1:nxy) = beta*y(i,1:nxy)
      end do
    endif

    hksz = a%get_hksz()
    j=1
    do i=1,n,hksz
      ir    = min(hksz,n-i+1) 
      mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
      k     = a%hkoffs(j) + 1
      call psb_z_hll_csmm_inner(i,ir,nxy,mxrwl,a%irn(i),&
           & alpha,a%ja(k),hksz,a%val(k),hksz,&
           & a%is_triangle(),a%is_unit(),&
           & x,ldx,zone,y,ldy,tra,ctra,info) 
      if (info /= psb_success_) goto 9999
      j = j + 1 
    end do


  else if (.not.tra) then 

    n = a%get_ncols()
    m = a%get_nrows()

    if (ldx<n) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3*ione,n/))
      goto 9999
    end if

    if (ldy<m) then 
      info = 36
      call psb_errpush(info,name,i_err=(/5*ione,m/))
      goto 9999
    end if


    hksz = a%get_hksz()
    j=1
    do i=1,m,hksz
      ir    = min(hksz,m-i+1) 
      mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
      k     = a%hkoffs(j) + 1
      call psb_z_hll_csmm_inner(i,ir,nxy,mxrwl,a%irn(i),&
           & alpha,a%ja(k),hksz,a%val(k),hksz,&
           & a%is_triangle(),a%is_unit(),&
           & x,ldx,beta,y,ldy,tra,ctra,info) 
      if (info /= psb_success_) goto 9999
      j = j + 1 
    end do

  end if


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_z_hll_csmm_inner(ir,m,nc,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit,x,ldx,beta,y,ldy,tra,ctra,info) 
    integer(psb_ipk_), intent(in)   :: ir,m,n,nc,ldj,ldv,ja(ldj,*),irn(*),ldx,ldy
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(ldx,*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(ldy,*)
    logical, intent(in)             :: is_triangle, is_unit, tra, ctra
    integer(psb_ipk_), intent(out)  :: info

    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)    :: acc(4), tmp(nc)

    info = psb_success_
    if (tra) then 

      if (beta == zone) then 
        do i=1,m
          do j=1, irn(i)
            jc = ja(i,j)
            y(jc,1:nc) = y(jc,1:nc) + alpha*val(i,j)*x(ir+i-1,1:nc)
          end do
        end do
      else
        info = -10

      end if

    else if (ctra) then 

      if (beta == zone) then 
        do i=1,m
          do j=1, irn(i)
            jc = ja(i,j)
            y(jc,1:nc) = y(jc,1:nc) + alpha*conjg(val(i,j))*x(ir+i-1,1:nc)
          end do
        end do
      else
        info = -10
      end if

    else if (.not.(tra.or.ctra)) then 

      if (alpha == zzero) then 
        if (beta == zzero) then 
          do i=1,m
            y(ir+i-1,1:nc) = zzero
          end do
        else
          do i=1,m
            y(ir+i-1,1:nc) =  beta*y(ir+i-1,1:nc) 
          end do
        end if

      else
        if (beta == zzero) then 
          do i=1,m
            tmp(1:nc) = zzero
            do j=1, irn(i)
              tmp(1:nc) = tmp(1:nc) + val(i,j)*x(ja(i,j),1:nc)
            end do
            y(ir+i-1,1:nc) = alpha*tmp(1:nc) 
          end do
        else
          do i=1,m
            tmp(1:nc) = zzero
            do j=1, irn(i)
              tmp(1:nc) = tmp(1:nc) + val(i,j)*x(ja(i,j),1:nc)
            end do
            y(ir+i-1,1:nc) = alpha*tmp(1:nc) + beta*y(ir+i-1,1:nc)
          end do
        endif
      end if
    end if

    if (is_unit) then 
      do i=1, min(m,n)
        y(ir+i-1,1:nc) = y(ir+i-1,1:nc) + alpha*x(ir+i-1,1:nc)
      end do
    end if

  end subroutine psb_z_hll_csmm_inner
end subroutine psb_z_hll_csmm
