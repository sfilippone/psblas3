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
  

subroutine psb_d_ell_csmm(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csmm
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)       :: info
  character, optional, intent(in)      :: trans

  character :: trans_
  integer(psb_ipk_)   :: i,j,k,m,n, nnz, ir, jc, nxy
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra, ctra
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='d_ell_csmm'
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

  allocate(acc(nxy), stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='allocate')
    goto 9999
  end if

  call  psb_d_ell_csmm_inner(m,n,nxy,alpha,size(a%ja,2,kind=psb_ipk_),&
       & a%ja,size(a%ja,1,kind=psb_ipk_),a%val,size(a%val,1,kind=psb_ipk_), &
       & a%is_triangle(),a%is_unit(),x,size(x,1,kind=psb_ipk_), &
       & beta,y,size(y,1,kind=psb_ipk_),tra,ctra,acc) 


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains
  subroutine psb_d_ell_csmm_inner(m,n,nxy,alpha,nc,ja,ldj,val,ldv,&
       & is_triangle,is_unit,x,ldx,beta,y,ldy,tra,ctra,acc) 
    integer(psb_ipk_), intent(in)    :: m,n,ldx,ldy,nxy,nc,ldj,ldv
    integer(psb_ipk_), intent(in)    :: ja(ldj,*)
    real(psb_dpk_), intent(in)      :: alpha, beta, x(ldx,*),val(ldv,*)
    real(psb_dpk_), intent(inout)   :: y(ldy,*)
    logical, intent(in)              :: is_triangle,is_unit,tra, ctra

    real(psb_dpk_), intent(inout)   :: acc(*)
    integer(psb_ipk_)   :: i,j,k, ir, jc


    if (alpha == dzero) then
      if (beta == dzero) then
        do i = 1, m
          y(i,1:nxy) = dzero
        enddo
      else
        do  i = 1, m
          y(i,1:nxy) = beta*y(i,1:nxy)
        end do
      endif
      return
    end if

    if (.not.(tra.or.ctra)) then 

      if (beta == dzero) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = acc(1:nxy)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) - val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = acc(1:nxy)
          end do

        else 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = alpha*acc(1:nxy)
          end do

        end if


      else if (beta == done) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nxy)  = y(i,1:nxy)
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = acc(1:nxy)
          end do

        else if (alpha == -done) then 

          do i=1,m
            acc(1:nxy)  = y(i,1:nxy)
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) - val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = acc(1:nxy)
          end do

        else 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = y(i,1:nxy) + alpha*acc(1:nxy)
          end do

        end if

      else if (beta == -done) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = -y(i,1:nxy) + acc(1:nxy)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = -y(i,1:nxy) -acc(1:nxy)
          end do

        else 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = -y(i,1:nxy) + alpha*acc(1:nxy)
          end do

        end if

      else 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = beta*y(i,1:nxy) + acc(1:nxy)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = beta*y(i,1:nxy) - acc(1:nxy)
          end do

        else 

          do i=1,m 
            acc(1:nxy)  = dzero
            do j=1,nc
              acc(1:nxy)  = acc(1:nxy) + val(i,j) * x(ja(i,j),1:nxy)          
            enddo
            y(i,1:nxy) = beta*y(i,1:nxy) + alpha*acc(1:nxy)
          end do

        end if

      end if

    else if (tra) then 

      if (beta == dzero) then 
        do i=1, m
          y(i,1:nxy) = dzero
        end do
      else if (beta == done) then 
        ! Do nothing
      else if (beta == -done) then 
        do i=1, m
          y(i,1:nxy) = -y(i,1:nxy) 
        end do
      else
        do i=1, m
          y(i,1:nxy) = beta*y(i,1:nxy) 
        end do
      end if

      if (alpha == done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) +  val(i,j)*x(i,1:nxy)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) -  val(i,j)*x(i,1:nxy)
          end do
        enddo

      else                    

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) + alpha*val(i,j)*x(i,1:nxy)
          end do
        enddo

      end if

    else if (ctra) then 

      if (beta == dzero) then 
        do i=1, m
          y(i,1:nxy) = dzero
        end do
      else if (beta == done) then 
        ! Do nothing
      else if (beta == -done) then 
        do i=1, m
          y(i,1:nxy) = -y(i,1:nxy) 
        end do
      else
        do i=1, m
          y(i,1:nxy) = beta*y(i,1:nxy) 
        end do
      end if

      if (alpha == done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) +  (val(i,j))*x(i,1:nxy)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) -  (val(i,j))*x(i,1:nxy)
          end do
        enddo

      else                    

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir,1:nxy) = y(ir,1:nxy) + alpha*(val(i,j))*x(i,1:nxy)
          end do
        enddo

      end if

    endif

    if (is_unit) then 
      do i=1, min(m,n)
        y(i,1:nxy) = y(i,1:nxy) + alpha*x(i,1:nxy)
      end do
    end if

  end subroutine psb_d_ell_csmm_inner
end subroutine psb_d_ell_csmm
