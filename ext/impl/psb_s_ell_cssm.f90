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
  

subroutine psb_s_ell_cssm(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_s_ell_mat_mod, psb_protect_name => psb_s_ell_cssm
  implicit none 
  class(psb_s_ell_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_spk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)       :: info
  character, optional, intent(in)      :: trans

  character :: trans_
  integer(psb_ipk_)   :: i,j,k,m,n, nnz, ir, jc, nxy
  real(psb_spk_), allocatable :: tmp(:,:), acc(:)
  logical            :: tra, ctra
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='s_ell_cssm'
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
  m = a%get_nrows()

  if (.not. (a%is_triangle().and.a%is_sorted())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(x,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3*ione,m/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5*ione,m/))
    goto 9999
  end if

  nxy = min(size(x,2),size(y,2))

  if (alpha == szero) then
    if (beta == szero) then
      do i = 1, m
        y(i,:) = szero
      enddo
    else
      do  i = 1, m
        y(i,:) = beta*y(i,:)
      end do
    endif
    return
  end if

  if (beta == szero) then 
    allocate(acc(nxy), stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    call inner_ellsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nxy,&
         & size(a%ja,2,kind=psb_ipk_),a%irn,a%idiag,&
         & a%ja,size(a%ja,1,kind=psb_ipk_),a%val,size(a%val,1,kind=psb_ipk_),&
         & x,size(x,1,kind=psb_ipk_),y,size(y,1,kind=psb_ipk_),acc,info) 

    if (info /= 0) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (alpha == sone) then 
      ! do nothing
    else if (alpha == -sone) then 
      do  i = 1, m
        y(i,:) = -y(i,:)
      end do
    else
      do  i = 1, m
        y(i,:) = alpha*y(i,:)
      end do
    end if
  else 
    allocate(tmp(m,nxy),acc(nxy), stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    call inner_ellsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nxy,&
         & size(a%ja,2,kind=psb_ipk_),a%irn,a%idiag,&
         & a%ja,size(a%ja,1,kind=psb_ipk_),a%val,size(a%val,1,kind=psb_ipk_),&
         & x,size(x,1,kind=psb_ipk_),tmp,size(tmp,1,kind=psb_ipk_),acc,info) 

    if (info == 0) &
         & call psb_geaxpby(m,nxy,alpha,tmp,beta,y,info)

    if (info /= 0) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif

  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains 

  subroutine inner_ellsm(tra,ctra,lower,unit,n,nc,nxy,irn,idiag,ja,ldj,val,ldv,&
       & x,ldx,y,ldy,acc,info) 
    implicit none 
    logical, intent(in)                  :: tra,ctra,lower,unit
    integer(psb_ipk_), intent(in)        :: n,nc,ldj,ldv,nxy,ldx,ldy
    integer(psb_ipk_), intent(in)        :: irn(*),idiag(*), ja(ldj,*)
    real(psb_spk_), intent(in)          :: val(ldv,*)
    real(psb_spk_), intent(in)          :: x(ldx,nxy)
    real(psb_spk_), intent(out)         :: y(ldy,nxy), acc(nxy)
    integer(psb_ipk_), intent(out)       :: info

    integer(psb_ipk_) :: i,j,k,m, ir, jc

    !
    ! The only error condition here is if
    ! the matrix is non-unit and some idiag value is illegal.
    !
    info = 0 

    if (.not.(tra.or.ctra)) then 

      if (lower) then 

        if (unit) then 
          do i=1, n
            acc = szero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = szero 
            do j=1,idiag(i)-1
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = (x(i,:) - acc)/val(i,idiag(i))
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 

          do i=n, 1, -1 
            acc = szero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do

        else if (.not.unit) then 

          do i=n, 1, -1 
            acc = szero 
            do j=idiag(i)+1, irn(i)
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = (x(i,:) - acc)/val(i,idiag(i))
          end do

        end if

      end if

    else if (tra) then 

      do i=1, n
        y(i,:) = x(i,:)
      end do

      if (lower) then 

        if (unit) then 

          do i=n, 1, -1
            acc = y(i,:) 
            do j=1,irn(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - val(i,j)*acc 
            end do

          end do

        else if (.not.unit) then 

          do i=n, 1, -1
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = y(i,:)/val(i,idiag(i))
            acc  = y(i,:) 
            do j=1,idiag(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - val(i,j)*acc 
            end do
          end do

        end if

      else if (.not.lower) then 

        if (unit) then 

          do i=1, n
            acc  = y(i,:) 
            do j=1, irn(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - val(i,j)*acc 
            end do
          end do

        else if (.not.unit) then 

          do i=1, n
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = y(i,:)/val(i,idiag(i))
            acc  = y(i,:) 
            do j=idiag(i)+1, irn(i) 
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - val(i,j)*acc 
            end do
          end do

        end if

      end if

    else if (ctra) then 

      do i=1, n
        y(i,:) = x(i,:)
      end do

      if (lower) then 

        if (unit) then 

          do i=n, 1, -1
            acc = y(i,:) 
            do j=1,irn(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - (val(i,j))*acc 
            end do

          end do

        else if (.not.unit) then 

          do i=n, 1, -1
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = y(i,:)/(val(i,idiag(i)))
            acc  = y(i,:) 
            do j=1,idiag(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - (val(i,j))*acc 
            end do
          end do

        end if

      else if (.not.lower) then 

        if (unit) then 

          do i=1, n
            acc  = y(i,:) 
            do j=1, irn(i)
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - (val(i,j))*acc 
            end do
          end do

        else if (.not.unit) then 

          do i=1, n
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i,:) = y(i,:)/(val(i,idiag(i)))
            acc  = y(i,:) 
            do j=idiag(i)+1, irn(i) 
              jc    = ja(i,j)
              y(jc,:) = y(jc,:) - (val(i,j))*acc 
            end do
          end do

        end if

      end if
    end if
  end subroutine inner_ellsm
end subroutine psb_s_ell_cssm
