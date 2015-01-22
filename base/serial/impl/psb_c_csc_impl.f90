!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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

! == ===================================
!
!
!
! Computational routines
!
!
!
!
!
!
! == ===================================

subroutine psb_c_csc_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csmv
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc
  complex(psb_spk_) :: acc
  logical   :: tra
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

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


  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')

  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if


  if (size(x,1)<n) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = n; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if


  if (alpha == czero) then
    if (beta == czero) then
      do i = 1, m
        y(i) = czero
      enddo
    else
      do  i = 1, m
        y(i) = beta*y(i)
      end do
    endif
    return
  end if

  if (tra) then 

    if (beta == czero) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = alpha*acc
        end do

      end if


    else if (beta == cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = y(i) + alpha*acc
        end do

      end if

    else if (beta == -cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = -y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = -y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = -y(i) + alpha*acc
        end do

      end if

    else 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) - acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) + alpha*acc
        end do

      end if

    end if

  else if (.not.tra) then 

    if (beta == czero) then 
      do i=1, m
        y(i) = czero
      end do
    else if (beta == cone) then 
      ! Do nothing
    else if (beta == -cone) then 
      do i=1, m
        y(i) = -y(i) 
      end do
    else
      do i=1, m
        y(i) = beta*y(i) 
      end do
    end if

    if (alpha == cone) then

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir) = y(ir) +  a%val(j)*x(i)
        end do
      enddo

    else if (alpha == -cone) then

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir) = y(ir) -  a%val(j)*x(i)
        end do
      enddo

    else                    

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir) = y(ir) + alpha*a%val(j)*x(i)
        end do
      enddo

    end if

  endif

  if (a%is_unit()) then 
    do i=1, min(m,n)
      y(i) = y(i) + alpha*x(i)
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_csmv

subroutine psb_c_csc_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csmm
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_), allocatable  :: acc(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if

  if (size(x,1)<n) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = n; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if


  nc = min(size(x,2) , size(y,2) )

  allocate(acc(nc), stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='allocate')
    goto 9999
  end if

  if (alpha == czero) then
    if (beta == czero) then
      do i = 1, m
        y(i,:) = czero
      enddo
    else
      do  i = 1, m
        y(i,:) = beta*y(i,:)
      end do
    endif
    return
  end if

  if (tra) then 

    if (beta == czero) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = alpha*acc
        end do

      end if


    else if (beta == cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) + alpha*acc
        end do

      end if

    else if (beta == -cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) + alpha*acc
        end do

      end if

    else 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) - acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + a%val(j) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + alpha*acc
        end do

      end if

    end if

  else if (.not.tra) then 

    if (beta == czero) then 
      do i=1, m
        y(i,:) = czero
      end do
    else if (beta == cone) then 
      ! Do nothing
    else if (beta == -cone) then 
      do i=1, m
        y(i,:) = -y(i,:) 
      end do
    else
      do i=1, m
        y(i,:) = beta*y(i,:) 
      end do
    end if

    if (alpha == cone) then

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir,:) = y(ir,:) +  a%val(j)*x(i,:)
        end do
      enddo

    else if (alpha == -cone) then

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir,:) = y(ir,:) -  a%val(j)*x(i,:)
        end do
      enddo

    else                    

      do i=1,n
        do j=a%icp(i), a%icp(i+1)-1
          ir = a%ia(j)
          y(ir,:) = y(ir,:) + alpha*a%val(j)*x(i,:)
        end do
      enddo

    end if

  endif

  if (a%is_unit()) then 
    do i=1, min(m,n)
      y(i,:) = y(i,:) + alpha*x(i,:)
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_csmm


subroutine psb_c_csc_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_cssv
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m, nnz, ir, jc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: tmp(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_cssv'
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

  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
  m = a%get_nrows()

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(x,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if


  
  if (alpha == czero) then
    if (beta == czero) then
      do i = 1, m
        y(i) = czero
      enddo
    else
      do  i = 1, m
        y(i) = beta*y(i)
      end do
    endif
    return
  end if

  if (beta == czero) then 
    call inner_cscsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & a%icp,a%ia,a%val,x,y) 
    if (alpha == cone) then 
      ! do nothing
    else if (alpha == -cone) then 
      do  i = 1, m
        y(i) = -y(i)
      end do
    else
      do  i = 1, m
        y(i) = alpha*y(i)
      end do
    end if
  else 
    allocate(tmp(m), stat=info) 
    if (info /= psb_success_) then 
      return
    end if
    tmp(1:m) = x(1:m)
    call inner_cscsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & a%icp,a%ia,a%val,tmp,y) 
    do  i = 1, m
      y(i) = alpha*tmp(i) + beta*y(i)
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains 

  subroutine inner_cscsv(tra,lower,unit,n,icp,ia,val,x,y) 
    implicit none 
    logical, intent(in)                 :: tra,lower,unit  
    integer(psb_ipk_), intent(in)                 :: icp(*), ia(*),n
    complex(psb_spk_), intent(in)          :: val(*)
    complex(psb_spk_), intent(in)          :: x(*)
    complex(psb_spk_), intent(out)         :: y(*)

    integer(psb_ipk_) :: i,j,k,m, ir, jc
    complex(psb_spk_) :: acc

    if (tra) then 

      if (lower) then 
        if (unit) then 
          do i=n, 1, -1 
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + val(j)*y(ia(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=n, 1, -1 
            acc = czero 
            do j=icp(i)+1, icp(i+1)-1
              acc = acc + val(j)*y(ia(j))
            end do
            y(i) = (x(i) - acc)/val(icp(i))
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 
          do i=1, n
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + val(j)*y(ia(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = czero 
            do j=icp(i), icp(i+1)-2
              acc = acc + val(j)*y(ia(j))
            end do
            y(i) = (x(i) - acc)/val(icp(i+1)-1)
          end do
        end if

      end if

    else if (.not.tra) then 

      do i=1, n
        y(i) = x(i)
      end do

      if (lower) then 

        if (unit) then 
          do i=1, n
            acc  = y(i) 
            do j=icp(i), icp(i+1)-1
              jc    = ia(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, n
            y(i) = y(i)/val(icp(i))
            acc  = y(i) 
            do j=icp(i)+1, icp(i+1)-1
              jc    = ia(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 
          do i=n, 1, -1
            acc = y(i) 
            do j=icp(i), icp(i+1)-1
              jc    = ia(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=n, 1, -1
            y(i) = y(i)/val(icp(i+1)-1)
            acc  = y(i) 
            do j=icp(i), icp(i+1)-2
              jc    = ia(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_cscsv

end subroutine psb_c_csc_cssv



subroutine psb_c_csc_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_cssm
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: tmp(:,:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_base_csmm'
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


  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')
  m   = a%get_nrows()

  if (size(x,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = psb_err_input_asize_small_i_
    ierr(1) = 5; ierr(2) = m; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  nc  = min(size(x,2) , size(y,2)) 

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if


  if (alpha == czero) then
    if (beta == czero) then
      do i = 1, m
        y(i,:) = czero
      enddo
    else
      do  i = 1, m
        y(i,:) = beta*y(i,:)
      end do
    endif
    return
  end if

  if (beta == czero) then 
    call inner_cscsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%icp,a%ia,a%val,x,size(x,1,kind=psb_ipk_),y,size(y,1,kind=psb_ipk_),info) 
    do  i = 1, m
      y(i,1:nc) = alpha*y(i,1:nc)
    end do
  else 
    allocate(tmp(m,nc), stat=info) 
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='allocate')
      goto 9999
    end if

    tmp(1:m,:) = x(1:m,1:nc)
    call inner_cscsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%icp,a%ia,a%val,tmp,size(tmp,1,kind=psb_ipk_),y,size(y,1,kind=psb_ipk_),info) 
    do  i = 1, m
      y(i,1:nc) = alpha*tmp(i,1:nc) + beta*y(i,1:nc)
    end do
  end if

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='inner_cscsm')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return


contains 

  subroutine inner_cscsm(tra,lower,unit,nr,nc,&
       & icp,ia,val,x,ldx,y,ldy,info) 
    implicit none 
    logical, intent(in)                 :: tra,lower,unit
    integer(psb_ipk_), intent(in)                 :: nr,nc,ldx,ldy,icp(*),ia(*)
    complex(psb_spk_), intent(in)          :: val(*), x(ldx,*)
    complex(psb_spk_), intent(out)         :: y(ldy,*)
    integer(psb_ipk_), intent(out)                :: info
    integer(psb_ipk_) :: i,j,k,m, ir, jc
    complex(psb_spk_), allocatable  :: acc(:)

    info = psb_success_
    allocate(acc(nc), stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      return
    end if


    if (tra) then 

      if (lower) then 
        if (unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + val(j)*y(ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=icp(i)+1, icp(i+1)-1
              acc = acc + val(j)*y(ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/val(icp(i))
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 
          do i=1, nr
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + val(j)*y(ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=1, nr
            acc = czero 
            do j=icp(i), icp(i+1)-2
              acc = acc + val(j)*y(ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/val(icp(i+1)-1)
          end do
        end if

      end if

    else if (.not.tra) then 

      do i=1, nr
        y(i,1:nc) = x(i,1:nc)
      end do

      if (lower) then

        if (unit) then  
          do i=1, nr
            acc  = y(i,1:nc) 
            do j=icp(i), icp(i+1)-1
              jc    = ia(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, nr
            y(i,1:nc) = y(i,1:nc)/val(icp(i))
            acc    = y(i,1:nc) 
            do j=icp(i)+1, icp(i+1)-1
              jc      = ia(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 
          do i=nr, 1, -1
            acc = y(i,1:nc) 
            do j=icp(i), icp(i+1)-1
              jc    = ia(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1
            y(i,1:nc) = y(i,1:nc)/val(icp(i+1)-1)
            acc  = y(i,1:nc) 
            do j=icp(i), icp(i+1)-2
              jc    = ia(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_cscsm

end subroutine psb_c_csc_cssm

function psb_c_csc_maxval(a) result(res)
  use psb_error_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_maxval
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  real(psb_spk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_maxval'
  logical, parameter :: debug=.false.


  if (a%is_unit()) then 
    res = sone
  else
    res = szero
  end if

  nnz = a%get_nzeros()
  if (allocated(a%val)) then 
    nnz = min(nnz,size(a%val))
    res = maxval(abs(a%val(1:nnz)))
  end if
end function psb_c_csc_maxval

function psb_c_csc_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csnm1

  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  real(psb_spk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_spk_) :: acc
  real(psb_spk_), allocatable :: vt(:)
  logical            :: tra, is_unit
  integer(psb_ipk_)  :: err_act
  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name='c_csc_csnm1'
  logical, parameter :: debug=.false.


  res = szero 
  m = a%get_nrows()
  n = a%get_ncols()
  is_unit = a%is_unit()
  do j=1, n
    if (is_unit) then 
      acc = sone
    else
      acc = szero
    end if
    do k=a%icp(j),a%icp(j+1)-1
      acc = acc + abs(a%val(k))
    end do
    res = max(res,acc)
  end do
  
  return

end function psb_c_csc_csnm1

subroutine psb_c_csc_colsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_colsum
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(out)             :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: vt(:)
  logical            :: tra, is_unit
  integer(psb_ipk_)  :: err_act, info
  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_ncols()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    ierr(1) = 1; ierr(2) = size(d); ierr(3) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  is_unit = a%is_unit()
  do i = 1, a%get_ncols()
    if (is_unit) then 
      d(i) = cone
    else
      d(i) = czero
    end if

    do j=a%icp(i),a%icp(i+1)-1  
      d(i) = d(i) + (a%val(j))
    end do
  end do
  
  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_colsum

subroutine psb_c_csc_aclsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_aclsum
  class(psb_c_csc_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_spk_) :: acc
  real(psb_spk_), allocatable :: vt(:)
  logical            :: tra, is_unit
  integer(psb_ipk_)  :: err_act, info
  integer(psb_ipk_)  :: ierr(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_ncols()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    ierr(1) = 1; ierr(2) = size(d); ierr(3) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  is_unit = a%is_unit()
  do i = 1, a%get_ncols()
    if (is_unit) then 
      d(i) = sone
    else
      d(i) = szero
    end if

    do j=a%icp(i),a%icp(i+1)-1  
      d(i) = d(i) + abs(a%val(j))
    end do
  end do
  
  if (a%is_unit()) then 
    do i=1, a%get_ncols()
      d(i) = d(i) + sone
    end do
  end if

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_aclsum

subroutine psb_c_csc_rowsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_rowsum
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_ncols()
  n = a%get_nrows()
  if (size(d) < n) then 
    info=psb_err_input_asize_small_i_
    ierr(1) = 1; ierr(2) = size(d); ierr(3) = n
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (a%is_unit()) then 
    d = cone
  else
    d = czero
  end if

  do i=1, m
    do j=a%icp(i),a%icp(i+1)-1
      k = a%ia(j)
      d(k) = d(k) + (a%val(k))
    end do
  end do

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_rowsum

subroutine psb_c_csc_arwsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_arwsum
  class(psb_c_csc_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_spk_) :: acc
  real(psb_spk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_ncols()
  n = a%get_nrows()
  if (size(d) < n) then 
    info=psb_err_input_asize_small_i_
    ierr(1) = 1; ierr(2) = size(d); ierr(3) = n
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (a%is_unit()) then 
    d = sone
  else
    d = szero
  end if

  do i=1, m
    do j=a%icp(i),a%icp(i+1)-1
      k = a%ia(j)
      d(k) = d(k) + abs(a%val(k))
    end do
  end do

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_arwsum


subroutine psb_c_csc_get_diag(a,d,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_get_diag
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act, mnm, i, j, k
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  mnm = min(a%get_nrows(),a%get_ncols())
  if (size(d) < mnm) then 
    info=psb_err_input_asize_invalid_i_
    ierr(1) = 2; ierr(2) = size(d); 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if


  if (a%is_unit()) then 
    d(1:mnm) = cone 
  else
    do i=1, mnm
      d(i) = czero
      do k=a%icp(i),a%icp(i+1)-1
        j=a%ia(k)
        if ((j == i) .and.(j <= mnm )) then 
          d(i) = a%val(k)
        endif
      enddo
    end do
  endif
  do i=mnm+1,size(d) 
    d(i) = czero
  end do
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_get_diag


subroutine psb_c_csc_scal(d,a,info,side) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_scal
  use psb_string_mod
  implicit none 
  class(psb_c_csc_sparse_mat), intent(inout) :: a
  complex(psb_spk_), intent(in)                 :: d(:)
  integer(psb_ipk_), intent(out)             :: info
  character, intent(in), optional            :: side

  integer(psb_ipk_) :: err_act,mnm, i, j, n
  type(psb_c_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name='scal'
  character :: side_
  logical   :: left 
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  side_ = 'L'
  if (present(side)) then 
    side_ = psb_toupper(side)
  end if

  if (a%is_unit()) then 
    call a%make_nonunit()
  end if

  left = (side_ == 'L')
  
  if (left) then 
    n = a%get_ncols()
    if (size(d) < n) then 
      info=psb_err_input_asize_invalid_i_
      ierr(1) = 2; ierr(2) = size(d); 
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if
    
    do i=1, a%get_nzeros()
      a%val(i) = a%val(i) * d(a%ia(i))
    enddo
  else
    n = a%get_nrows()
    if (size(d) < n) then 
      info=psb_err_input_asize_invalid_i_
      ierr(1) = 2; ierr(2) = size(d); 
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if

    do j=1, n
      do i = a%icp(j), a%icp(j+1) -1 
        a%val(i) = a%val(i) * d(j)
      end do
    enddo
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_scal


subroutine psb_c_csc_scals(d,a,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_scals
  implicit none 
  class(psb_c_csc_sparse_mat), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act,mnm, i, j, m
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_unit()) then 
    call a%make_nonunit()
  end if

  do i=1,a%get_nzeros()
    a%val(i) = a%val(i) * d
  enddo

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_scals


! == =================================== 
!
!
!
! Data management
!
!
!
!
!
! == ===================================   

subroutine psb_c_csc_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csgetptn
  implicit none

  class(psb_c_csc_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  nz = 0

  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif

  if ((imax<imin).or.(jmax_<jmin_)) return

  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif
  if ((append_).and.(present(nzin))) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .false.
  endif
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .false.
  endif
  if ((rscale_.or.cscale_).and.(present(iren))) then 
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='iren (rscale.or.cscale)')
    goto 9999
  end if

  call csc_getptn(imin,imax,jmin_,jmax_,a,nz,ia,ja,nzin_,append_,info,iren)
  
  if (rscale_) then 
    do i=nzin_+1, nzin_+nz
      ia(i) = ia(i) - imin + 1
    end do
  end if
  if (cscale_) then 
    do i=nzin_+1, nzin_+nz
      ja(i) = ja(i) - jmin_ + 1
    end do
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  subroutine csc_getptn(imin,imax,jmin,jmax,a,nz,ia,ja,nzin,append,info,&
       & iren)

    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none

    class(psb_c_csc_sparse_mat), intent(in)    :: a
    integer(psb_ipk_) :: imin,imax,jmin,jmax
    integer(psb_ipk_), intent(inout)               :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    integer(psb_ipk_), intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional                    :: iren(:)
    integer(psb_ipk_) :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,m,isz, nrd, ncd
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20) :: name='coo_getrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%get_nzeros()
    irw = imin
    lrw = min(imax,a%get_nrows())
    icl = jmin
    lcl = min(jmax,a%get_ncols())
    if (irw<0) then 
      info = psb_err_pivot_too_small_
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    nrd = max(1,a%get_nrows())
    ncd = max(1,a%get_ncols())
    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
         & ((nza+ncd-1)/ncd)*(lcl+1-icl),&
         & ((nza+nrd-1)/nrd)*(lrw+1-irw))
    nz = 0 

    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)

    if (info /= psb_success_) return
    isz = min(size(ia),size(ja))
    if (present(iren)) then 
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+ione,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,ja,info)
              isz = min(size(ia),size(ja))
            end if
            nz    = nz + 1
            ia(nzin_)  = iren(a%ia(j))
            ja(nzin_)  = iren(i)
          end if
        enddo
      end do
    else
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+ione,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,ja,info)
              isz = min(size(ia),size(ja))
            end if
            nz    = nz + 1
            ia(nzin_)  = (a%ia(j))
            ja(nzin_)  = (i)
          end if
        enddo
      end do
    end if
    
  end subroutine csc_getptn
  
end subroutine psb_c_csc_csgetptn




subroutine psb_c_csc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csgetrow
  implicit none

  class(psb_c_csc_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  nz = 0
    
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif


  if ((imax<imin).or.(jmax_<jmin_)) return

  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif
  if ((append_).and.(present(nzin))) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .false.
  endif
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .false.
  endif
  if ((rscale_.or.cscale_).and.(present(iren))) then 
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='iren (rscale.or.cscale)')
    goto 9999
  end if

  call csc_getrow(imin,imax,jmin_,jmax_,a,nz,ia,ja,val,nzin_,append_,info,&
       & iren)
  
  if (rscale_) then 
    do i=nzin_+1, nzin_+nz
      ia(i) = ia(i) - imin + 1
    end do
  end if
  if (cscale_) then 
    do i=nzin_+1, nzin_+nz
      ja(i) = ja(i) - jmin_ + 1
    end do
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  subroutine csc_getrow(imin,imax,jmin,jmax,a,nz,ia,ja,val,nzin,append,info,&
       & iren)

    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none

    class(psb_c_csc_sparse_mat), intent(in)    :: a
    integer(psb_ipk_) :: imin,imax,jmin,jmax
    integer(psb_ipk_), intent(inout)               :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer(psb_ipk_), intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional                    :: iren(:)
    integer(psb_ipk_) :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,isz,m, nrd, ncd
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20) :: name='coo_getrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    m = a%get_nrows()
    nza = a%get_nzeros()
    irw = imin
    lrw = min(imax,a%get_nrows())
    icl = jmin
    lcl = min(jmax,a%get_ncols())
    if (irw<0) then 
      info = psb_err_pivot_too_small_
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    nrd = max(1,a%get_nrows())
    ncd = max(1,a%get_ncols())
    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
         & ((nza+ncd-1)/ncd)*(lcl+1-icl),&
         & ((nza+nrd-1)/nrd)*(lrw+1-irw))
    nz = 0 

    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,val,info)

    if (info /= psb_success_) return
    isz = min(size(ia),size(ja),size(val))
    if (present(iren)) then 
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+ione,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,ja,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,val,info)
              isz = min(size(ia),size(ja),size(val))
            end if
            nz    = nz + 1
            val(nzin_) = a%val(j)
            ia(nzin_)  = iren(a%ia(j))
            ja(nzin_)  = iren(i)
          end if
        enddo
      end do
    else
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+ione,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,ja,info)
              call psb_ensure_size(int(1.25*nzin_)+ione,val,info)
              isz = min(size(ia),size(ja),size(val))
            end if
            nz    = nz + 1
            val(nzin_) = a%val(j)
            ia(nzin_)  = (a%ia(j))
            ja(nzin_)  = (i)
          end if
        enddo
      end do
    end if
  end subroutine csc_getrow

end subroutine psb_c_csc_csgetrow



subroutine psb_c_csc_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csput_a
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)


  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_csput_a'
  logical, parameter :: debug=.false.
  integer(psb_ipk_) :: nza, i,j,k, nzl, isza, debug_level, debug_unit

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info = psb_success_

  if (nz <= 0) then 
    info = psb_err_iarg_neg_
    ierr(1)=1
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = psb_err_input_asize_invalid_i_
    ierr(1)=2
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = psb_err_input_asize_invalid_i_
    ierr(1)=3
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(val) < nz) then 
    info = psb_err_input_asize_invalid_i_
    ierr(1)=4
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (nz == 0) return

  nza  = a%get_nzeros()

  if (a%is_bld()) then 
    ! Build phase should only ever be in COO
    info = psb_err_invalid_mat_state_

  else  if (a%is_upd()) then 
    call  psb_c_csc_srch_upd(nz,ia,ja,val,a,&
         & imin,imax,jmin,jmax,info,gtl)

    if (info < 0) then 
      info = psb_err_internal_error_
    else if (info > 0) then 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & ': Discarded entries not  belonging to us.'                    
      info = psb_success_
    end if

  else 
    ! State is wrong.
    info = psb_err_invalid_mat_state_
  end if
  if (info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


contains

  subroutine psb_c_csc_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_sort_mod
    implicit none 

    class(psb_c_csc_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: nz, imin,imax,jmin,jmax
    integer(psb_ipk_), intent(in) :: ia(:),ja(:)
    complex(psb_spk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in), optional  :: gtl(:)
    integer(psb_ipk_) :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nr,nc,nnz,dupl,ng, nar, nac
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20)    :: name='c_csc_srch_upd'

    info = psb_success_
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    dupl = a%get_dupl()

    if (.not.a%is_sorted()) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = a%get_nzeros()
    nar = a%get_nrows()
    nac = a%get_ncols()

    if (present(gtl)) then 
      ng = size(gtl)

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            ic = gtl(ic)
            if ((ic > 0).and.(ic <= nac)) then 
              i1 = a%icp(ic)
              i2 = a%icp(ic+1)
              nr=i2-i1

              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
              if (ip>0) then 
                a%val(i1+ip-1) = val(i)
              else
                info = max(info,3)
              end if
            else
              info = max(info,2)
            end if
          else
            info = max(info,1)
          end if
        end do

      case(psb_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            ic = gtl(ic)
            if ((ic > 0).and.(ic <= nac)) then 
              i1 = a%icp(ic)
              i2 = a%icp(ic+1)
              nr=i2-i1

              ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
              if (ip>0) then 
                a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
              else
                info = max(info,3)
              end if
            else
              info = max(info,2)
            end if
          else
            info = max(info,1)
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    else

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 

          if ((ic > 0).and.(ic <= nac)) then 
            i1 = a%icp(ic)
            i2 = a%icp(ic+1)
            nr=i2-i1

            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
            if (ip>0) then 
              a%val(i1+ip-1) = val(i)
            else
              info = max(info,3)
            end if
          else
            info = max(info,2)
          end if

        end do

      case(psb_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ic > 0).and.(ic <= nac)) then 
            i1 = a%icp(ic)
            i2 = a%icp(ic+1)
            nr=i2-i1

            ip = psb_ibsrch(ir,nr,a%ia(i1:i2-1))    
            if (ip>0) then 
              a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
            else
              info = max(info,3)
            end if
          else
              info = max(info,2)
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine psb_c_csc_srch_upd

end subroutine psb_c_csc_csput_a



subroutine psb_c_cp_csc_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_cp_csc_from_coo
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)                        :: info

  type(psb_c_coo_sparse_mat)   :: tmp
  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  ! We need to make a copy because mv_from will have to
  ! sort in column-major order.
  call tmp%cp_from_coo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

end subroutine psb_c_cp_csc_from_coo



subroutine psb_c_cp_csc_to_coo(a,b,info) 
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_cp_csc_to_coo
  implicit none 

  class(psb_c_csc_sparse_mat), intent(in)  :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                      :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,irw, err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat

  do i=1, nc
    do j=a%icp(i),a%icp(i+1)-1
      b%ia(j)  = a%ia(j)
      b%ja(j)  = i
      b%val(j) = a%val(j)
    end do
  end do
  
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)


end subroutine psb_c_cp_csc_to_coo


subroutine psb_c_mv_csc_to_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_mv_csc_to_coo
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                        :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,irw, err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat
  call b%set_nzeros(a%get_nzeros())
  call move_alloc(a%ia,b%ia)
  call move_alloc(a%val,b%val)
  call psb_realloc(nza,b%ja,info)
  if (info /= psb_success_) return
  do i=1, nc
    do j=a%icp(i),a%icp(i+1)-1
      b%ja(j)  = i
    end do
  end do
  call a%free()
  call b%fix(info)

end subroutine psb_c_mv_csc_to_coo



subroutine psb_c_mv_csc_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_mv_csc_from_coo
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                        :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, err_act, nc, icl
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name='c_mv_csc_from_coo'

  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  call b%fix(info, idir=psb_col_major_)
  if (info /= psb_success_) return

  nr  = b%get_nrows()
  nc  = b%get_ncols()
  nza = b%get_nzeros()
  
  a%psb_c_base_sparse_mat = b%psb_c_base_sparse_mat

  ! Dirty trick: call move_alloc to have the new data allocated just once.
  call move_alloc(b%ja,itemp)
  call move_alloc(b%ia,a%ia)
  call move_alloc(b%val,a%val)
  call psb_realloc(max(nr+1,nc+1),a%icp,info)
  call b%free()

  if (nza <= 0) then 
    a%icp(:) = 1
  else
    a%icp(1) = 1
    if (nc < itemp(nza)) then 
      write(debug_unit,*) trim(name),': CLSHR=.false. : ',&
           &nc,itemp(nza),' Expect trouble!'
      info = 12
    end if

    j = 1 
    i = 1
    icl = itemp(j) 

    outer: do 
      inner: do 
        if (i >= icl) exit inner
        if (i > nc) then 
          write(debug_unit,*) trim(name),&
               & 'Strange situation: i>nr ',i,nc,j,nza,icl
          exit outer
        end if
        a%icp(i+1) = a%icp(i) 
        i = i + 1
      end do inner
      j = j + 1
      if (j > nza) exit
      if (itemp(j) /= icl) then 
        a%icp(i+1) = j
        icl = itemp(j) 
        i = i + 1
      endif
      if (i > nc) exit
    enddo outer
    !
    ! Cleanup empty rows at the end
    !
    if (j /= (nza+1)) then 
      write(debug_unit,*) trim(name),': Problem from loop :',j,nza
      info = 13
    endif
    do 
      if (i > nc) exit
      a%icp(i+1) = j
      i = i + 1
    end do

  endif


end subroutine psb_c_mv_csc_from_coo


subroutine psb_c_mv_csc_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_mv_csc_to_fmt
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout)  :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
    ! Need to fix trivial copies! 
  type is (psb_c_csc_sparse_mat) 
    b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat
    call move_alloc(a%icp, b%icp)
    call move_alloc(a%ia,  b%ia)
    call move_alloc(a%val, b%val)
    call a%free()
    
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_c_mv_csc_to_fmt
!!$

subroutine psb_c_cp_csc_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_cp_csc_to_fmt
  implicit none 

  class(psb_c_csc_sparse_mat), intent(in)   :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                       :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nz, nr, i,j,irw, err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_


  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    b%psb_c_base_sparse_mat = a%psb_c_base_sparse_mat
    nc = a%get_ncols()
    nz = a%get_nzeros()
    if (info == 0) call psb_safe_cpy( a%icp(1:nc+1), b%icp , info)
    if (info == 0) call psb_safe_cpy( a%ia(1:nz),    b%ia  , info)
    if (info == 0) call psb_safe_cpy( a%val(1:nz),   b%val , info)

  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_c_cp_csc_to_fmt


subroutine psb_c_mv_csc_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_mv_csc_from_fmt
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout)  :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                         :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    a%psb_c_base_sparse_mat = b%psb_c_base_sparse_mat
    call move_alloc(b%icp, a%icp)
    call move_alloc(b%ia,  a%ia)
    call move_alloc(b%val, a%val)
    call b%free()

  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

end subroutine psb_c_mv_csc_from_fmt



subroutine psb_c_cp_csc_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_cp_csc_from_fmt
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(in)   :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nz, nr, i,j,irw, err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    a%psb_c_base_sparse_mat = b%psb_c_base_sparse_mat
    nc = b%get_ncols()
    nz = b%get_nzeros()
    if (info == 0) call psb_safe_cpy( b%icp(1:nc+1), a%icp , info)
    if (info == 0) call psb_safe_cpy( b%ia(1:nz),    a%ia  , info)
    if (info == 0) call psb_safe_cpy( b%val(1:nz),   a%val , info)

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
  
end subroutine psb_c_cp_csc_from_fmt


subroutine psb_c_csc_mold(a,b,info) 
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_mold
  use psb_error_mod
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in)                  :: a
  class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csc_mold'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  info = 0 
  if (allocated(b)) then 
    call b%free()
    deallocate(b,stat=info)
  end if
  if (info == 0) allocate(psb_c_csc_sparse_mat :: b, stat=info)

  if (info /= 0) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_mold

subroutine  psb_c_csc_reallocate_nz(nz,a) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_reallocate_nz
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_c_csc_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  call psb_realloc(nz,a%ia,info)
  if (info == psb_success_) call psb_realloc(nz,a%val,info)
  if (info == psb_success_) call psb_realloc(max(nz,a%get_nrows()+1,a%get_ncols()+1),a%icp,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_reallocate_nz



subroutine psb_c_csc_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_csgetblk
  implicit none

  class(psb_c_csc_sparse_mat), intent(in) :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act, nzin, nzout
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical :: append_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(append)) then 
    append_ = append
  else
    append_ = .false.
  endif
  if (append_) then 
    nzin = a%get_nzeros()
  else
    nzin = 0
  endif

  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
       & nzin=nzin, rscale=rscale, cscale=cscale)

  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_csgetblk

subroutine psb_c_csc_reinit(a,clear)
  use psb_error_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_reinit
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='reinit'
  logical  :: clear_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_


  if (present(clear)) then 
    clear_ = clear
  else
    clear_ = .true.
  end if

  if (a%is_bld() .or. a%is_upd()) then 
    ! do nothing
    return
  else if (a%is_asb()) then 
    if (clear_) a%val(:) = czero
    call a%set_upd()
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_reinit

subroutine  psb_c_csc_trim(a)
  use psb_realloc_mod
  use psb_error_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_trim
  implicit none 
  class(psb_c_csc_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info, nz, n
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  n   = a%get_ncols()
  nz  = a%get_nzeros()
  if (info == psb_success_) call psb_realloc(n+1,a%icp,info)
  if (info == psb_success_) call psb_realloc(nz,a%ia,info)
  if (info == psb_success_) call psb_realloc(nz,a%val,info)

  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_trim

subroutine  psb_c_csc_allocate_mnnz(m,n,a,nz) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_allocate_mnnz
  implicit none 
  integer(psb_ipk_), intent(in) :: m,n
  class(psb_c_csc_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(in), optional :: nz
  integer(psb_ipk_) :: err_act, info, nz_
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (m < 0) then 
    info = psb_err_iarg_neg_
    ierr(1) = ione; ierr(2) = izero; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif
  if (n < 0) then 
    info = psb_err_iarg_neg_
    ierr(1) = 2; ierr(2) = izero; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif
  if (present(nz)) then 
    nz_ = nz
  else
    nz_ = max(7*m,7*n,1)
  end if
  if (nz_ < 0) then 
    info = psb_err_iarg_neg_
    ierr(1) = 3; ierr(2) = izero; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif

  if (info == psb_success_) call psb_realloc(n+1,a%icp,info)
  if (info == psb_success_) call psb_realloc(nz_,a%ia,info)
  if (info == psb_success_) call psb_realloc(nz_,a%val,info)
  if (info == psb_success_) then 
    a%icp=0
    call a%set_nrows(m)
    call a%set_ncols(n)
    call a%set_bld()
    call a%set_triangle(.false.)
    call a%set_unit(.false.)
    call a%set_dupl(psb_dupl_def_)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csc_allocate_mnnz

subroutine psb_c_csc_print(iout,a,iv,head,ivr,ivc)
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => psb_c_csc_print
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_c_csc_sparse_mat), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='c_csc_print'
  logical, parameter :: debug=.false.

  character(len=*), parameter  :: datatype='complex'
  character(len=80)                 :: frmtv 
  integer(psb_ipk_) :: i,j, nmx, ni, nr, nc, nz

  if (present(head)) then 
    write(iout,'(a)') '%%MatrixMarket matrix coordinate complex general'
    write(iout,'(a,a)') '% ',head 
    write(iout,'(a)') '%'    
    write(iout,'(a,a)') '% COO'
  endif

  nr = a%get_nrows()
  nc = a%get_ncols()
  nz = a%get_nzeros()
  nmx = max(nr,nc,1)
  ni  = floor(log10(1.0*nmx)) + 1

  if (datatype=='real') then 
    write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
  else 
    write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),2(es26.18,1x),2(i',ni,',1x))'
  end if
  write(iout,*) nr, nc, nz 
  if(present(iv)) then 
    do i=1, nc
      do j=a%icp(i),a%icp(i+1)-1 
        write(iout,frmtv) iv(a%ia(j)),iv(i),a%val(j)
      end do
    enddo
  else      
    if (present(ivr).and..not.present(ivc)) then 
      do i=1, nc
        do j=a%icp(i),a%icp(i+1)-1 
          write(iout,frmtv) ivr(a%ia(j)),i,a%val(j)
        end do
      enddo
    else if (present(ivr).and.present(ivc)) then 
      do i=1, nc
        do j=a%icp(i),a%icp(i+1)-1 
          write(iout,frmtv) ivr(a%ia(j)),ivc(i),a%val(j)
        end do
      enddo
    else if (.not.present(ivr).and.present(ivc)) then 
      do i=1, nc
        do j=a%icp(i),a%icp(i+1)-1 
          write(iout,frmtv) (a%ia(j)),ivc(i),a%val(j)
        end do
      enddo
    else if (.not.present(ivr).and..not.present(ivc)) then 
      do i=1, nc
        do j=a%icp(i),a%icp(i+1)-1 
          write(iout,frmtv) (a%ia(j)),(i),a%val(j)
        end do
      enddo
    endif
  endif

end subroutine psb_c_csc_print

subroutine psb_ccscspspmm(a,b,c,info)
  use psb_c_mat_mod
  use psb_serial_mod, psb_protect_name => psb_ccscspspmm

  implicit none 

  class(psb_c_csc_sparse_mat), intent(in) :: a,b
  type(psb_c_csc_sparse_mat), intent(out)  :: c
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_) :: nze, ma,na,mb,nb, nzc, nza, nzb,nzeb
  character(len=20) :: name
  integer(psb_ipk_) :: err_act
  name='psb_cscspspmm'
  call psb_erractionsave(err_act)
  info = psb_success_

  ma = a%get_nrows()
  na = a%get_ncols()
  mb = b%get_nrows()
  nb = b%get_ncols()


  if ( mb /= na ) then 
    write(psb_err_unit,*) 'Mismatch in SPSPMM: ',ma,na,mb,nb
  endif
  nza = a%get_nzeros()
  nzb = b%get_nzeros()
  nzc = 2*(nza+nzb)
  nze = ma*(((nza+ma-1)/ma)*((nzb+mb-1)/mb) )
  nzeb = (((nza+na-1)/na)*((nzb+nb-1)/nb))*nb
  ! Estimate number of nonzeros on output.
  ! Turns out this is often a large  overestimate.
  call c%allocate(ma,nb,min(nzc,nze,nzeb))


  call csc_spspmm(a,b,c,info)

  call c%set_asb()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains
  
  subroutine csc_spspmm(a,b,c,info)
    implicit none 
    type(psb_c_csc_sparse_mat), intent(in)  :: a,b
    type(psb_c_csc_sparse_mat), intent(inout) :: c
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_)              :: ma,na,mb,nb
    integer(psb_ipk_), allocatable :: icol(:), idxs(:), iaux(:)
    complex(psb_spk_), allocatable    :: col(:)
    type(psb_int_heap)             :: heap
    integer(psb_ipk_)              :: i,j,k,irw,icl,icf, iret, &
         & nzc,nnzre, isz, ipb, irwsz, nrc, nze
    complex(psb_spk_)                 :: cfb


    info = psb_success_
    ma = a%get_nrows()
    na = a%get_ncols()
    mb = b%get_nrows()
    nb = b%get_ncols()

    nze = min(size(c%val),size(c%ia))
    isz = max(ma,na,mb,nb)
    call psb_realloc(isz,col,info)
    if (info == 0) call psb_realloc(isz,idxs,info)
    if (info == 0) call psb_realloc(isz,icol,info)
    if (info /= 0) return 
    col  = dzero
    icol = 0
    nzc  = 1    
    do j = 1,nb
      c%icp(j) = nzc
      nrc = 0 
      do k = b%icp(j), b%icp(j+1)-1
        icl = b%ia(k)
        cfb = b%val(k)
        irwsz = a%icp(icl+1)-a%icp(icl)
        do i = a%icp(icl),a%icp(icl+1)-1
          irw = a%ia(i)
          if (icol(irw)<j) then 
            nrc = nrc + 1 
            idxs(nrc) = irw
            icol(irw) = j
          end if
          col(irw)  = col(irw)  + cfb*a%val(i)
        end do
      end do
      if (nrc > 0 ) then 
        if ((nzc+nrc)>nze) then 
          nze = max(nb*((nzc+j-1)/j),nzc+2*nrc)
          call psb_realloc(nze,c%val,info)
          if (info == 0) call psb_realloc(nze,c%ia,info)
          if (info /= 0) return
        end if
        call psb_msort(idxs(1:nrc))
        do i=1, nrc
          irw        = idxs(i) 
          c%ia(nzc)  = irw 
          c%val(nzc) = col(irw) 
          col(irw)   = dzero
          nzc        = nzc + 1 
        end do
      end if
    end do

    c%icp(nb+1) = nzc

  end subroutine csc_spspmm

end subroutine psb_ccscspspmm
