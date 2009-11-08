
!=====================================
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
!=====================================

subroutine c_csc_csmv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csmv_impl
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  complex(psb_spk_) :: acc
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='c_csc_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0 

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif



  tra  = (psb_toupper(trans_)=='T')
  ctra = (psb_toupper(trans_)=='C')

  if (tra.or.ctra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
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

  else if (ctra) then 

    if (beta == czero) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = alpha*acc
        end do

      end if


    else if (beta == cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = y(i) + alpha*acc
        end do

      end if

    else if (beta == -cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = -y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = -y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = -y(i) + alpha*acc
        end do

      end if

    else 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) - acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j))          
          enddo
          y(i) = beta*y(i) + alpha*acc
        end do

      end if

    end if

  else if (.not.(tra.or.ctra)) then 

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

  if (a%is_triangle().and.a%is_unit()) then 
    do i=1, min(m,n)
      y(i) = y(i) + alpha*x(i)
    end do
  end if

  call psb_erractionrestore(err_act)
  return
9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine c_csc_csmv_impl

subroutine c_csc_csmm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csmm_impl
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_), allocatable  :: acc(:)
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='c_csc_csmm'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  tra  = (psb_toupper(trans_)=='T')
  ctra = (psb_toupper(trans_)=='C')

  if (tra.or.ctra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if

  nc = min(size(x,2) , size(y,2) )

  allocate(acc(nc), stat=info)
  if(info /= 0) then
    info=4010
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

  else if (ctra) then 

    if (beta == czero) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = alpha*acc
        end do

      end if


    else if (beta == cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = y(i,:) + alpha*acc
        end do

      end if

    else if (beta == -cone) then 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = -y(i,:) + alpha*acc
        end do

      end if

    else 

      if (alpha == cone) then 
        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + acc
        end do

      else if (alpha == -cone) then 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) - acc
        end do

      else 

        do i=1,m 
          acc  = czero
          do j=a%icp(i), a%icp(i+1)-1
            acc  = acc + conjg(a%val(j)) * x(a%ia(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + alpha*acc
        end do

      end if

    end if

  else if (.not.(tra.or.ctra)) then 

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

  if (a%is_triangle().and.a%is_unit()) then 
    do i=1, min(m,n)
      y(i,:) = y(i,:) + alpha*x(i,:)
    end do
  end if

  call psb_erractionrestore(err_act)
  return
9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine c_csc_csmm_impl


subroutine c_csc_cssv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_cssv_impl
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: tmp(:)
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='c_csc_cssv'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)
  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  tra  = (psb_toupper(trans_)=='T')
  ctra = (psb_toupper(trans_)=='C')

  m = a%get_nrows()

  if (.not. (a%is_triangle())) then 
    info = 1121
    call psb_errpush(info,name)
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
    call inner_cscsv(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),&
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
    if (info /= 0) then 
      return
    end if
    tmp(1:m) = x(1:m)
    call inner_cscsv(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & a%icp,a%ia,a%val,tmp,y) 
    do  i = 1, m
      y(i) = alpha*tmp(i) + beta*y(i)
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains 

  subroutine inner_cscsv(tra,ctra,lower,unit,n,icp,ia,val,x,y) 
    implicit none 
    logical, intent(in)                 :: tra,ctra,lower,unit  
    integer, intent(in)                 :: icp(*), ia(*),n
    complex(psb_spk_), intent(in)          :: val(*)
    complex(psb_spk_), intent(in)          :: x(*)
    complex(psb_spk_), intent(out)         :: y(*)

    integer :: i,j,k,m, ir, jc
    complex(psb_spk_) :: acc

    if (tra) then 

      if (lower) then 
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
      else if (.not.lower) then 

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

      end if

    else if (ctra) then 

      if (lower) then 
        if (unit) then 
          do i=1, n
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + conjg(val(j))*y(ia(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = czero 
            do j=icp(i), icp(i+1)-2
              acc = acc + conjg(val(j))*y(ia(j))
            end do
            y(i) = (x(i) - acc)/conjg(val(icp(i+1)-1))
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=n, 1, -1 
            acc = czero 
            do j=icp(i), icp(i+1)-1
              acc = acc + conjg(val(j))*y(ia(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=n, 1, -1 
            acc = czero 
            do j=icp(i)+1, icp(i+1)-1
              acc = acc + conjg(val(j))*y(ia(j))
            end do
            y(i) = (x(i) - acc)/conjg(val(icp(i)))
          end do
        end if

      end if

    else if (.not.(tra.or.ctra)) then 

      do i=1, n
        y(i) = x(i)
      end do

      if (lower) then 
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
      else if (.not.lower) then 

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

      end if
    end if
  end subroutine inner_cscsv

end subroutine c_csc_cssv_impl



subroutine c_csc_cssm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_cssm_impl
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  complex(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  complex(psb_spk_) :: acc
  complex(psb_spk_), allocatable :: tmp(:,:)
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='c_base_csmm'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  tra  = (psb_toupper(trans_)=='T')
  ctra = (psb_toupper(trans_)=='C')

  m   = a%get_nrows()
  nc  = min(size(x,2) , size(y,2)) 

  if (.not. (a%is_triangle())) then 
    info = 1121
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
    call inner_cscsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%icp,a%ia,a%val,x,size(x,1),y,size(y,1),info) 
    do  i = 1, m
      y(i,1:nc) = alpha*y(i,1:nc)
    end do
  else 
    allocate(tmp(m,nc), stat=info) 
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='allocate')
      goto 9999
    end if

    tmp(1:m,:) = x(1:m,1:nc)
    call inner_cscsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%icp,a%ia,a%val,tmp,size(tmp,1),y,size(y,1),info) 
    do  i = 1, m
      y(i,1:nc) = alpha*tmp(i,1:nc) + beta*y(i,1:nc)
    end do
  end if

  if(info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='inner_cscsm')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


contains 

  subroutine inner_cscsm(tra,ctra,lower,unit,nr,nc,&
       & icp,ia,val,x,ldx,y,ldy,info) 
    implicit none 
    logical, intent(in)                 :: tra,ctra,lower,unit
    integer, intent(in)                 :: nr,nc,ldx,ldy,icp(*),ia(*)
    complex(psb_spk_), intent(in)       :: val(*), x(ldx,*)
    complex(psb_spk_), intent(out)      :: y(ldy,*)
    integer, intent(out)                :: info
    integer :: i,j,k,m, ir, jc
    complex(psb_spk_), allocatable  :: acc(:)

    info = 0
    allocate(acc(nc), stat=info)
    if(info /= 0) then
      info=4010
      return
    end if


    if (tra) then 

      if (lower) then 
        if (unit) then 
          do i=1, nr
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-1
              acc = acc + a%val(j)*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=1, nr
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-2
              acc = acc + a%val(j)*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/a%val(a%icp(i+1)-1)
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-1
              acc = acc + a%val(j)*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=a%icp(i)+1, a%icp(i+1)-1
              acc = acc + a%val(j)*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/a%val(a%icp(i))
          end do
        end if

      end if

    else if (ctra) then 

      if (lower) then 
        if (unit) then 
          do i=1, nr
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-1
              acc = acc + conjg(a%val(j))*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=1, nr
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-2
              acc = acc + conjg(a%val(j))*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/conjg(a%val(a%icp(i+1)-1))
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=a%icp(i), a%icp(i+1)-1
              acc = acc + conjg(a%val(j))*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1 
            acc = czero 
            do j=a%icp(i)+1, a%icp(i+1)-1
              acc = acc + conjg(a%val(j))*y(a%ia(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/conjg(a%val(a%icp(i)))
          end do
        end if

      end if

    else if (.not.(tra.or.ctra)) then 

      do i=1, nr
        y(i,1:nc) = x(i,1:nc)
      end do

      if (lower) then 
        if (unit) then 
          do i=nr, 1, -1
            acc = y(i,1:nc) 
            do j=a%icp(i), a%icp(i+1)-1
              jc    = a%ia(j)
              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1
            y(i,1:nc) = y(i,1:nc)/a%val(a%icp(i+1)-1)
            acc  = y(i,1:nc) 
            do j=a%icp(i), a%icp(i+1)-2
              jc    = a%ia(j)
              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
            end do
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=1, nr
            acc  = y(i,1:nc) 
            do j=a%icp(i), a%icp(i+1)-1
              jc    = a%ia(j)
              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, nr
            y(i,1:nc) = y(i,1:nc)/a%val(a%icp(i))
            acc    = y(i,1:nc) 
            do j=a%icp(i)+1, a%icp(i+1)-1
              jc      = a%ia(j)
              y(jc,1:nc) = y(jc,1:nc) - a%val(j)*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_cscsm

end subroutine c_csc_cssm_impl

function c_csc_csnmi_impl(a) result(res)
  use psb_error_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csnmi_impl
  implicit none 
  class(psb_c_csc_sparse_mat), intent(in) :: a
  real(psb_spk_)         :: res

  integer   :: i,j,k,m,n, nr, ir, jc, nc, info
  real(psb_spk_), allocatable  :: acc(:) 
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='c_csnmi'
  logical, parameter :: debug=.false.


  res = szero 
  nr = a%get_nrows()
  nc = a%get_ncols()
  allocate(acc(nr),stat=info)
  if (info /= 0) then 
    return
  end if
  acc(:) = szero
  do i=1, nc
    do j=a%icp(i),a%icp(i+1)-1  
      acc(a%ia(j)) = acc(a%ia(j)) + abs(a%val(j))
    end do
  end do
  do i=1, nr
    res = max(res,acc(i))
  end do
  deallocate(acc)

end function c_csc_csnmi_impl

!===================================== 
!
!
!
! Data management
!
!
!
!
!
!=====================================   

subroutine c_csc_csgetptn_impl(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csgetptn_impl
  implicit none

  class(psb_c_csc_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0

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
    info = 583
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

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
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
    integer                              :: imin,imax,jmin,jmax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    integer, intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer                              :: info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,m,isz
    integer  :: debug_level, debug_unit
    character(len=20) :: name='coo_getrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%get_nzeros()
    irw = imin
    lrw = min(imax,a%get_nrows())
    icl = jmin
    lcl = min(jmax,a%get_ncols())
    if (irw<0) then 
      info = 2
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif


    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
         & ((nza*(lcl+1-icl))/a%get_ncols()) )
    nz = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)

    if (info /= 0) return
    isz = min(size(ia),size(ja))
    if (present(iren)) then 
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
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
              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
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
  
end subroutine c_csc_csgetptn_impl




subroutine c_csc_csgetrow_impl(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csgetrow_impl
  implicit none

  class(psb_c_csc_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0

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
    info = 583
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

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
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
    integer                              :: imin,imax,jmin,jmax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
    integer, intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer                              :: info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw, icl, lcl,isz,m
    integer  :: debug_level, debug_unit
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
      info = 2
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    nzt = min((a%icp(lcl+1)-a%icp(icl)),&
         & ((nza*(lrw+1-irw))/m) )
    nz = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info==0) call psb_ensure_size(nzin_+nzt,ja,info)
    if (info==0) call psb_ensure_size(nzin_+nzt,val,info)

    if (info /= 0) return
    isz = min(size(ia),size(ja),size(val))
    if (present(iren)) then 
      do i=icl, lcl
        do j=a%icp(i), a%icp(i+1) - 1
          if ((imin <= a%ia(j)).and.(a%ia(j)<=imax)) then 
            nzin_ = nzin_ + 1
            if (nzin_>isz) then 
              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
              call psb_ensure_size(int(1.25*nzin_)+1,val,info)
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
              call psb_ensure_size(int(1.25*nzin_)+1,ia,info)
              call psb_ensure_size(int(1.25*nzin_)+1,ja,info)
              call psb_ensure_size(int(1.25*nzin_)+1,val,info)
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

end subroutine c_csc_csgetrow_impl



subroutine c_csc_csput_impl(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_csc_csput_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='c_csc_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)

  info = 0
  nza  = a%get_nzeros()

  if (a%is_bld()) then 
    ! Build phase should only ever be in COO
    info = 1121

  else  if (a%is_upd()) then 
    call  c_csc_srch_upd(nz,ia,ja,val,a,&
         & imin,imax,jmin,jmax,info,gtl)
    
    if (info /= 0) then  

      info = 1121
    end if

  else 
    ! State is wrong.
    info = 1121
  end if
  if (info /= 0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


contains

  subroutine c_csc_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_sort_mod
    implicit none 

    class(psb_c_csc_sparse_mat), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax
    integer, intent(in) :: ia(:),ja(:)
    complex(psb_spk_), intent(in) :: val(:)
    integer, intent(out) :: info
    integer, intent(in), optional  :: gtl(:)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nr,nc,nnz,dupl,ng, nar, nac
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='c_csc_srch_upd'

    info = 0
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
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ir,' in: ',i1,i2,&
                     & ' : ',a%ia(i1:i2-1)
                info = i
                return
              end if

            else

              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
            end if
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
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ir,' in: ',i1,i2,&
                     & ' : ',a%ia(i1:i2-1)
                info = i
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
            end if

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
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Was searching ',ir,' in: ',i1,i2,&
                   & ' : ',a%ia(i1:i2-1)
              info = i
              return
            end if

          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding col that does not belong to us.'
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
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Was searching ',ir,' in: ',i1,i2,&
                   & ' : ',a%ia(i1:i2-1)
              info = i
              return
            end if
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding col that does not belong to us.'
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine c_csc_srch_upd

end subroutine c_csc_csput_impl



subroutine c_cp_csc_from_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_cp_csc_from_coo_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psb_c_coo_sparse_mat)   :: tmp
  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0
  ! This is to have fix_coo called behind the scenes
  call tmp%cp_from_coo(b,info)
  if (info ==0) call a%mv_from_coo(tmp,info)

end subroutine c_cp_csc_from_coo_impl



subroutine c_cp_csc_to_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_cp_csc_to_coo_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(in)  :: a
  class(psb_c_coo_sparse_mat), intent(out) :: b
  integer, intent(out)                      :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, idl,err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  call b%psb_c_base_sparse_mat%cp_from(a%psb_c_base_sparse_mat)

  do i=1, nc
    do j=a%icp(i),a%icp(i+1)-1
      b%ia(j)  = a%ia(j)
      b%ja(j)  = i
      b%val(j) = a%val(j)
    end do
  end do
  
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)


end subroutine c_cp_csc_to_coo_impl


subroutine c_mv_csc_to_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_mv_csc_to_coo_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(out)   :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, idl,err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%psb_c_base_sparse_mat%mv_from(a%psb_c_base_sparse_mat)
  call b%set_nzeros(a%get_nzeros())
  call move_alloc(a%ia,b%ia)
  call move_alloc(a%val,b%val)
  call psb_realloc(nza,b%ja,info)
  if (info /= 0) return
  do i=1, nc
    do j=a%icp(i),a%icp(i+1)-1
      b%ja(j)  = i
    end do
  end do
  call a%free()
  call b%fix(info)

end subroutine c_mv_csc_to_coo_impl



subroutine c_mv_csc_from_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_mv_csc_from_coo_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc, icl
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  call b%fix(info, idir=1)
  if (info /= 0) return

  nr  = b%get_nrows()
  nc  = b%get_ncols()
  nza = b%get_nzeros()
  
  call a%psb_c_base_sparse_mat%mv_from(b%psb_c_base_sparse_mat)

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
               & 'Strange situation: i>nr ',i,nc,j,nza,icl,idl
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


end subroutine c_mv_csc_from_coo_impl


subroutine c_mv_csc_to_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_mv_csc_to_fmt_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(out)  :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
    ! Need to fix trivial copies! 
  type is (psb_c_csc_sparse_mat) 
    call b%psb_c_base_sparse_mat%mv_from(a%psb_c_base_sparse_mat)
    call move_alloc(a%icp, b%icp)
    call move_alloc(a%ia,  b%ia)
    call move_alloc(a%val, b%val)
    call a%free()
    
  class default
    call tmp%mv_from_fmt(a,info)
    if (info == 0) call b%mv_from_coo(tmp,info)
  end select

end subroutine c_mv_csc_to_fmt_impl
!!$

subroutine c_cp_csc_to_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_cp_csc_to_fmt_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(in)   :: a
  class(psb_c_base_sparse_mat), intent(out) :: b
  integer, intent(out)                       :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0


  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    call b%psb_c_base_sparse_mat%cp_from(a%psb_c_base_sparse_mat)
    b%icp = a%icp
    b%ia  = a%ia
    b%val = a%val

  class default
    call tmp%cp_from_fmt(a,info)
    if (info == 0) call b%mv_from_coo(tmp,info)
  end select

end subroutine c_cp_csc_to_fmt_impl


subroutine c_mv_csc_from_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_mv_csc_from_fmt_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout)  :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    call a%psb_c_base_sparse_mat%mv_from(b%psb_c_base_sparse_mat)
    call move_alloc(b%icp, a%icp)
    call move_alloc(b%ia,  a%ia)
    call move_alloc(b%val, a%val)
    call b%free()

  class default
    call tmp%mv_from_fmt(b,info)
    if (info == 0) call a%mv_from_coo(tmp,info)
  end select

end subroutine c_mv_csc_from_fmt_impl



subroutine c_cp_csc_from_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_c_base_mat_mod
  use psb_c_csc_mat_mod, psb_protect_name => c_cp_csc_from_fmt_impl
  implicit none 

  class(psb_c_csc_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_c_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_c_csc_sparse_mat) 
    call a%psb_c_base_sparse_mat%cp_from(b%psb_c_base_sparse_mat)
    a%icp = b%icp
    a%ia  = b%ia
    a%val = b%val

  class default
    call tmp%cp_from_fmt(b,info)
    if (info == 0) call a%mv_from_coo(tmp,info)
  end select
end subroutine c_cp_csc_from_fmt_impl

