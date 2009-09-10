
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine d_cxx_csmv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_csmv_impl
  implicit none 
  class(psbn_d_cxx_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_cxx_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = 0 

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    write(0,*) 'Error: csmv called on an unassembled mat'
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  tra = ((trans_=='T').or.(trans_=='t'))

  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if


  if (alpha == dzero) then
    if (beta == dzero) then
      do i = 1, m
        y(i) = dzero
      enddo
    else
      do  i = 1, m
        y(i) = beta*y(i)
      end do
    endif
    return
  end if

  if (.not.tra) then 

    if (beta == dzero) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = alpha*acc
        end do

      end if


    else if (beta == done) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = y(i) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = y(i) + alpha*acc
        end do

      end if

    else if (beta == -done) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = -y(i) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = -y(i) -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = -y(i) + alpha*acc
        end do

      end if

    else 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = beta*y(i) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = beta*y(i) - acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j))          
          enddo
          y(i) = beta*y(i) + alpha*acc
        end do

      end if

    end if

  else if (tra) then 

    if (beta == dzero) then 
      do i=1, m
        y(i) = dzero
      end do
    else if (beta == done) then 
      ! Do nothing
    else if (beta == -done) then 
      do i=1, m
        y(i) = -y(i) 
      end do
    else
      do i=1, m
        y(i) = beta*y(i) 
      end do
    end if

    if (alpha.eq.done) then

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
          y(ir) = y(ir) +  a%val(j)*x(i)
        end do
      enddo

    else if (alpha.eq.-done) then

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
          y(ir) = y(ir) -  a%val(j)*x(i)
        end do
      enddo

    else                    

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
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
9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine d_cxx_csmv_impl

subroutine d_cxx_csmm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_csmm_impl
  implicit none 
  class(psbn_d_cxx_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_cxx_csmm'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  tra = ((trans_=='T').or.(trans_=='t'))
  if (.not.a%is_asb()) then 
    info = 1121
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

  nc = min(size(x,2) , size(y,2) )

  allocate(acc(nc), stat=info)
  if(info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='allocate')
    goto 9999
  end if

  if (alpha == dzero) then
    if (beta == dzero) then
      do i = 1, m
        y(i,:) = dzero
      enddo
    else
      do  i = 1, m
        y(i,:) = beta*y(i,:)
      end do
    endif
    return
  end if

  if (.not.tra) then 

    if (beta == dzero) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = alpha*acc
        end do

      end if


    else if (beta == done) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = y(i,:) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = y(i,:) + alpha*acc
        end do

      end if

    else if (beta == -done) then 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = -y(i,:) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = -y(i,:) -acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = -y(i,:) + alpha*acc
        end do

      end if

    else 

      if (alpha == done) then 
        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + acc
        end do

      else if (alpha == -done) then 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = beta*y(i,:) - acc
        end do

      else 

        do i=1,m 
          acc  = dzero
          do j=a%irp(i), a%irp(i+1)-1
            acc  = acc + a%val(j) * x(a%ja(j),:)          
          enddo
          y(i,:) = beta*y(i,:) + alpha*acc
        end do

      end if

    end if

  else if (tra) then 

    if (beta == dzero) then 
      do i=1, m
        y(i,:) = dzero
      end do
    else if (beta == done) then 
      ! Do nothing
    else if (beta == -done) then 
      do i=1, m
        y(i,:) = -y(i,:) 
      end do
    else
      do i=1, m
        y(i,:) = beta*y(i,:) 
      end do
    end if

    if (alpha.eq.done) then

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
          y(ir,:) = y(ir,:) +  a%val(j)*x(i,:)
        end do
      enddo

    else if (alpha.eq.-done) then

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
          y(ir,:) = y(ir,:) -  a%val(j)*x(i,:)
        end do
      enddo

    else                    

      do i=1,n
        do j=a%irp(i), a%irp(i+1)-1
          ir = a%ja(j)
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
9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine d_cxx_csmm_impl


subroutine d_cxx_cssv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_cssv_impl
  implicit none 
  class(psbn_d_cxx_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_cxx_cssv'
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

  tra = ((trans_=='T').or.(trans_=='t'))
  m = a%get_nrows()

  if (.not. (a%is_triangle())) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if


  if (alpha == dzero) then
    if (beta == dzero) then
      do i = 1, m
        y(i) = dzero
      enddo
    else
      do  i = 1, m
        y(i) = beta*y(i)
      end do
    endif
    return
  end if

  if (beta == dzero) then 
    call inner_cxxsv(tra,a,x,y)
    do  i = 1, m
      y(i) = alpha*y(i)
    end do
  else 
    allocate(tmp(m), stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation error in CXXSV '
      return
    end if
    tmp(1:m) = x(1:m)
    call inner_cxxsv(tra,a,tmp,y)
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

  subroutine inner_cxxsv(tra,a,x,y) 
    implicit none 
    logical, intent(in)                 :: tra  
    class(psbn_d_cxx_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: x(:)
    real(psb_dpk_), intent(out)         :: y(:)

    integer :: i,j,k,m, ir, jc
    real(psb_dpk_) :: acc

    if (.not.tra) then 

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          do i=1, a%get_nrows()
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.a%is_unit()) then 
          do i=1, a%get_nrows()
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-2
              acc = acc + a%val(j)*x(a%ja(j))
            end do
            y(i) = (x(i) - acc)/a%val(a%irp(i+1)-1)
          end do
        end if
      else if (a%is_upper()) then 

        if (a%is_unit()) then 
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.a%is_unit()) then 
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do j=a%irp(i)+1, a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j))
            end do
            y(i) = (x(i) - acc)/a%val(a%irp(i))
          end do
        end if

      end if

    else if (tra) then 

      do i=1, a%get_nrows()
        y(i) = x(i)
      end do

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          do i=a%get_nrows(), 1, -1
            acc = y(i) 
            do j=a%irp(i), a%irp(i+1)-1
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
            end do
          end do
        else if (.not.a%is_unit()) then 
          do i=a%get_nrows(), 1, -1
            y(i) = y(i)/a%val(a%irp(i+1)-1)
            acc  = y(i) 
            do j=a%irp(i), a%irp(i+1)-2
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
            end do
          end do
        end if
      else if (a%is_upper()) then 

        if (a%is_unit()) then 
          do i=1, a%get_nrows()
            acc  = y(i) 
            do j=a%irp(i), a%irp(i+1)-1
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
            end do
          end do
        else if (.not.a%is_unit()) then 
          do i=1, a%get_nrows()
            y(i) = y(i)/a%val(a%irp(i))
            acc  = y(i) 
            do j=a%irp(i)+1, a%irp(i+1)-1
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_cxxsv

end subroutine d_cxx_cssv_impl



subroutine d_cxx_cssm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_cssm_impl
  implicit none 
  class(psbn_d_cxx_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:,:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_base_csmm'
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


  tra = ((trans_=='T').or.(trans_=='t'))
  m   = a%get_nrows()
  nc  = min(size(x,2) , size(y,2)) 

  if (.not. (a%is_triangle())) then 
    write(0,*) 'Called SM on a non-triangular mat!'
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if


  if (alpha == dzero) then
    if (beta == dzero) then
      do i = 1, m
        y(i,:) = dzero
      enddo
    else
      do  i = 1, m
        y(i,:) = beta*y(i,:)
      end do
    endif
    return
  end if

  if (beta == dzero) then 
    call inner_cxxsm(tra,a,x,y,info)
    do  i = 1, m
      y(i,:) = alpha*y(i,:)
    end do
  else 
    allocate(tmp(m,nc), stat=info) 
    if(info /= 0) then
      info=4010
      call psb_errpush(info,name,a_err='allocate')
      goto 9999
    end if

    tmp(1:m,:) = x(1:m,:)
    call inner_cxxsm(tra,a,tmp,y,info)
    do  i = 1, m
      y(i,:) = alpha*tmp(i,:) + beta*y(i,:)
    end do
  end if

  if(info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='inner_cxxsm')
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

  subroutine inner_cxxsm(tra,a,x,y,info) 
    implicit none 
    logical, intent(in)                 :: tra  
    class(psbn_d_cxx_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: x(:,:)
    real(psb_dpk_), intent(out)         :: y(:,:)
    integer, intent(out)                :: info
    integer :: i,j,k,m, ir, jc
    real(psb_dpk_), allocatable  :: acc(:)

    info = 0
    allocate(acc(size(x,2)), stat=info)
    if(info /= 0) then
      info=4010
      return
    end if


    if (.not.tra) then 

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          do i=1, a%get_nrows()
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do
        else if (.not.a%is_unit()) then 
          do i=1, a%get_nrows()
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-2
              acc = acc + a%val(j)*x(a%ja(j),:)
            end do
            y(i,:) = (x(i,:) - acc)/a%val(a%irp(i+1)-1)
          end do
        end if
      else if (a%is_upper()) then 

        if (a%is_unit()) then 
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do j=a%irp(i), a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do
        else if (.not.a%is_unit()) then 
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do j=a%irp(i)+1, a%irp(i+1)-1
              acc = acc + a%val(j)*x(a%ja(j),:)
            end do
            y(i,:) = (x(i,:) - acc)/a%val(a%irp(i))
          end do
        end if

      end if

    else if (tra) then 

      do i=1, a%get_nrows()
        y(i,:) = x(i,:)
      end do

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          do i=a%get_nrows(), 1, -1
            acc = y(i,:) 
            do j=a%irp(i), a%irp(i+1)-1
              jc    = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
            end do
          end do
        else if (.not.a%is_unit()) then 
          do i=a%get_nrows(), 1, -1
            y(i,:) = y(i,:)/a%val(a%irp(i+1)-1)
            acc  = y(i,:) 
            do j=a%irp(i), a%irp(i+1)-2
              jc    = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
            end do
          end do
        end if
      else if (a%is_upper()) then 

        if (a%is_unit()) then 
          do i=1, a%get_nrows()
            acc  = y(i,:) 
            do j=a%irp(i), a%irp(i+1)-1
              jc    = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
            end do
          end do
        else if (.not.a%is_unit()) then 
          do i=1, a%get_nrows()
            y(i,:) = y(i,:)/a%val(a%irp(i))
            acc    = y(i,:) 
            do j=a%irp(i)+1, a%irp(i+1)-1
              jc      = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_cxxsm

end subroutine d_cxx_cssm_impl

function d_cxx_csnmi_impl(a) result(res)
  use psb_error_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_csnmi_impl
  implicit none 
  class(psbn_d_cxx_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nr, ir, jc, nc
  real(psb_dpk_) :: acc
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_csnmi'
  logical, parameter :: debug=.false.


  res = dzero 
 
  do i = 1, a%get_nrows()
    acc = dzero
    do j=a%irp(i),a%irp(i+1)-1  
      acc = acc + abs(a%val(j))
    end do
    res = max(res,acc)
  end do

end function d_cxx_csnmi_impl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !
  ! Data management
  !
  !
  !
  !
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  



subroutine d_cxx_csput_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cxx_csput_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_cxx_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)

  info = 0
  nza  = a%get_nzeros()

  if (a%is_bld()) then 
    ! Build phase should only ever be in COO
    info = 1121

  else  if (a%is_upd()) then 
    call  d_cxx_srch_upd(nz,ia,ja,val,a,&
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

  subroutine d_cxx_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    class(psbn_d_cxx_sparse_mat), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax
    integer, intent(in) :: ia(:),ja(:)
    real(psb_dpk_), intent(in) :: val(:)
    integer, intent(out) :: info
    integer, intent(in), optional  :: gtl(:)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,nnz,dupl,ng
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='d_cxx_srch_upd'

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

    if (present(gtl)) then 
      ng = size(gtl)

      select case(dupl)
      case(psbn_dupl_ovwrt_,psbn_dupl_err_)
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
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%irp(ir)
              i2 = a%irp(ir+1)
              nc=i2-i1

              ip = psb_ibsrch(ic,nc,a%ja(i1:i2-1))    
              if (ip>0) then 
                a%val(i1+ip-1) = val(i)
              else
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ic,' in: ',i1,i2,&
                     & ' : ',a%ja(i1:i2-1)
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

      case(psbn_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            ic = gtl(ic)
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%irp(ir)
              i2 = a%irp(ir+1)
              nc = i2-i1
              ip = psb_ibsrch(ic,nc,a%ja(i1:i2-1))
              if (ip>0) then 
                a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
              else
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ic,' in: ',i1,i2,&
                     & ' : ',a%ja(i1:i2-1)
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
      case(psbn_dupl_ovwrt_,psbn_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 

          if ((ir > 0).and.(ir <= a%m)) then 

            i1 = a%irp(ir)
            i2 = a%irp(ir+1)
            nc=i2-i1

            ip = psb_ibsrch(ic,nc,a%ja(i1:i2-1))    
            if (ip>0) then 
              a%val(i1+ip-1) = val(i)
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Was searching ',ic,' in: ',i1,i2,&
                   & ' : ',a%ja(i1:i2-1)
              info = i
              return
            end if

          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if

        end do

      case(psbn_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 
            i1 = a%irp(ir)
            i2 = a%irp(ir+1)
            nc = i2-i1
            ip = psb_ibsrch(ic,nc,a%ja(i1:i2-1))
            if (ip>0) then 
              a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
            else
              info = i
              return
            end if
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine d_cxx_srch_upd

end subroutine d_cxx_csput_impl



subroutine d_cp_cxx_from_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cp_cxx_from_coo_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psbn_d_coo_sparse_mat)   :: tmp
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

end subroutine d_cp_cxx_from_coo_impl



subroutine d_cp_cxx_to_coo_impl(a,b,info) 
  use psb_const_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cp_cxx_to_coo_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(in)  :: a
  class(psbn_d_coo_sparse_mat), intent(out) :: b
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

  call b%allocate(nr,nr,nza)

  do i=1, nr
    do j=a%irp(i),a%irp(i+1)-1
      b%ia(j)  = i
      b%ja(j)  = a%ja(j)
      b%val(j) = a%val(j)
    end do
  end do
  
  call b%set_nzeros(a%get_nzeros())
  call b%set_nrows(a%get_nrows())
  call b%set_ncols(a%get_ncols())
  call b%set_dupl(a%get_dupl())
  call b%set_state(a%get_state())
  call b%set_triangle(a%is_triangle())
  call b%set_upper(a%is_upper())
  call b%fix(info)


end subroutine d_cp_cxx_to_coo_impl


subroutine d_mv_cxx_to_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_mv_cxx_to_coo_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(out)   :: b
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

  call b%set_nzeros(a%get_nzeros())
  call b%set_nrows(a%get_nrows())
  call b%set_ncols(a%get_ncols())
  call b%set_dupl(a%get_dupl())
  call b%set_state(a%get_state())
  call b%set_triangle(a%is_triangle())
  call b%set_upper(a%is_upper())

  call move_alloc(a%ja,b%ja)
  call move_alloc(a%val,b%val)
  call psb_realloc(nza,b%ia,info)
  if (info /= 0) return
  do i=1, nr
    do j=a%irp(i),a%irp(i+1)-1
      b%ia(j)  = i
    end do
  end do
  call a%free()
  call b%fix(info)


end subroutine d_mv_cxx_to_coo_impl



subroutine d_mv_cxx_from_coo_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_mv_cxx_from_coo_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  call b%fix(info)
  if (info /= 0) return

  nr  = b%get_nrows()
  nc  = b%get_ncols()
  nza = b%get_nzeros()

  call a%set_nrows(b%get_nrows())
  call a%set_ncols(b%get_ncols())
  call a%set_dupl(b%get_dupl())
  call a%set_state(b%get_state())
  call a%set_triangle(b%is_triangle())
  call a%set_upper(b%is_upper())
  ! Dirty trick: call move_alloc to have the new data allocated just once.
  call move_alloc(b%ia,itemp)
  call move_alloc(b%ja,a%ja)
  call move_alloc(b%val,a%val)
  call psb_realloc(max(nr+1,nc+1),a%irp,info)
  call b%free()

  if (nza <= 0) then 
    a%irp(:) = 1
  else
    a%irp(1) = 1
    if (nr < itemp(nza)) then 
      write(debug_unit,*) trim(name),': RWSHR=.false. : ',&
           &nr,itemp(nza),' Expect trouble!'
      info = 12
    end if

    j = 1 
    i = 1
    irw = itemp(j) 

    outer: do 
      inner: do 
        if (i >= irw) exit inner
        if (i>nr) then 
          write(debug_unit,*) trim(name),&
               & 'Strange situation: i>nr ',i,nr,j,nza,irw,idl
          exit outer
        end if
        a%irp(i+1) = a%irp(i) 
        i = i + 1
      end do inner
      j = j + 1
      if (j > nza) exit
      if (itemp(j) /= irw) then 
        a%irp(i+1) = j
        irw = itemp(j) 
        i = i + 1
      endif
      if (i>nr) exit
    enddo outer
    !
    ! Cleanup empty rows at the end
    !
    if (j /= (nza+1)) then 
      write(debug_unit,*) trim(name),': Problem from loop :',j,nza
      info = 13
    endif
    do 
      if (i>nr) exit
      a%irp(i+1) = j
      i = i + 1
    end do

  endif


end subroutine d_mv_cxx_from_coo_impl


subroutine d_mv_cxx_to_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_mv_cxx_to_fmt_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  class(psbn_d_base_sparse_mat), intent(out)  :: b
  integer, intent(out)                        :: info

  !locals
  type(psbn_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  class is (psbn_d_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
  class default
    call tmp%mv_from_fmt(a,info)
    if (info == 0) call b%mv_from_coo(tmp,info)
  end select

end subroutine d_mv_cxx_to_fmt_impl


subroutine d_cp_cxx_to_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cp_cxx_to_fmt_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(in)   :: a
  class(psbn_d_base_sparse_mat), intent(out) :: b
  integer, intent(out)                       :: info

  !locals
  type(psbn_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0


  select type (b)
  class is (psbn_d_coo_sparse_mat) 
    call a%cp_to_coo(b,info)
  class default
    call tmp%cp_from_fmt(a,info)
    if (info == 0) call b%mv_from_coo(tmp,info)
  end select

end subroutine d_cp_cxx_to_fmt_impl


subroutine d_mv_cxx_from_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_mv_cxx_from_fmt_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout)  :: a
  class(psbn_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info

  !locals
  type(psbn_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  class is (psbn_d_coo_sparse_mat) 
    call a%mv_from_coo(b,info)
  class default
    call tmp%mv_from_fmt(b,info)
    if (info == 0) call a%mv_from_coo(tmp,info)
  end select

end subroutine d_mv_cxx_from_fmt_impl



subroutine d_cp_cxx_from_fmt_impl(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod
  use psbn_d_cxx_mat_mod, psb_protect_name => d_cp_cxx_from_fmt_impl
  implicit none 

  class(psbn_d_cxx_sparse_mat), intent(inout) :: a
  class(psbn_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psbn_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, idl,err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = 0

  select type (b)
  class is (psbn_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)
  class default
    call tmp%cp_from_fmt(b,info)
    if (info == 0) call a%mv_from_coo(tmp,info)
  end select
end subroutine d_cp_cxx_from_fmt_impl
