
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

subroutine psb_d_csr_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csmv
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='d_csr_csmv'
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
    call psb_errpush(info,name,i_err=(/3,n,0,0,0/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
    goto 9999
  end if


  call psb_d_csr_csmv_inner(m,n,alpha,a%irp,a%ja,a%val,&
       & a%is_triangle(),a%is_unit(),&
       & x,beta,y,tra,ctra) 

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
  subroutine psb_d_csr_csmv_inner(m,n,alpha,irp,ja,val,is_triangle,is_unit,&
       & x,beta,y,tra,ctra) 
    integer, intent(in)             :: m,n,irp(*),ja(*)
    real(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(*)
    real(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)             :: is_triangle,is_unit,tra, ctra


    integer   :: i,j,k, ir, jc
    real(psb_dpk_) :: acc

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


    if ((.not.tra).and.(.not.ctra)) then 

      if (beta == dzero) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = -acc
          end do

        else 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = alpha*acc
          end do

        end if


      else if (beta == done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = y(i) -acc
          end do

        else 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = y(i) + alpha*acc
          end do

        end if

      else if (beta == -done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = -y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = -y(i) -acc
          end do

        else 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = -y(i) + alpha*acc
          end do

        end if

      else 

        if (alpha == done) then 
          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = beta*y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
            enddo
            y(i) = beta*y(i) - acc
          end do

        else 

          do i=1,m 
            acc  = dzero
            do j=irp(i), irp(i+1)-1
              acc  = acc + val(j) * x(ja(j))          
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

      if (alpha == done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) +  val(j)*x(i)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) -  val(j)*x(i)
          end do
        enddo

      else                    

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) + alpha*val(j)*x(i)
          end do
        enddo

      end if

    else if (ctra) then 

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

      if (alpha == done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) +  (val(j))*x(i)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) -  (val(j))*x(i)
          end do
        enddo

      else                    

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir) = y(ir) + alpha*(val(j))*x(i)
          end do
        enddo

      end if

    endif

    if (is_triangle.and.is_unit) then 
      do i=1, min(m,n)
        y(i) = y(i) + alpha*x(i)
      end do
    end if


  end subroutine psb_d_csr_csmv_inner


end subroutine psb_d_csr_csmv

subroutine psb_d_csr_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csmm
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='d_csr_csmm'
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
    call psb_errpush(info,name,i_err=(/3,n,0,0,0/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
    goto 9999
  end if

  nc = min(size(x,2) , size(y,2) )

  allocate(acc(nc), stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='allocate')
    goto 9999
  end if
  
  call  psb_d_csr_csmm_inner(m,n,nc,alpha,a%irp,a%ja,a%val, &
       & a%is_triangle(),a%is_unit(),x,size(x,1), &
       & beta,y,size(y,1),tra,ctra,acc) 


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
  subroutine psb_d_csr_csmm_inner(m,n,nc,alpha,irp,ja,val,&
       & is_triangle,is_unit,x,ldx,beta,y,ldy,tra,ctra,acc) 
    integer, intent(in)             :: m,n,ldx,ldy,nc,irp(*),ja(*)
    real(psb_dpk_), intent(in)      :: alpha, beta, x(ldx,*),val(*)
    real(psb_dpk_), intent(inout)   :: y(ldy,*)
    logical, intent(in)             :: is_triangle,is_unit,tra,ctra

    real(psb_dpk_), intent(inout)   :: acc(*)
    integer   :: i,j,k, ir, jc


    if (alpha == dzero) then
      if (beta == dzero) then
        do i = 1, m
          y(i,1:nc) = dzero
        enddo
      else
        do  i = 1, m
          y(i,1:nc) = beta*y(i,1:nc)
        end do
      endif
      return
    end if

    if ((.not.tra).and.(.not.ctra)) then 
      if (beta == dzero) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = acc(1:nc)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = -acc(1:nc)
          end do

        else 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = alpha*acc(1:nc)
          end do

        end if


      else if (beta == done) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = y(i,1:nc) + acc(1:nc)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = y(i,1:nc) -acc(1:nc)
          end do

        else 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = y(i,1:nc) + alpha*acc(1:nc)
          end do

        end if

      else if (beta == -done) then 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = -y(i,1:nc) + acc(1:nc)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = -y(i,1:nc) -acc(1:nc)
          end do

        else 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = -y(i,1:nc) + alpha*acc(1:nc)
          end do

        end if

      else 

        if (alpha == done) then 
          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = beta*y(i,1:nc) + acc(1:nc)
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = beta*y(i,1:nc) - acc(1:nc)
          end do

        else 

          do i=1,m 
            acc(1:nc)  = dzero
            do j=irp(i), irp(i+1)-1
              acc(1:nc)  = acc(1:nc) + val(j) * x(ja(j),1:nc)          
            enddo
            y(i,1:nc) = beta*y(i,1:nc) + alpha*acc(1:nc)
          end do

        end if

      end if

    else if (tra) then 

      if (beta == dzero) then 
        do i=1, m
          y(i,1:nc) = dzero
        end do
      else if (beta == done) then 
        ! Do nothing
      else if (beta == -done) then 
        do i=1, m
          y(i,1:nc) = -y(i,1:nc) 
        end do
      else
        do i=1, m
          y(i,1:nc) = beta*y(i,1:nc) 
        end do
      end if

      if (alpha == done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) +  val(j)*x(i,1:nc)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) -  val(j)*x(i,1:nc)
          end do
        enddo

      else                    

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) + alpha*val(j)*x(i,1:nc)
          end do
        enddo

      end if

    else if (ctra) then 

      if (beta == dzero) then 
        do i=1, m
          y(i,1:nc) = dzero
        end do
      else if (beta == done) then 
        ! Do nothing
      else if (beta == -done) then 
        do i=1, m
          y(i,1:nc) = -y(i,1:nc) 
        end do
      else
        do i=1, m
          y(i,1:nc) = beta*y(i,1:nc) 
        end do
      end if

      if (alpha == done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) +  (val(j))*x(i,1:nc)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) -  (val(j))*x(i,1:nc)
          end do
        enddo

      else                    

        do i=1,n
          do j=irp(i), irp(i+1)-1
            ir = ja(j)
            y(ir,1:nc) = y(ir,1:nc) + alpha*(val(j))*x(i,1:nc)
          end do
        enddo

      end if

    endif

    if (is_triangle.and.is_unit) then 
      do i=1, min(m,n)
        y(i,1:nc) = y(i,1:nc) + alpha*x(i,1:nc)
      end do
    end if

  end subroutine psb_d_csr_csmm_inner

end subroutine psb_d_csr_csmm


subroutine psb_d_csr_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_cssv
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:)
  logical   :: tra,ctra
  Integer :: err_act
  character(len=20)  :: name='d_csr_cssv'
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

  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')
  m = a%get_nrows()

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(x)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/3,m,0,0,0/))
    goto 9999
  end if

  if (size(y)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
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

    call inner_csrsv(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & a%irp,a%ja,a%val,x,y) 
    if (alpha == done) then 
      ! do nothing
    else if (alpha == -done) then 
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

    call inner_csrsv(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & a%irp,a%ja,a%val,x,tmp) 

    call psb_geaxpby(m,alpha,tmp,beta,y,info)

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

  subroutine inner_csrsv(tra,ctra,lower,unit,n,irp,ja,val,x,y) 
    implicit none 
    logical, intent(in)            :: tra,ctra,lower,unit  
    integer, intent(in)            :: irp(*), ja(*),n
    real(psb_dpk_), intent(in)  :: val(*)
    real(psb_dpk_), intent(in)  :: x(*)
    real(psb_dpk_), intent(out) :: y(*)

    integer :: i,j,k,m, ir, jc
    real(psb_dpk_) :: acc

    if ((.not.tra).and.(.not.ctra)) then 

      if (lower) then 
        if (unit) then 
          do i=1, n
            acc = dzero 
            do j=irp(i), irp(i+1)-1
              acc = acc + val(j)*y(ja(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = dzero 
            do j=irp(i), irp(i+1)-2
              acc = acc + val(j)*y(ja(j))
            end do
            y(i) = (x(i) - acc)/val(irp(i+1)-1)
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=n, 1, -1 
            acc = dzero 
            do j=irp(i), irp(i+1)-1
              acc = acc + val(j)*y(ja(j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=n, 1, -1 
            acc = dzero 
            do j=irp(i)+1, irp(i+1)-1
              acc = acc + val(j)*y(ja(j))
            end do
            y(i) = (x(i) - acc)/val(irp(i))
          end do
        end if

      end if

    else if (tra) then 

      do i=1, n
        y(i) = x(i)
      end do

      if (lower) then 
        if (unit) then 
          do i=n, 1, -1
            acc = y(i) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=n, 1, -1
            y(i) = y(i)/val(irp(i+1)-1)
            acc  = y(i) 
            do j=irp(i), irp(i+1)-2
              jc    = ja(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=1, n
            acc  = y(i) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, n
            y(i) = y(i)/val(irp(i))
            acc  = y(i) 
            do j=irp(i)+1, irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - val(j)*acc 
            end do
          end do
        end if

      end if

    else if (ctra) then 

      do i=1, n
        y(i) = x(i)
      end do

      if (lower) then 
        if (unit) then 
          do i=n, 1, -1
            acc = y(i) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - (val(j))*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=n, 1, -1
            y(i) = y(i)/val(irp(i+1)-1)
            acc  = y(i) 
            do j=irp(i), irp(i+1)-2
              jc    = ja(j)
              y(jc) = y(jc) - (val(j))*acc 
            end do
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=1, n
            acc  = y(i) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - (val(j))*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, n
            y(i) = y(i)/val(irp(i))
            acc  = y(i) 
            do j=irp(i)+1, irp(i+1)-1
              jc    = ja(j)
              y(jc) = y(jc) - (val(j))*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_csrsv

end subroutine psb_d_csr_cssv



subroutine psb_d_csr_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_string_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_cssm
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:,:)
  logical   :: tra, ctra
  Integer :: err_act
  character(len=20)  :: name='d_csr_cssm'
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


  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')

  m   = a%get_nrows()
  nc  = min(size(x,2) , size(y,2)) 

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
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
    call inner_csrsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%irp,a%ja,a%val,x,size(x,1),y,size(y,1),info) 
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

    call inner_csrsm(tra,ctra,a%is_lower(),a%is_unit(),a%get_nrows(),nc,&
         & a%irp,a%ja,a%val,x,size(x,1),tmp,size(tmp,1),info) 
    do  i = 1, m
      y(i,1:nc) = alpha*tmp(i,1:nc) + beta*y(i,1:nc)
    end do
  end if

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='inner_csrsm')
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

  subroutine inner_csrsm(tra,ctra,lower,unit,nr,nc,&
       & irp,ja,val,x,ldx,y,ldy,info) 
    implicit none 
    logical, intent(in)              :: tra,ctra,lower,unit
    integer, intent(in)              :: nr,nc,ldx,ldy,irp(*),ja(*)
    real(psb_dpk_), intent(in)    :: val(*), x(ldx,*)
    real(psb_dpk_), intent(out)   :: y(ldy,*)
    integer, intent(out)             :: info
    integer :: i,j,k,m, ir, jc
    real(psb_dpk_), allocatable  :: acc(:)

    info = psb_success_
    allocate(acc(nc), stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      return
    end if


    if ((.not.tra).and.(.not.ctra)) then 
      if (lower) then 
        if (unit) then 
          do i=1, nr
            acc = dzero 
            do j=irp(i), irp(i+1)-1
              acc = acc + val(j)*y(ja(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=1, nr
            acc = dzero 
            do j=irp(i), irp(i+1)-2
              acc = acc + val(j)*y(ja(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/val(irp(i+1)-1)
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=nr, 1, -1 
            acc = dzero 
            do j=irp(i), irp(i+1)-1
              acc = acc + val(j)*y(ja(j),1:nc)
            end do
            y(i,1:nc) = x(i,1:nc) - acc
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1 
            acc = dzero 
            do j=irp(i)+1, irp(i+1)-1
              acc = acc + val(j)*y(ja(j),1:nc)
            end do
            y(i,1:nc) = (x(i,1:nc) - acc)/val(irp(i))
          end do
        end if

      end if

    else if (tra) then 

      do i=1, nr
        y(i,1:nc) = x(i,1:nc)
      end do

      if (lower) then 
        if (unit) then 
          do i=nr, 1, -1
            acc = y(i,1:nc) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1
            y(i,1:nc) = y(i,1:nc)/val(irp(i+1)-1)
            acc  = y(i,1:nc) 
            do j=irp(i), irp(i+1)-2
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=1, nr
            acc  = y(i,1:nc) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, nr
            y(i,1:nc) = y(i,1:nc)/val(irp(i))
            acc    = y(i,1:nc) 
            do j=irp(i)+1, irp(i+1)-1
              jc      = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - val(j)*acc 
            end do
          end do
        end if

      end if

    else if (ctra) then 

      do i=1, nr
        y(i,1:nc) = x(i,1:nc)
      end do

      if (lower) then 
        if (unit) then 
          do i=nr, 1, -1
            acc = y(i,1:nc) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - (val(j))*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=nr, 1, -1
            y(i,1:nc) = y(i,1:nc)/(val(irp(i+1)-1))
            acc  = y(i,1:nc) 
            do j=irp(i), irp(i+1)-2
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - (val(j))*acc 
            end do
          end do
        end if
      else if (.not.lower) then 

        if (unit) then 
          do i=1, nr
            acc  = y(i,1:nc) 
            do j=irp(i), irp(i+1)-1
              jc    = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - (val(j))*acc 
            end do
          end do
        else if (.not.unit) then 
          do i=1, nr
            y(i,1:nc) = y(i,1:nc)/(val(irp(i)))
            acc    = y(i,1:nc) 
            do j=irp(i)+1, irp(i+1)-1
              jc      = ja(j)
              y(jc,1:nc) = y(jc,1:nc) - (val(j))*acc 
            end do
          end do
        end if

      end if
    end if
  end subroutine inner_csrsm

end subroutine psb_d_csr_cssm

function psb_d_csr_maxval(a) result(res)
  use psb_error_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_maxval
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nnz, ir, jc, nc, info
  character(len=20)  :: name='d_csr_maxval'
  logical, parameter :: debug=.false.


  res = dzero
  nnz = a%get_nzeros()
  if (allocated(a%val)) then 
    nnz = min(nnz,size(a%val))
    res = maxval(abs(a%val(1:nnz)))
  end if
end function psb_d_csr_maxval

function psb_d_csr_csnmi(a) result(res)
  use psb_error_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csnmi
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
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

end function psb_d_csr_csnmi

function psb_d_csr_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csnm1

  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_csr_csnm1'
  logical, parameter :: debug=.false.


  res = dzero
  nnz = a%get_nzeros()
  m = a%get_nrows()
  n = a%get_ncols()
  allocate(vt(n),stat=info)
  if (info /= 0) return
  vt(:) = dzero
  do i=1, m
    do j=a%irp(i),a%irp(i+1)-1
      k = a%ja(j)
      vt(k) = vt(k) + abs(a%val(j))
    end do
  end do
  res = maxval(vt(1:n))
  deallocate(vt,stat=info)

  return

end function psb_d_csr_csnm1

subroutine psb_d_csr_rowsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_rowsum
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_nrows()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  do i = 1, a%get_nrows()
    d(i) = dzero
    do j=a%irp(i),a%irp(i+1)-1  
      d(i) = d(i) + (a%val(j))
    end do
  end do

  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_rowsum

subroutine psb_d_csr_arwsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_arwsum
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_nrows()
  if (size(d) < m) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = m
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


  do i = 1, a%get_nrows()
    d(i) = dzero
    do j=a%irp(i),a%irp(i+1)-1  
      d(i) = d(i) + abs(a%val(j))
    end do
  end do

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_arwsum

subroutine psb_d_csr_colsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_colsum
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_nrows()
  n = a%get_ncols()
  if (size(d) < n) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = n
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  d   = dzero

  do i=1, m
    do j=a%irp(i),a%irp(i+1)-1
      k = a%ja(j)
      d(k) = d(k) + (a%val(j))
    end do
  end do

  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_colsum

subroutine psb_d_csr_aclsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_aclsum
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  Integer :: err_act, info, int_err(5)
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  m = a%get_nrows()
  n = a%get_ncols()
  if (size(d) < n) then 
    info=psb_err_input_asize_small_i_
    int_err(1) = 1
    int_err(2) = size(d)
    int_err(3) = n
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  d   = dzero

  do i=1, m
    do j=a%irp(i),a%irp(i+1)-1
      k = a%ja(j)
      d(k) = d(k) + abs(a%val(j))
    end do
  end do

  return
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_aclsum

subroutine psb_d_csr_get_diag(a,d,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_get_diag
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act, mnm, i, j, k
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  mnm = min(a%get_nrows(),a%get_ncols())
  if (size(d) < mnm) then 
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
    goto 9999
  end if


  if (a%is_triangle().and.a%is_unit()) then 
    d(1:mnm) = done
  else
    do i=1, mnm
      d(i) = dzero
      do k=a%irp(i),a%irp(i+1)-1
        j=a%ja(k)
        if ((j == i) .and.(j <= mnm )) then 
          d(i) = a%val(k)
        endif
      enddo
    end do
  end if
  do i=mnm+1,size(d) 
    d(i) = czero
  end do

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_get_diag


subroutine psb_d_csr_scal(d,a,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_scal
  implicit none 
  class(psb_d_csr_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  m = a%get_nrows()
  if (size(d) < m) then 
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name,i_err=(/2,size(d),0,0,0/))
    goto 9999
  end if

  do i=1, m 
    do j = a%irp(i), a%irp(i+1) -1 
      a%val(j) = a%val(j) * d(i)
    end do
  enddo

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_scal


subroutine psb_d_csr_scals(d,a,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_scals
  implicit none 
  class(psb_d_csr_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer, intent(out)            :: info

  Integer :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)


  do i=1,a%get_nzeros()
    a%val(i) = a%val(i) * d
  enddo

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_scals




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


subroutine  psb_d_csr_reallocate_nz(nz,a) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_reallocate_nz
  implicit none 
  integer, intent(in) :: nz
  class(psb_d_csr_sparse_mat), intent(inout) :: a
  Integer :: err_act, info
  character(len=20)  :: name='d_csr_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  call psb_realloc(nz,a%ja,info)
  if (info == psb_success_) call psb_realloc(nz,a%val,info)
  if (info == psb_success_) call psb_realloc(&
       & max(nz,a%get_nrows()+1,a%get_ncols()+1),a%irp,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_alloc_dealloc_,name)
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

end subroutine psb_d_csr_reallocate_nz

subroutine psb_d_csr_mold(a,b,info) 
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_mold
  use psb_error_mod
  implicit none 
  class(psb_d_csr_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer, intent(out)                    :: info
  Integer :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_d_csr_sparse_mat :: b, stat=info)

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_csr_mold

subroutine  psb_d_csr_allocate_mnnz(m,n,a,nz) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_allocate_mnnz
  implicit none 
  integer, intent(in) :: m,n
  class(psb_d_csr_sparse_mat), intent(inout) :: a
  integer, intent(in), optional :: nz
  Integer :: err_act, info, nz_
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (m < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
    goto 9999
  endif
  if (n < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
    goto 9999
  endif
  if (present(nz)) then 
    nz_ = nz
  else
    nz_ = max(7*m,7*n,1)
  end if
  if (nz_ < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
    goto 9999
  endif

  if (info == psb_success_) call psb_realloc(m+1,a%irp,info)
  if (info == psb_success_) call psb_realloc(nz_,a%ja,info)
  if (info == psb_success_) call psb_realloc(nz_,a%val,info)
  if (info == psb_success_) then 
    a%irp=0
    call a%set_nrows(m)
    call a%set_ncols(n)
    call a%set_bld()
    call a%set_triangle(.false.)
    call a%set_unit(.false.)
    call a%set_dupl(psb_dupl_def_)
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

end subroutine psb_d_csr_allocate_mnnz


subroutine psb_d_csr_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csgetptn
  implicit none

  class(psb_d_csr_sparse_mat), intent(in) :: a
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
  info = psb_success_

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

  if ((imax<imin).or.(jmax_<jmin_)) then 
    nz = 0
    return
  end if

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

  call csr_getptn(imin,imax,jmin_,jmax_,a,nz,ia,ja,nzin_,append_,info,iren)
  
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  subroutine csr_getptn(imin,imax,jmin,jmax,a,nz,ia,ja,nzin,append,info,&
       & iren)

    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none

    class(psb_d_csr_sparse_mat), intent(in)    :: a
    integer                              :: imin,imax,jmin,jmax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    integer, intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer                              :: info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw
    integer  :: debug_level, debug_unit
    character(len=20) :: name='csr_getptn'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%get_nzeros()
    irw = imin
    lrw = min(imax,a%get_nrows())
    if (irw<0) then 
      info = psb_err_pivot_too_small_
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    nzt = a%irp(lrw+1)-a%irp(irw)
    nz = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)

    if (info /= psb_success_) return
    
    if (present(iren)) then 
      do i=irw, lrw
        do j=a%irp(i), a%irp(i+1) - 1
          if ((jmin <= a%ja(j)).and.(a%ja(j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            ia(nzin_)  = iren(i)
            ja(nzin_)  = iren(a%ja(j))
          end if
        enddo
      end do
    else
      do i=irw, lrw
        do j=a%irp(i), a%irp(i+1) - 1
          if ((jmin <= a%ja(j)).and.(a%ja(j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            ia(nzin_)  = (i)
            ja(nzin_)  = (a%ja(j))
          end if
        enddo
      end do
    end if

  end subroutine csr_getptn
  
end subroutine psb_d_csr_csgetptn


subroutine psb_d_csr_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csgetrow
  implicit none

  class(psb_d_csr_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
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
  info = psb_success_
  
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

  if ((imax<imin).or.(jmax_<jmin_)) then 
    nz = 0
    return
  end if

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

  call csr_getrow(imin,imax,jmin_,jmax_,a,nz,ia,ja,val,nzin_,append_,info,&
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  subroutine csr_getrow(imin,imax,jmin,jmax,a,nz,ia,ja,val,nzin,append,info,&
       & iren)

    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none

    class(psb_d_csr_sparse_mat), intent(in)    :: a
    integer                              :: imin,imax,jmin,jmax
    integer, intent(out)                 :: nz
    integer, allocatable, intent(inout)  :: ia(:), ja(:)
    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer, intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer                              :: info
    integer, optional                    :: iren(:)
    integer  :: nzin_, nza, idx,i,j,k, nzt, irw, lrw
    integer  :: debug_level, debug_unit
    character(len=20) :: name='coo_getrow'

    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    nza = a%get_nzeros()
    irw = imin
    lrw = min(imax,a%get_nrows())
    if (irw<0) then 
      info = psb_err_pivot_too_small_
      return
    end if

    if (append) then 
      nzin_ = nzin
    else
      nzin_ = 0
    endif

    nzt = a%irp(lrw+1)-a%irp(irw)
    nz = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,val,info)

    if (info /= psb_success_) return
    
    if (present(iren)) then 
      do i=irw, lrw
        do j=a%irp(i), a%irp(i+1) - 1
          if ((jmin <= a%ja(j)).and.(a%ja(j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            val(nzin_) = a%val(j)
            ia(nzin_)  = iren(i)
            ja(nzin_)  = iren(a%ja(j))
          end if
        enddo
      end do
    else
      do i=irw, lrw
        do j=a%irp(i), a%irp(i+1) - 1
          if ((jmin <= a%ja(j)).and.(a%ja(j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            val(nzin_) = a%val(j)
            ia(nzin_)  = (i)
            ja(nzin_)  = (a%ja(j))
          end if
        enddo
      end do
    end if

  end subroutine csr_getrow

end subroutine psb_d_csr_csgetrow

subroutine psb_d_csr_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csgetblk
  implicit none

  class(psb_d_csr_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(in)                  :: imin,imax
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  Integer :: err_act, nzin, nzout
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_csgetblk



subroutine psb_d_csr_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_csput
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_csr_csput'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)


  call psb_erractionsave(err_act)
  info = psb_success_

  if (nz <= 0) then 
    info = psb_err_iarg_neg_
    int_err(1)=1
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=2
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=3
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(val) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=4
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (nz == 0) return

  nza  = a%get_nzeros()

  if (a%is_bld()) then 
    ! Build phase should only ever be in COO
    info = psb_err_invalid_mat_state_

  else  if (a%is_upd()) then 
    call  psb_d_csr_srch_upd(nz,ia,ja,val,a,&
         & imin,imax,jmin,jmax,info,gtl)

    if (info /= psb_success_) then  

      info = psb_err_invalid_mat_state_
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


contains

  subroutine psb_d_csr_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_sort_mod
    implicit none 

    class(psb_d_csr_sparse_mat), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax
    integer, intent(in) :: ia(:),ja(:)
    real(psb_dpk_), intent(in) :: val(:)
    integer, intent(out) :: info
    integer, intent(in), optional  :: gtl(:)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nr,nc,nnz,dupl,ng
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='d_csr_srch_upd'

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
    nr  = a%get_nrows()
    nc  = a%get_ncols()

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
            if ((ir > 0).and.(ir <= nr)) then 
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
            if ((ir > 0).and.(ir <= nr)) then 
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
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 

          if ((ir > 0).and.(ir <= nr)) then 

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

      case(psb_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= nr)) then 
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

  end subroutine psb_d_csr_srch_upd

end subroutine psb_d_csr_csput


subroutine psb_d_csr_reinit(a,clear)
  use psb_error_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_reinit
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  Integer :: err_act, info
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
    if (clear_) a%val(:) = dzero
    call a%set_upd()
  else
    info = psb_err_invalid_mat_state_
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

end subroutine psb_d_csr_reinit

subroutine  psb_d_csr_trim(a)
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_trim
  implicit none 
  class(psb_d_csr_sparse_mat), intent(inout) :: a
  Integer :: err_act, info, nz, m 
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  m   = a%get_nrows()
  nz  = a%get_nzeros()
  if (info == psb_success_) call psb_realloc(m+1,a%irp,info)

  if (info == psb_success_) call psb_realloc(nz,a%ja,info)
  if (info == psb_success_) call psb_realloc(nz,a%val,info)

  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_csr_trim

subroutine psb_d_csr_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_string_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_print
  implicit none 

  integer, intent(in)               :: iout
  class(psb_d_csr_sparse_mat), intent(in) :: a   
  integer, intent(in), optional     :: iv(:)
  integer, intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer, intent(in), optional     :: ivr(:), ivc(:)

  Integer :: err_act
  character(len=20)  :: name='d_csr_print'
  logical, parameter :: debug=.false.
  character(len=*), parameter  :: datatype='real'
  character(len=80)                 :: frmtv 
  integer  :: irs,ics,i,j, nmx, ni, nr, nc, nz

  if (present(eirs)) then 
    irs = eirs
  else
    irs = 0
  endif
  if (present(eics)) then 
    ics = eics
  else
    ics = 0
  endif

  if (present(head)) then 
    write(iout,'(a)') '%%MatrixMarket matrix coordinate real general'
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
    do i=1, nr
      do j=a%irp(i),a%irp(i+1)-1 
        write(iout,frmtv) iv(i),iv(a%ja(j)),a%val(j)
      end do
    enddo
  else      
    if (present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=a%irp(i),a%irp(i+1)-1 
          write(iout,frmtv) ivr(i),(a%ja(j)),a%val(j)
        end do
      enddo
    else if (present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=a%irp(i),a%irp(i+1)-1 
          write(iout,frmtv) ivr(i),ivc(a%ja(j)),a%val(j)
        end do
      enddo
    else if (.not.present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=a%irp(i),a%irp(i+1)-1 
          write(iout,frmtv) (i),ivc(a%ja(j)),a%val(j)
        end do
      enddo
    else if (.not.present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=a%irp(i),a%irp(i+1)-1 
          write(iout,frmtv) (i),(a%ja(j)),a%val(j)
        end do
      enddo
    endif
  endif

end subroutine psb_d_csr_print


subroutine psb_d_cp_csr_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_cp_csr_from_coo
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer, intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  ! This is to have fix_coo called behind the scenes
  call tmp%cp_from_coo(b,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

end subroutine psb_d_cp_csr_from_coo



subroutine psb_d_cp_csr_to_coo(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_cp_csr_to_coo
  implicit none 

  class(psb_d_csr_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                      :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)

  do i=1, nr
    do j=a%irp(i),a%irp(i+1)-1
      b%ia(j)  = i
      b%ja(j)  = a%ja(j)
      b%val(j) = a%val(j)
    end do
  end do
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)


end subroutine psb_d_cp_csr_to_coo


subroutine psb_d_mv_csr_to_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_mv_csr_to_coo
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, nc,i,j,irw, err_act
  Integer, Parameter  :: maxtry=8
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%psb_d_base_sparse_mat%mv_from(a%psb_d_base_sparse_mat)
  call b%set_nzeros(a%get_nzeros())
  call move_alloc(a%ja,b%ja)
  call move_alloc(a%val,b%val)
  call psb_realloc(nza,b%ia,info)
  if (info /= psb_success_) return
  do i=1, nr
    do j=a%irp(i),a%irp(i+1)-1
      b%ia(j)  = i
    end do
  end do
  call a%free()
  call b%fix(info)


end subroutine psb_d_mv_csr_to_coo



subroutine psb_d_mv_csr_from_coo(a,b,info) 
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_mv_csr_from_coo
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)                        :: info

  integer, allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name='mv_from_coo'

  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  call b%fix(info)
  if (info /= psb_success_) return

  nr  = b%get_nrows()
  nc  = b%get_ncols()
  nza = b%get_nzeros()
  
  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)

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
               & 'Strange situation: i>nr ',i,nr,j,nza,irw
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


end subroutine psb_d_mv_csr_from_coo


subroutine psb_d_mv_csr_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_mv_csr_to_fmt
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout)  :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
    ! Need to fix trivial copies! 
  type is (psb_d_csr_sparse_mat) 
    call b%psb_d_base_sparse_mat%mv_from(a%psb_d_base_sparse_mat)
    call move_alloc(a%irp, b%irp)
    call move_alloc(a%ja,  b%ja)
    call move_alloc(a%val, b%val)
    call a%free()
    
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_csr_to_fmt


subroutine psb_d_cp_csr_to_fmt(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_cp_csr_to_fmt
  implicit none 

  class(psb_d_csr_sparse_mat), intent(in)   :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                       :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_


  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_d_csr_sparse_mat) 
    call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)
    call psb_safe_cpy( a%irp, b%irp , info)
    call psb_safe_cpy( a%ja , b%ja  , info)
    call psb_safe_cpy( a%val, b%val , info)

  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_cp_csr_to_fmt


subroutine psb_d_mv_csr_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_mv_csr_from_fmt
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)                         :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nza, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  type is (psb_d_csr_sparse_mat) 
    call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
    call move_alloc(b%irp, a%irp)
    call move_alloc(b%ja,  a%ja)
    call move_alloc(b%val, a%val)
    call b%free()

  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_csr_from_fmt



subroutine psb_d_cp_csr_from_fmt(a,b,info) 
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_cp_csr_from_fmt
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer, intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  Integer             :: nz, nr, i,j,irw, err_act, nc
  Integer, Parameter  :: maxtry=8
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_csr_sparse_mat) 
    call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
    call psb_safe_cpy( b%irp, a%irp , info)
    call psb_safe_cpy( b%ja , a%ja  , info)
    call psb_safe_cpy( b%val, a%val , info)

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_csr_from_fmt


subroutine psb_d_csr_cp_from(a,b)
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_cp_from
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout) :: a
  type(psb_d_csr_sparse_mat), intent(in)   :: b


  Integer :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_success_

  call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
  call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
  call psb_safe_cpy( b%irp, a%irp , info)
  call psb_safe_cpy( b%ja , a%ja  , info)
  call psb_safe_cpy( b%val, a%val , info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_csr_cp_from

subroutine psb_d_csr_mv_from(a,b)
  use psb_error_mod
  use psb_d_csr_mat_mod, psb_protect_name => psb_d_csr_mv_from
  implicit none 

  class(psb_d_csr_sparse_mat), intent(inout)  :: a
  type(psb_d_csr_sparse_mat), intent(inout) :: b


  Integer :: err_act, info
  character(len=20)  :: name='mv_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
  call move_alloc(b%irp, a%irp)
  call move_alloc(b%ja,  a%ja)
  call move_alloc(b%val, a%val)
  call b%free()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_d_csr_mv_from


