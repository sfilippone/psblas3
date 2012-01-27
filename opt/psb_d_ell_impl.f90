
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

subroutine psb_d_ell_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csmv
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_csmv'
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
    info = 36
    call psb_errpush(info,name,i_err=(/3,n,0,0,0/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
    goto 9999
  end if


  call psb_d_ell_csmv_inner(m,n,alpha,size(a%ja,2),&
       & a%ja,size(a%ja,1),a%val,size(a%val,1),&
       & a%is_triangle(),a%is_unit(),&
       & x,beta,y,tra) 

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

  subroutine psb_d_ell_csmv_inner(m,n,alpha,nc,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,tra) 
    integer(psb_ipk_), intent(in)             :: m,n,nc,ldj,ldv,ja(ldj,*)
    real(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    real(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)             :: is_triangle,is_unit,tra


    integer(psb_ipk_) :: i,j,k, ir, jc, m4
    real(psb_dpk_) :: acc(4) 


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

        m4 = mod(m,4) 
        do i=1,m4
          acc(1) = dzero 
          do j=1,nc
            acc(1)  = acc(1) + val(i,j) * x(ja(i,j))          
          enddo
          y(i) = alpha*acc(1)
        end do


        if (alpha == done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4)  = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = acc(1:4)
          end do

        else if (alpha == -done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) - val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = acc(1:4)
          end do

        else 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4)  = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = alpha * acc(1:4)
          end do

        end if


      else if (beta == done) then 
        

        m4 = mod(m,4) 
        do i=1,m4
          acc(1) = dzero 
          do j=1,nc
            acc(1)  = acc(1) + val(i,j) * x(ja(i,j))          
          enddo
          y(i) = y(i) + alpha*acc(1)
        end do

        if (alpha == done) then 
          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = y(i:i+3) + acc(1:4)
          end do

        else if (alpha == -done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) - val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = y(i:i+3) + acc(1:4)
          end do

        else 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = y(i:i+3) + alpha*acc(1:4)
          end do

        end if

      else if (beta == -done) then 

        m4 = mod(m,4) 
        do i=1,m4
          acc(1) = dzero 
          do j=1,nc
            acc(1)  = acc(1) + val(i,j) * x(ja(i,j))          
          enddo
          y(i) = - y(i) + alpha*acc(1)
        end do

        if (alpha == done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = -y(i:i+3) + acc(1:4)
          end do

        else if (alpha == -done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) - val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = -y(i:i+3) + acc(1:4)
          end do

        else 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = -y(i:i+3) + alpha*acc(1:4)
          end do

        end if

      else 

        m4 = mod(m,4) 
        do i=1,m4
          acc(1) = dzero 
          do j=1,nc
            acc(1)  = acc(1) + val(i,j) * x(ja(i,j))          
          enddo
          y(i) = beta*y(i) + alpha*acc(1)
        end do

        if (alpha == done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = beta*y(i:i+3) + acc(1:4)
          end do

        else if (alpha == -done) then 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) - val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = beta*y(i:i+3) + acc(1:4)
          end do

        else 

          do i=m4+1,m,4
            acc  = dzero
            do j=1,nc
              acc(1:4) = acc(1:4) + val(i:i+3,j) * x(ja(i:i+3,j))          
            enddo
            y(i:i+3) = beta*y(i:i+3) + alpha*acc(1:4)
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

      !
      ! Need to think about this.
      ! Transpose does not mix well with ELLPACK.
      ! 
      if (alpha == done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir) = y(ir) +  val(i,j)*x(i)
          end do
        enddo

      else if (alpha == -done) then

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir) = y(ir) -  val(i,j)*x(i)
          end do
        enddo

      else                    

        do i=1,n
          do j=1,nc
            ir = ja(i,j)
            y(ir) = y(ir) +  alpha*val(i,j)*x(i)
          end do
        enddo

      end if

    endif

    if (is_triangle.and.is_unit) then 
      do i=1, min(m,n)
        y(i) = y(i) + alpha*x(i)
      end do
    end if


  end subroutine psb_d_ell_csmv_inner


end subroutine psb_d_ell_csmv

subroutine psb_d_ell_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csmm
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nxy
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
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
  tra = (psb_toupper(trans_) == 'T').or.(psb_toupper(trans_)=='C')

  if (tra) then 
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

  nxy = min(size(x,2) , size(y,2) )

  allocate(acc(nxy), stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='allocate')
    goto 9999
  end if

  call  psb_d_ell_csmm_inner(m,n,nxy,alpha,size(a%ja,2),&
       & a%ja,size(a%ja,1),a%val,size(a%val,1), &
       & a%is_triangle(),a%is_unit(),x,size(x,1), &
       & beta,y,size(y,1),tra,acc) 


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
  subroutine psb_d_ell_csmm_inner(m,n,nxy,alpha,nc,ja,ldj,val,ldv,&
       & is_triangle,is_unit,x,ldx,beta,y,ldy,tra,acc) 
    integer(psb_ipk_), intent(in)             :: m,n,ldx,ldy,nxy,nc,ldj,ldv
    integer(psb_ipk_), intent(in)             :: ja(ldj,*)
    real(psb_dpk_), intent(in)      :: alpha, beta, x(ldx,*),val(ldv,*)
    real(psb_dpk_), intent(inout)   :: y(ldy,*)
    logical, intent(in)             :: is_triangle,is_unit,tra

    real(psb_dpk_), intent(inout)   :: acc(*)
    integer(psb_ipk_) :: i,j,k, ir, jc


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

    if (.not.tra) then 

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

    endif

    if (is_triangle.and.is_unit) then 
      do i=1, min(m,n)
        y(i,1:nxy) = y(i,1:nxy) + alpha*x(i,1:nxy)
      end do
    end if

  end subroutine psb_d_ell_csmm_inner

end subroutine psb_d_ell_csmm


subroutine psb_d_ell_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_cssv
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_cssv'
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

    call inner_ellsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & size(a%ja,2),a%irn,a%idiag,a%ja,size(a%ja,1),a%val,size(a%val,1),&
         & x,y,info) 

    if (info /= 0) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif

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
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    call inner_ellsv(tra,a%is_lower(),a%is_unit(),a%get_nrows(),&
         & size(a%ja,2),a%irn,a%idiag,a%ja,size(a%ja,1),a%val,size(a%val,1),&
         & x,tmp,info) 

    if (info == 0) &
         & call psb_geaxpby(m,alpha,tmp,beta,y,info)

    if (info /= 0) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif

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

  subroutine inner_ellsv(tra,lower,unit,n,nc,irn,idiag,ja,ldj,val,ldv,x,y,info) 
    implicit none 
    logical, intent(in)                 :: tra,lower,unit
    integer(psb_ipk_), intent(in)                 :: n,nc,ldj,ldv
    integer(psb_ipk_), intent(in)                 :: irn(*),idiag(*), ja(ldj,*)
    real(psb_dpk_), intent(in)          :: val(ldv,*)
    real(psb_dpk_), intent(in)          :: x(*)
    real(psb_dpk_), intent(out)         :: y(*)
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_ipk_) :: i,j,k,m, ir, jc
    real(psb_dpk_) :: acc

    !
    ! The only error condition here is if
    ! the matrix is non-unit and some idiag value is illegal.
    !
    info = 0 
    if (.not.tra) then 

      if (lower) then 

        if (unit) then 
          do i=1, n
            acc = dzero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j))
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = dzero 
            do j=1,idiag(i)-1
              acc = acc + val(i,j)*y(ja(i,j))
            end do
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i) = (x(i) - acc)/val(i,idiag(i))
          end do
        end if

      else if (.not.lower) then 

        if (unit) then 

          do i=n, 1, -1 
            acc = dzero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j))
            end do
            y(i) = x(i) - acc
          end do

        else if (.not.unit) then 

          do i=n, 1, -1 
            acc = dzero 
            do j=idiag(i)+1, irn(i)
              acc = acc + val(i,j)*y(ja(i,j))
            end do
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i) = (x(i) - acc)/val(i,idiag(i))
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

            do j=1,irn(i)
              jc    = ja(i,j)
              y(jc) = y(jc) - val(i,j)*acc 
            end do

          end do

        else if (.not.unit) then 

          do i=n, 1, -1
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i) = y(i)/val(i,idiag(i))
            acc  = y(i) 
            do j=1,idiag(i)
              jc    = ja(i,j)
              y(jc) = y(jc) - val(i,j)*acc 
            end do
          end do

        end if

      else if (.not.lower) then 

        if (unit) then 

          do i=1, n
            acc  = y(i) 
            do j=1, irn(i)
              jc    = ja(i,j)
              y(jc) = y(jc) - val(i,j)*acc 
            end do
          end do

        else if (.not.unit) then 

          do i=1, n
            if (idiag(i) <= 0) then 
              info = -1 
              return
            endif
            y(i) = y(i)/val(i,idiag(i))
            acc  = y(i) 
            do j=idiag(i)+1, irn(i) 
              jc    = ja(i,j)
              y(jc) = y(jc) - val(i,j)*acc 
            end do
          end do

        end if

      end if
    end if
  end subroutine inner_ellsv

end subroutine psb_d_ell_cssv



subroutine psb_d_ell_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_cssm
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nxy
  real(psb_dpk_), allocatable :: tmp(:,:), acc(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_cssm'
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
    info = 36
    call psb_errpush(info,name,i_err=(/3,m,0,0,0/))
    goto 9999
  end if

  if (size(y,1)<m) then 
    info = 36
    call psb_errpush(info,name,i_err=(/5,m,0,0,0/))
    goto 9999
  end if

  nxy = min(size(x,2),size(y,2))

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
    allocate(acc(nxy), stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    call inner_ellsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nxy,&
         & size(a%ja,2),a%irn,a%idiag,a%ja,size(a%ja,1),a%val,size(a%val,1),&
         & x,size(x,1),y,size(y,1),acc,info) 

    if (info /= 0) then 
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif

    if (alpha == done) then 
      ! do nothing
    else if (alpha == -done) then 
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

    call inner_ellsm(tra,a%is_lower(),a%is_unit(),a%get_nrows(),nxy,&
         & size(a%ja,2),a%irn,a%idiag,a%ja,size(a%ja,1),a%val,size(a%val,1),&
         & x,size(x,1),tmp,size(tmp,1),acc,info) 

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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains 

  subroutine inner_ellsm(tra,lower,unit,n,nc,nxy,irn,idiag,ja,ldj,val,ldv,&
       & x,ldx,y,ldy,acc,info) 
    implicit none 
    logical, intent(in)                 :: tra,lower,unit
    integer(psb_ipk_), intent(in)                 :: n,nc,ldj,ldv,nxy,ldx,ldy
    integer(psb_ipk_), intent(in)                 :: irn(*),idiag(*), ja(ldj,*)
    real(psb_dpk_), intent(in)          :: val(ldv,*)
    real(psb_dpk_), intent(in)          :: x(ldx,nxy)
    real(psb_dpk_), intent(out)         :: y(ldy,nxy), acc(nxy)
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_ipk_) :: i,j,k,m, ir, jc

    !
    ! The only error condition here is if
    ! the matrix is non-unit and some idiag value is illegal.
    !
    info = 0 

    if (.not.tra) then 

      if (lower) then 

        if (unit) then 
          do i=1, n
            acc = dzero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do
        else if (.not.unit) then 
          do i=1, n
            acc = dzero 
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
            acc = dzero 
            do j=1,irn(i)
              acc = acc + val(i,j)*y(ja(i,j),:)
            end do
            y(i,:) = x(i,:) - acc
          end do

        else if (.not.unit) then 

          do i=n, 1, -1 
            acc = dzero 
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
    end if
  end subroutine inner_ellsm

end subroutine psb_d_ell_cssm


function psb_d_ell_csnmi(a) result(res)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csnmi
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nr, ir, jc, nc
  real(psb_dpk_) :: acc
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_csnmi'
  logical, parameter :: debug=.false.


  res = dzero 
 
  do i = 1, a%get_nrows()
    acc = sum(abs(a%val(i,:)))
    res = max(res,acc)
  end do

end function psb_d_ell_csnmi


function psb_d_ell_csnm1(a) result(res)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csnm1

  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_csnm1'
  logical, parameter :: debug=.false.


  res = -done 
  nnz = a%get_nzeros()
  m = a%get_nrows()
  n = a%get_ncols()
  allocate(vt(n),stat=info)
  if (info /= 0) return
  vt(:) = dzero
  do i=1, m
    do j=1,a%irn(i)
      k = a%ja(i,j)
      vt(k) = vt(k) + abs(a%val(i,j))
    end do
  end do
  res = maxval(vt(1:n))
  deallocate(vt,stat=info)

  return

end function psb_d_ell_csnm1


subroutine psb_d_ell_rowsum(d,a) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_rowsum
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
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
    do j=1,a%irn(i) 
      d(i) = d(i) + (a%val(i,j))
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

end subroutine psb_d_ell_rowsum

subroutine psb_d_ell_arwsum(d,a) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_arwsum
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
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
    do j=1,a%irn(i)
      d(i) = d(i) + abs(a%val(i,j))
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

end subroutine psb_d_ell_arwsum

subroutine psb_d_ell_colsum(d,a) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_colsum
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
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
    do j=1,a%irn(i)
      k = a%ja(i,j)
      d(k) = d(k) + (a%val(i,j))
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

end subroutine psb_d_ell_colsum

subroutine psb_d_ell_aclsum(d,a) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_aclsum
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
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
    do j=1,a%irn(i)
      k = a%ja(i,j)
      d(k) = d(k) + abs(a%val(i,j))
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

end subroutine psb_d_ell_aclsum


subroutine psb_d_ell_get_diag(a,d,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_get_diag
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act, mnm, i, j, k
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
      if (1<=a%idiag(i).and.(a%idiag(i)<=size(a%ja,2))) &
           &  d(i) = a%val(i,a%idiag(i))
    end do
  end if
  do i=mnm+1,size(d) 
    d(i) = dzero
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

end subroutine psb_d_ell_get_diag


subroutine psb_d_ell_scal(d,a,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_scal
  implicit none 
  class(psb_d_ell_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act,mnm, i, j, m
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
    a%val(i,:) = a%val(i,:) * d(i)
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

end subroutine psb_d_ell_scal


subroutine psb_d_ell_scals(d,a,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_scals
  implicit none 
  class(psb_d_ell_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)


  a%val(:,:) = a%val(:,:) * d

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_d_ell_scals




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


subroutine  psb_d_ell_reallocate_nz(nz,a) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_reallocate_nz
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_d_ell_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: m, nzrm
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='d_ell_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  !
  ! What should this really do??? 
  ! 
  m    = a%get_nrows()
  nzrm = (nz+m-1)/m
  call psb_realloc(m,nzrm,a%ja,info)
  if (info == psb_success_) call psb_realloc(m,nzrm,a%val,info)
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

end subroutine psb_d_ell_reallocate_nz

subroutine psb_d_ell_mold(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_mold
  implicit none 
  class(psb_d_ell_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_d_ell_sparse_mat :: b, stat=info)

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

end subroutine psb_d_ell_mold

subroutine  psb_d_ell_allocate_mnnz(m,n,a,nz) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_allocate_mnnz
  implicit none 
  integer(psb_ipk_), intent(in) :: m,n
  class(psb_d_ell_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(in), optional :: nz
  integer(psb_ipk_) :: err_act, info, nz_
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
    nz_ = (nz + m -1 )/m
  else
    nz_ = (max(7*m,7*n,1)+m-1)/m
  end if
  if (nz_ < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
    goto 9999
  endif

  if (info == psb_success_) call psb_realloc(m,a%irn,info)
  if (info == psb_success_) call psb_realloc(m,a%idiag,info)
  if (info == psb_success_) call psb_realloc(m,nz_,a%ja,info)
  if (info == psb_success_) call psb_realloc(m,nz_,a%val,info)
  if (info == psb_success_) then 
    a%irn   = 0
    a%idiag = 0
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

end subroutine psb_d_ell_allocate_mnnz


subroutine psb_d_ell_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csgetptn
  implicit none

  class(psb_d_ell_sparse_mat), intent(in) :: a
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

  call ell_getptn(imin,imax,jmin_,jmax_,a,nz,ia,ja,nzin_,append_,info,iren)
  
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

  subroutine ell_getptn(imin,imax,jmin,jmax,a,nz,ia,ja,nzin,append,info,&
       & iren)

    implicit none

    class(psb_d_ell_sparse_mat), intent(in)    :: a
    integer(psb_ipk_) :: imin,imax,jmin,jmax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    integer(psb_ipk_), intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional                    :: iren(:)
    integer(psb_ipk_) :: nzin_, nza, idx,i,j,k, nzt, irw, lrw
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20) :: name='ell_getptn'

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

    nzt = sum(a%irn(irw:lrw))
    nz  = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)

    if (info /= psb_success_) return
    
    if (present(iren)) then 
      do i=irw, lrw
        do j=1,a%irn(i)
          if ((jmin <= a%ja(i,j)).and.(a%ja(i,j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            ia(nzin_)  = iren(i)
            ja(nzin_)  = iren(a%ja(i,j))
          end if
        enddo
      end do
    else
      do i=irw, lrw
        do j=1,a%irn(i)
          if ((jmin <= a%ja(i,j)).and.(a%ja(i,j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            ia(nzin_)  = (i)
            ja(nzin_)  = (a%ja(i,j))
          end if
        enddo
      end do
    end if

  end subroutine ell_getptn
  
end subroutine psb_d_ell_csgetptn


subroutine psb_d_ell_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csgetrow
  implicit none

  class(psb_d_ell_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i
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

  call ell_getrow(imin,imax,jmin_,jmax_,a,nz,ia,ja,val,nzin_,append_,info,&
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

  subroutine ell_getrow(imin,imax,jmin,jmax,a,nz,ia,ja,val,nzin,append,info,&
       & iren)

    implicit none

    class(psb_d_ell_sparse_mat), intent(in)    :: a
    integer(psb_ipk_) :: imin,imax,jmin,jmax
    integer(psb_ipk_), intent(out)                 :: nz
    integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
    real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
    integer(psb_ipk_), intent(in)                  :: nzin
    logical, intent(in)                  :: append
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional                    :: iren(:)
    integer(psb_ipk_) :: nzin_, nza, idx,i,j,k, nzt, irw, lrw
    integer(psb_ipk_) :: debug_level, debug_unit
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

    nzt = sum(a%irn(irw:lrw))
    nz  = 0 


    call psb_ensure_size(nzin_+nzt,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,ja,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nzt,val,info)

    if (info /= psb_success_) return
    
    if (present(iren)) then 
      do i=irw, lrw
        do j=1,a%irn(i)
          if ((jmin <= a%ja(i,j)).and.(a%ja(i,j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            val(nzin_) = a%val(i,j)
            ia(nzin_)  = iren(i)
            ja(nzin_)  = iren(a%ja(i,j))
          end if
        enddo
      end do
    else
      do i=irw, lrw
        do j=1,a%irn(i)
          if ((jmin <= a%ja(i,j)).and.(a%ja(i,j)<=jmax)) then 
            nzin_ = nzin_ + 1
            nz    = nz + 1
            val(nzin_) = a%val(i,j)
            ia(nzin_)  = (i)
            ja(nzin_)  = (a%ja(i,j))
          end if
        enddo
      end do
    end if

  end subroutine ell_getrow

end subroutine psb_d_ell_csgetrow

subroutine psb_d_ell_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csgetblk
  implicit none

  class(psb_d_ell_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act, nzin, nzout
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

end subroutine psb_d_ell_csgetblk



subroutine psb_d_ell_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csput
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)


  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_csput'
  logical, parameter :: debug=.false.
  integer(psb_ipk_) :: nza, i,j,k, nzl, isza, int_err(5)


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
    call  psb_d_ell_srch_upd(nz,ia,ja,val,a,&
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

  subroutine psb_d_ell_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    implicit none 

    class(psb_d_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in) :: nz, imin,imax,jmin,jmax
    integer(psb_ipk_), intent(in) :: ia(:),ja(:)
    real(psb_dpk_), intent(in) :: val(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in), optional  :: gtl(:)
    integer(psb_ipk_) :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nr,nc,nnz,dupl,ng
    integer(psb_ipk_) :: debug_level, debug_unit
    character(len=20)    :: name='d_ell_srch_upd'

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
              nc = a%irn(ir)
              ip = psb_ibsrch(ic,nc,a%ja(i,1:nc))    
              if (ip>0) then 
                a%val(i,ip) = val(i)
              else
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ic,' in: ',nc,&
                     & ' : ',a%ja(i,1:nc)
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
              nc = a%irn(ir)
              ip = psb_ibsrch(ic,nc,a%ja(i,1:nc))    
              if (ip>0) then 
                a%val(i,ip) = a%val(i,ip) + val(i)
              else
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ic,' in: ',nc,&
                     & ' : ',a%ja(i,1:nc)
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

            nc = a%irn(ir)
            ip = psb_ibsrch(ic,nc,a%ja(i,1:nc))    
            if (ip>0) then 
              a%val(i,ip) = val(i)
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Was searching ',ic,' in: ',nc,&
                   & ' : ',a%ja(i,1:nc)
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
            nc = a%irn(ir)
            ip = psb_ibsrch(ic,nc,a%ja(i,1:nc))    
            if (ip>0) then 
              a%val(i,ip) = a%val(i,ip) + val(i)
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

  end subroutine psb_d_ell_srch_upd

end subroutine psb_d_ell_csput



subroutine psb_d_ell_reinit(a,clear)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_reinit
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  integer(psb_ipk_) :: err_act, info
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
    if (clear_) a%val(:,:) = dzero
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

end subroutine psb_d_ell_reinit

subroutine  psb_d_ell_trim(a)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_trim
  implicit none 
  class(psb_d_ell_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info, nz, m, nzm
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  m   = a%get_nrows()
  nzm = maxval(a%irn(1:m))
  
  call psb_realloc(m,a%irn,info)
  if (info == psb_success_) call psb_realloc(m,a%idiag,info)
  if (info == psb_success_) call psb_realloc(m,nzm,a%ja,info)
  if (info == psb_success_) call psb_realloc(m,nzm,a%val,info)

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

end subroutine psb_d_ell_trim

subroutine psb_d_ell_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_print
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_d_ell_sparse_mat), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  integer(psb_ipk_), intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_ell_print'
  logical, parameter :: debug=.false.

  character(len=80)                 :: frmtv 
  integer(psb_ipk_) :: irs,ics,i,j, nmx, ni, nr, nc, nz

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

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nz  = a%get_nzeros()
  nmx = max(nr,nc,1)
  ni  = floor(log10(1.0*nmx)) + 1

  write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
  write(iout,*) nr, nc, nz 
  if(present(iv)) then 
    do i=1, nr
      do j=1,a%irn(i)
        write(iout,frmtv) iv(i),iv(a%ja(i,j)),a%val(i,j)
      end do
    enddo
  else      
    if (present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmtv) ivr(i),(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmtv) ivr(i),ivc(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (.not.present(ivr).and.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmtv) (i),ivc(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    else if (.not.present(ivr).and..not.present(ivc)) then 
      do i=1, nr
        do j=1,a%irn(i)
          write(iout,frmtv) (i),(a%ja(i,j)),a%val(i,j)
        end do
      enddo
    endif
  endif

end subroutine psb_d_ell_print


subroutine psb_d_cp_ell_from_coo(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_cp_ell_from_coo
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

   info = psb_success_
   ! This is to have fix_coo called behind the scenes
   call b%cp_to_coo(tmp,info)
   if (info == psb_success_) call a%mv_from_coo(tmp,info)

end subroutine psb_d_cp_ell_from_coo



subroutine psb_d_cp_ell_to_coo(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_cp_ell_to_coo
  implicit none 

  class(psb_d_ell_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                      :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,k,irw, idl,err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)

  k=0
  do i=1, nr
    do j=1,a%irn(i)
      k = k + 1
      b%ia(k)  = i
      b%ja(k)  = a%ja(i,j)
      b%val(k) = a%val(i,j)
    end do
  end do
  call b%set_nzeros(a%get_nzeros())
  call b%fix(info)

end subroutine psb_d_cp_ell_to_coo

subroutine psb_d_mv_ell_to_coo(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_mv_ell_to_coo
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                       :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,k,irw, idl,err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  ! Taking  a path slightly slower but with less memory footprint
  deallocate(a%idiag)
  call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)

  call psb_realloc(nza,b%ia,info)
  if (info == 0)   call psb_realloc(nza,b%ja,info)
  if (info /= 0) goto 9999
  k=0
  do i=1, nr
    do j=1,a%irn(i)
      k = k + 1
      b%ia(k)  = i
      b%ja(k)  = a%ja(i,j)
    end do
  end do
  deallocate(a%ja, stat=info)

  if (info == 0) call psb_realloc(nza,b%val,info)
  if (info /= 0) goto 9999

  k=0
  do i=1, nr
    do j=1,a%irn(i)
      k = k + 1
      b%val(k)  = a%val(i,j)
    end do
  end do
  call a%free()
  call b%set_nzeros(nza)
  call b%fix(info)
  return

9999 continue
  info = psb_err_alloc_dealloc_
  return
end subroutine psb_d_mv_ell_to_coo


subroutine psb_d_mv_ell_from_coo(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_mv_ell_from_coo
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                        :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,k, idl,err_act, nc, nzm, ir, ic
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  call b%fix(info)
  if (info /= psb_success_) return

  nr  = b%get_nrows()
  nc  = b%get_ncols()
  nza = b%get_nzeros()
  if (b%is_sorted()) then 
    ! If it is sorted then we can lessen memory impact 
    call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)

    ! First compute the number of nonzeros in each row.
    call psb_realloc(nr,a%irn,info) 
    if (info /= 0) goto 9999
    a%irn = 0
    do i=1, nza
      a%irn(b%ia(i)) = a%irn(b%ia(i)) + 1
    end do
    nzm = 0 
    do i=1, nr
      nzm = max(nzm,a%irn(i))
      a%irn(i) = 0
    end do
    ! Second: copy the column indices.
    call psb_realloc(nr,a%idiag,info) 
    if (info == 0) call psb_realloc(nr,nzm,a%ja,info) 
    if (info /= 0) goto 9999
    do i=1, nza
      ir = b%ia(i)
      ic = b%ja(i)
      j  = a%irn(ir) + 1 
      a%ja(ir,j) = ic
      a%irn(ir)  = j
    end do
    ! Third copy the other stuff
    deallocate(b%ia,b%ja,stat=info) 
    if (info == 0) call psb_realloc(nr,a%idiag,info)
    if (info == 0) call psb_realloc(nr,nzm,a%val,info)
    if (info /= 0) goto 9999
    k = 0 
    do i=1, nr
      a%idiag(i) = 0 
      do j=1, a%irn(i)
        k = k + 1 
        a%val(i,j) = b%val(k)
        if (i==a%ja(i,j)) a%idiag(i)=j
      end do
      do j=a%irn(i)+1, nzm
        a%ja(i,j) = i
        a%val(i,j) = dzero
      end do
    end do

  else 
    ! If b is not sorted, the only way is to copy. 
    call a%cp_from_coo(b,info)
    if (info /= 0) goto 9999
  end if

  call b%free()

  return

9999 continue
  info = psb_err_alloc_dealloc_
  return


end subroutine psb_d_mv_ell_from_coo


subroutine psb_d_mv_ell_to_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_mv_ell_to_fmt
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout)  :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_to_coo(b,info)
    ! Need to fix trivial copies! 
  type is (psb_d_ell_sparse_mat) 
    call b%psb_d_base_sparse_mat%mv_from(a%psb_d_base_sparse_mat)
    call move_alloc(a%irn,   b%irn)
    call move_alloc(a%idiag, b%idiag)
    call move_alloc(a%ja,    b%ja)
    call move_alloc(a%val,   b%val)
    call a%free()
    
  class default
    call a%mv_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_ell_to_fmt


subroutine psb_d_cp_ell_to_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_cp_ell_to_fmt
  implicit none 

  class(psb_d_ell_sparse_mat), intent(in)   :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                       :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_


  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_to_coo(b,info)

  type is (psb_d_ell_sparse_mat) 
    call b%psb_d_base_sparse_mat%cp_from(a%psb_d_base_sparse_mat)
    call psb_safe_cpy( a%idiag, b%idiag , info)
    call psb_safe_cpy( a%irn,   b%irn , info)
    call psb_safe_cpy( a%ja ,   b%ja  , info)
    call psb_safe_cpy( a%val,   b%val , info)

  class default
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_cp_ell_to_fmt


subroutine psb_d_mv_ell_from_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_mv_ell_from_fmt
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                         :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  type is (psb_d_ell_sparse_mat) 
    call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
    call move_alloc(b%irn,   a%irn)
    call move_alloc(b%idiag, a%idiag)
    call move_alloc(b%ja,    a%ja)
    call move_alloc(b%val,   a%val)
    call b%free()

  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

end subroutine psb_d_mv_ell_from_fmt



subroutine psb_d_cp_ell_from_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_cp_ell_from_fmt
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nz, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  select type (b)
  type is (psb_d_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

  type is (psb_d_ell_sparse_mat) 
    call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
    call psb_safe_cpy( b%irn,   a%irn ,  info)
    call psb_safe_cpy( b%idiag, a%idiag, info)
    call psb_safe_cpy( b%ja ,   a%ja  ,  info)
    call psb_safe_cpy( b%val,   a%val ,  info)

  class default
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
end subroutine psb_d_cp_ell_from_fmt


subroutine psb_d_ell_cp_from(a,b)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_cp_from
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  type(psb_d_ell_sparse_mat), intent(in)   :: b


  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  info = psb_success_

  call a%allocate(b%get_nrows(),b%get_ncols(),b%get_nzeros())
  call a%psb_d_base_sparse_mat%cp_from(b%psb_d_base_sparse_mat)
  call psb_safe_cpy( b%irn,   a%irn ,  info)
  call psb_safe_cpy( b%idiag, a%idiag, info)
  call psb_safe_cpy( b%ja ,   a%ja  ,  info)
  call psb_safe_cpy( b%val,   a%val ,  info)

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

end subroutine psb_d_ell_cp_from

subroutine psb_d_ell_mv_from(a,b)
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_mv_from
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout)  :: a
  type(psb_d_ell_sparse_mat), intent(inout) :: b


  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='mv_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%psb_d_base_sparse_mat%mv_from(b%psb_d_base_sparse_mat)
  call move_alloc(b%idiag, a%idiag)
  call move_alloc(b%irn,   a%irn)
  call move_alloc(b%ja,    a%ja)
  call move_alloc(b%val,   a%val)
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

end subroutine psb_d_ell_mv_from


