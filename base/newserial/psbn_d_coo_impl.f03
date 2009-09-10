
subroutine d_coo_cssm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_cssm_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
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

  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  if (.not. (a%is_triangle())) then 
    write(0,*) 'Called SM on a non-triangular mat!'
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
  tra = ((trans_=='T').or.(trans_=='t'))
  m   = a%get_nrows()
  nc  = min(size(x,2) , size(y,2)) 


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
    call inner_coosm(tra,a,x,y,info)
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
    call inner_coosm(tra,a,tmp,y,info)
    do  i = 1, m
      y(i,:) = alpha*tmp(i,:) + beta*y(i,:)
    end do
  end if

  if(info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='inner_coosm')
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

  subroutine inner_coosm(tra,a,x,y,info) 
    implicit none 
    logical, intent(in)                 :: tra  
    class(psbn_d_coo_sparse_mat), intent(in) :: a
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


    if (.not.a%is_sorted()) then 
      info = 1121
      return
    end if

    nnz = a%get_nzeros()

    if (.not.tra) then 

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          j = 1
          do i=1, a%get_nrows()
            acc = dzero
            do 
              if (j > nnz) exit
              if (a%ia(j) > i) exit
              acc = acc + a%val(j)*x(a%ja(j),:)
              j   = j + 1
            end do
            y(i,:) = x(i,:) - acc
          end do
        else if (.not.a%is_unit()) then 
          j = 1
          do i=1, a%get_nrows()
            acc = dzero
            do 
              if (j > nnz) exit
              if (a%ia(j) > i) exit
              if (a%ja(j) == i) then 
                y(i,:) = (x(i,:) - acc)/a%val(j)
                j = j + 1
                exit
              end if
              acc = acc + a%val(j)*x(a%ja(j),:)
              j   = j + 1
            end do
          end do
        end if

      else if (a%is_upper()) then 
        if (a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              acc = acc + a%val(j)*x(a%ja(j),:)
              j   = j - 1
            end do
            y(i,:) = x(i,:) - acc
          end do

        else if (.not.a%is_unit()) then 

          j = nnz
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              if (a%ja(j) == i) then 
                y(i,:) = (x(i,:) - acc)/a%val(j)
                j = j - 1
                exit
              end if
              acc = acc + a%val(j)*x(a%ja(j),:)
              j   = j - 1
            end do
          end do
        end if

      end if

    else if (tra) then 

      do i=1, a%get_nrows()
        y(i,:) = x(i,:)
      end do

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1
            acc = y(i,:) 
            do
              if (j < 1) exit
              if (a%ia(j) < i) exit
              jc    = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
              j     = j - 1 
            end do
          end do
        else if (.not.a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1
            if (a%ja(j) == i) then 
              y(i,:) = y(i,:) /a%val(j)
              j    = j - 1
            end if
            acc  = y(i,:) 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              jc    = a%ja(j)
              y(jc,:) = y(jc,:) - a%val(j)*acc 
              j     = j - 1
            end do
          end do

        else if (a%is_upper()) then 
          if (a%is_unit()) then 
            j = 1
            do i=1, a%get_nrows()
              acc = y(i,:)
              do 
                if (j > nnz) exit
                if (a%ia(j) > i) exit
                jc    = a%ja(j)
                y(jc,:) = y(jc,:) - a%val(j)*acc 
                j   = j + 1
              end do
            end do
          else if (.not.a%is_unit()) then 
            j = 1
            do i=1, a%get_nrows()
              if (a%ja(j) == i) then 
                y(i,:) = y(i,:) /a%val(j)
                j    = j + 1
              end if
              acc = y(i,:)
              do 
                if (j > nnz) exit
                if (a%ia(j) > i) exit
                jc    = a%ja(j)
                y(jc,:) = y(jc,:) - a%val(j)*acc 
                j   = j + 1
              end do
            end do
          end if
        end if
      end if
    end if
  end subroutine inner_coosm

end subroutine d_coo_cssm_impl



subroutine d_coo_cssv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_const_mod
  use psb_error_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_cssv_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
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
  character(len=20)  :: name='d_coo_cssv_impl'
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
    call inner_coosv(tra,a,x,y,info)
    if (info /= 0) then 
      call psb_errpush(info,name)
      goto 9999
    end if
    do  i = 1, m
      y(i) = alpha*y(i)
    end do
  else 
    allocate(tmp(m), stat=info) 
    if (info /= 0) then 
      info=4010
      call psb_errpush(info,name,a_err='allocate')
      goto 9999
    end if

    tmp(1:m) = x(1:m)
    call inner_coosv(tra,a,tmp,y,info)
    if (info /= 0) then 
      call psb_errpush(info,name)
      goto 9999
    end if
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

  subroutine inner_coosv(tra,a,x,y,info) 
    implicit none 
    logical, intent(in)                 :: tra  
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: x(:)
    real(psb_dpk_), intent(out)         :: y(:)
    integer, intent(out) :: info

    integer :: i,j,k,m, ir, jc, nnz
    real(psb_dpk_) :: acc

    info = 0
    if (.not.a%is_sorted()) then 
      info = 1121
      return
    end if

    nnz = a%get_nzeros()

    if (.not.tra) then 

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          j = 1
          do i=1, a%get_nrows()
            acc = dzero
            do 
              if (j > nnz) exit
              if (a%ia(j) > i) exit
              acc = acc + a%val(j)*x(a%ja(j))
              j   = j + 1
            end do
            y(i) = x(i) - acc
          end do
        else if (.not.a%is_unit()) then 
          j = 1
          do i=1, a%get_nrows()
            acc = dzero
            do 
              if (j > nnz) exit
              if (a%ia(j) > i) exit
              if (a%ja(j) == i) then 
                y(i) = (x(i) - acc)/a%val(j)
                j = j + 1
                exit
              end if
              acc = acc + a%val(j)*x(a%ja(j))
              j   = j + 1
            end do
          end do
        end if

      else if (a%is_upper()) then 
        if (a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              acc = acc + a%val(j)*x(a%ja(j))
              j   = j - 1
            end do
            y(i) = x(i) - acc
          end do

        else if (.not.a%is_unit()) then 

          j = nnz
          do i=a%get_nrows(), 1, -1 
            acc = dzero 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              if (a%ja(j) == i) then 
                y(i) = (x(i) - acc)/a%val(j)
                j = j - 1
                exit
              end if
              acc = acc + a%val(j)*x(a%ja(j))
              j   = j - 1
            end do
          end do
        end if

      end if

    else if (tra) then 

      do i=1, a%get_nrows()
        y(i) = x(i)
      end do

      if (a%is_lower()) then 
        if (a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1
            acc = y(i) 
            do
              if (j < 1) exit
              if (a%ia(j) < i) exit
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
              j     = j - 1 
            end do
          end do
        else if (.not.a%is_unit()) then 
          j = nnz
          do i=a%get_nrows(), 1, -1
            if (a%ja(j) == i) then 
              y(i) = y(i) /a%val(j)
              j    = j - 1
            end if
            acc  = y(i) 
            do 
              if (j < 1) exit
              if (a%ia(j) < i) exit
              jc    = a%ja(j)
              y(jc) = y(jc) - a%val(j)*acc 
              j     = j - 1
            end do
          end do

        else if (a%is_upper()) then 
          if (a%is_unit()) then 
            j = 1
            do i=1, a%get_nrows()
              acc = y(i)
              do 
                if (j > nnz) exit
                if (a%ia(j) > i) exit
                jc    = a%ja(j)
                y(jc) = y(jc) - a%val(j)*acc 
                j   = j + 1
              end do
            end do
          else if (.not.a%is_unit()) then 
            j = 1
            do i=1, a%get_nrows()
              if (a%ja(j) == i) then 
                y(i) = y(i) /a%val(j)
                j    = j + 1
              end if
              acc = y(i)
              do 
                if (j > nnz) exit
                if (a%ia(j) > i) exit
                jc    = a%ja(j)
                y(jc) = y(jc) - a%val(j)*acc 
                j   = j + 1
              end do
            end do
          end if
        end if
      end if
    end if

  end subroutine inner_coosv


end subroutine d_coo_cssv_impl

subroutine d_coo_csmv_impl(alpha,a,x,beta,y,info,trans) 
  use psb_const_mod
  use psb_error_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_csMv_impl
  implicit none 

  class(psbn_d_coo_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_coo_csmv_impl'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if

  tra = ((trans_=='T').or.(trans_=='t'))



  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if
  nnz = a%get_nzeros()

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
  else 
    if (.not.a%is_unit()) then 
      if (beta == dzero) then
        do i = 1, m
          y(i) = dzero
        enddo
      else
        do  i = 1, m
          y(i) = beta*y(i)
        end do
      endif

    else  if (a%is_unit()) then 
      if (beta == dzero) then
        do i = 1, min(m,n)
          y(i) = alpha*x(i)
        enddo
        do i = min(m,n)+1, m
          y(i) = dzero
        enddo
      else
        do  i = 1, min(m,n) 
          y(i) = beta*y(i) + alpha*x(i)
        end do
        do i = min(m,n)+1, m
          y(i) = beta*y(i)
        enddo
      endif

    endif

  end if

  if (.not.tra) then 
    i    = 1
    j    = i
    if (nnz > 0) then 
      ir   = a%ia(1) 
      acc  = dzero
      do 
        if (i>nnz) then 
          y(ir) = y(ir) + alpha * acc
          exit
        endif
        if (a%ia(i) /= ir) then 
          y(ir) = y(ir) + alpha * acc
          ir    = a%ia(i) 
          acc   = dzero
        endif
        acc     = acc + a%val(i) * x(a%ja(i))
        i       = i + 1               
      enddo
    end if

  else if (tra) then 
    if (alpha.eq.done) then
      i    = 1
      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir) = y(ir) +  a%val(i)*x(jc)
      enddo

    else if (alpha.eq.-done) then

      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir) = y(ir) - a%val(i)*x(jc)
      enddo

    else                    

      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir) = y(ir) + alpha*a%val(i)*x(jc)
      enddo

    end if                  !.....end testing on alpha

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine d_coo_csmv_impl


subroutine d_coo_csmm_impl(alpha,a,x,beta,y,info,trans) 
  use psb_const_mod
  use psb_error_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_csmm_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer, intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_coo_csmm_impl'
  logical, parameter :: debug=.false.

  info = 0
  call psb_erractionsave(err_act)


  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif


  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if


  tra = ((trans_=='T').or.(trans_=='t'))

  if (tra) then 
    m = a%get_ncols()
    n = a%get_nrows()
  else
    n = a%get_ncols()
    m = a%get_nrows()
  end if
  nnz = a%get_nzeros()

  nc = min(size(x,2), size(y,2))
  allocate(acc(nc),stat=info)
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
  else 
    if (.not.a%is_unit()) then 
      if (beta == dzero) then
        do i = 1, m
          y(i,:) = dzero
        enddo
      else
        do  i = 1, m
          y(i,:) = beta*y(i,:)
        end do
      endif

    else  if (a%is_unit()) then 
      if (beta == dzero) then
        do i = 1, min(m,n)
          y(i,:) = alpha*x(i,:)
        enddo
        do i = min(m,n)+1, m
          y(i,:) = dzero
        enddo
      else
        do  i = 1, min(m,n) 
          y(i,:) = beta*y(i,:) + alpha*x(i,:)
        end do
        do i = min(m,n)+1, m
          y(i,:) = beta*y(i,:)
        enddo
      endif

    endif

  end if

  if (.not.tra) then 
    i    = 1
    j    = i
    if (nnz > 0) then 
      ir   = a%ia(1) 
      acc  = dzero
      do 
        if (i>nnz) then 
          y(ir,:) = y(ir,:) + alpha * acc
          exit
        endif
        if (a%ia(i) /= ir) then 
          y(ir,:) = y(ir,:) + alpha * acc
          ir    = a%ia(i) 
          acc   = dzero
        endif
        acc     = acc + a%val(i) * x(a%ja(i),:)
        i       = i + 1               
      enddo
    end if

  else if (tra) then 
    if (alpha.eq.done) then
      i    = 1
      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir,:) = y(ir,:) +  a%val(i)*x(jc,:)
      enddo

    else if (alpha.eq.-done) then

      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir,:) = y(ir,:) - a%val(i)*x(jc,:)
      enddo

    else                    

      do i=1,nnz
        ir = a%ja(i)
        jc = a%ia(i)
        y(ir,:) = y(ir,:) + alpha*a%val(i)*x(jc,:)
      enddo

    end if                  !.....end testing on alpha

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine d_coo_csmm_impl

function d_coo_csnmi_impl(a) result(res)
  use psb_error_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_csnmi_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer   :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:,:)
  logical   :: tra
  Integer :: err_act
  character(len=20)  :: name='d_base_csnmi'
  logical, parameter :: debug=.false.


  res = dzero 
  nnz = a%get_nzeros()
  i   = 1
  j   = i
  do while (i<=nnz) 
    do while ((a%ia(j) == a%ia(i)).and. (j <= nnz))
      j = j+1
    enddo
    acc = dzero
    do k=i, j-1
      acc = acc + abs(a%val(k))
    end do
    res = max(res,acc)
  end do

end function d_coo_csnmi_impl



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



subroutine d_coo_csput_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_coo_csput_impl
  implicit none 
    
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)


  Integer            :: err_act
  character(len=20)  :: name='d_coo_csput_impl'
  logical, parameter :: debug=.false.
  integer            :: nza, i,j,k, nzl, isza, int_err(5)

  info = 0
  call psb_erractionsave(err_act)
  
  if (nz <= 0) then 
    info = 10
    int_err(1)=1
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = 35
    int_err(1)=2
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = 35
    int_err(1)=3
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(val) < nz) then 
    info = 35
    int_err(1)=4
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (nz == 0) return


  nza  = a%get_nzeros()
  isza = a%get_size()
!!$  write(0,*) 'On entry to csput_impl: ',nza
  if (a%is_bld()) then 
    ! Build phase. Must handle reallocations in a sensible way.
    if (isza < (nza+nz)) then 
      call a%reallocate(max(nza+nz,int(1.5*isza)))
      isza = a%get_size()
!!$      write(0,*) 'isza: ',isza,nza+nz
    endif

    call psb_inner_ins(nz,ia,ja,val,nza,a%ia,a%ja,a%val,isza,&
         & imin,imax,jmin,jmax,info,gtl)
    call a%set_nzeros(nza)
    

  else  if (a%is_upd()) then 

    call  d_coo_srch_upd(nz,ia,ja,val,a,&
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

  subroutine psb_inner_ins(nz,ia,ja,val,nza,ia1,ia2,aspk,maxsz,&
       & imin,imax,jmin,jmax,info,gtl)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,maxsz
    integer, intent(in) :: ia(:),ja(:)
    integer, intent(inout) :: nza,ia1(:),ia2(:)
    real(psb_dpk_), intent(in) :: val(:)
    real(psb_dpk_), intent(inout) :: aspk(:)
    integer, intent(out) :: info
    integer, intent(in), optional  :: gtl(:)
    integer :: i,ir,ic,ng

    info = 0
    if (present(gtl)) then 
      ng = size(gtl) 

      do i=1, nz 
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
            nza = nza + 1 
            if (nza > maxsz) then 
              info = -91
              return
            endif
            ia1(nza) = ir
            ia2(nza) = ic
            aspk(nza) = val(i)
          end if
        end if
      end do
    else

      do i=1, nz 
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
          nza = nza + 1 
          if (nza > maxsz) then 
            write(0,*) 'err -92 ',nza,maxsz
            info = -92
            return
          endif
          ia1(nza) = ir
          ia2(nza) = ic
          aspk(nza) = val(i)
        end if
      end do
    end if

  end subroutine psb_inner_ins


  subroutine d_coo_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info,gtl)

    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax
    integer, intent(in) :: ia(:),ja(:)
    real(psb_dpk_), intent(in) :: val(:)
    integer, intent(out) :: info
    integer, intent(in), optional  :: gtl(:)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,nnz,dupl,ng
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='d_coo_srch_upd'

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
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            if ((ir > 0).and.(ir <= a%m)) then 
              ic = gtl(ic) 
              if (ir /= ilr) then 
                i1 = psb_ibsrch(ir,nnz,a%ia)
                i2 = i1
                do 
                  if (i2+1 > nnz) exit
                  if (a%ia(i2+1) /= a%ia(i2)) exit
                  i2 = i2 + 1
                end do
                do 
                  if (i1-1 < 1) exit
                  if (a%ia(i1-1) /= a%ia(i1)) exit
                  i1 = i1 - 1
                end do
                ilr = ir
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              ip = psb_issrch(ic,nc,a%ja(i1:i2))
              if (ip>0) then 
                a%val(i1+ip-1) = val(i)
              else
                info = i 
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
            endif
          end if
        end do
      case(psbn_dupl_add_)
        ! Add
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            ic = gtl(ic) 
            if ((ir > 0).and.(ir <= a%m)) then 

              if (ir /= ilr) then 
                i1 = psb_ibsrch(ir,nnz,a%ia)
                i2 = i1
                do 
                  if (i2+1 > nnz) exit
                  if (a%ia(i2+1) /= a%ia(i2)) exit
                  i2 = i2 + 1
                end do
                do 
                  if (i1-1 < 1) exit
                  if (a%ia(i1-1) /= a%ia(i1)) exit
                  i1 = i1 - 1
                end do
                ilr = ir
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              ip = psb_issrch(ic,nc,a%ja(i1:i2))
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
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 

            if (ir /= ilr) then 
              i1 = psb_ibsrch(ir,nnz,a%ia)
              i2 = i1
              do 
                if (i2+1 > nnz) exit
                if (a%ia(i2+1) /= a%ia(i2)) exit
                i2 = i2 + 1
              end do
              do 
                if (i1-1 < 1) exit
                if (a%ia(i1-1) /= a%ia(i1)) exit
                i1 = i1 - 1
              end do
              ilr = ir
            else
              i1 = 1
              i2 = 1
            end if
            nc = i2-i1+1
            ip = psb_issrch(ic,nc,a%ja(i1:i2))
            if (ip>0) then 
              a%val(i1+ip-1) = val(i)
            else
              info = i 
              return
            end if
          end if
        end do

      case(psbn_dupl_add_)
        ! Add
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 

            if (ir /= ilr) then 
              i1 = psb_ibsrch(ir,nnz,a%ia)
              i2 = i1
              do 
                if (i2+1 > nnz) exit
                if (a%ia(i2+1) /= a%ia(i2)) exit
                i2 = i2 + 1
              end do
              do 
                if (i1-1 < 1) exit
                if (a%ia(i1-1) /= a%ia(i1)) exit
                i1 = i1 - 1
              end do
              ilr = ir
            else
              i1 = 1
              i2 = 1
            end if
            nc = i2-i1+1
            ip = psb_issrch(ic,nc,a%ja(i1:i2))
            if (ip>0) then 
              a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
            else
              info = i 
              return
            end if
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine d_coo_srch_upd

end subroutine d_coo_csput_impl


subroutine d_cp_coo_to_coo_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_cp_coo_to_coo_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
  class(psbn_d_coo_sparse_mat), intent(out) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  info = 0
  call b%set_nzeros(a%get_nzeros())
  call b%set_nrows(a%get_nrows())
  call b%set_ncols(a%get_ncols())
  call b%set_dupl(a%get_dupl())
  call b%set_state(a%get_state())
  call b%set_triangle(a%is_triangle())
  call b%set_upper(a%is_upper()) 

  call b%reallocate(a%get_nzeros())

  b%ia(:)  = a%ia(:)
  b%ja(:)  = a%ja(:)
  b%val(:) = a%val(:)

  call b%fix(info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_cp_coo_to_coo_impl
  
subroutine d_cp_coo_from_coo_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_cp_coo_from_coo_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(in) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.
  integer :: m,n,nz


  call psb_erractionsave(err_act)
  info = 0
  call a%set_nzeros(b%get_nzeros())
  call a%set_nrows(b%get_nrows())
  call a%set_ncols(b%get_ncols())
  call a%set_dupl(b%get_dupl())
  call a%set_state(b%get_state())
  call a%set_triangle(b%is_triangle())
  call a%set_upper(b%is_upper())

  call a%reallocate(b%get_nzeros())

  a%ia(:)  = b%ia(:)
  a%ja(:)  = b%ja(:)
  a%val(:) = b%val(:)

  call a%fix(info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_cp_coo_from_coo_impl


subroutine d_cp_coo_to_fmt_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_cp_coo_to_fmt_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(in) :: a
  class(psbn_d_base_sparse_mat), intent(out) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  info = 0

  call b%cp_from_coo(a,info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_cp_coo_to_fmt_impl
  
subroutine d_cp_coo_from_fmt_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_cp_coo_from_fmt_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_base_sparse_mat), intent(in) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.
  integer :: m,n,nz


  call psb_erractionsave(err_act)
  info = 0

  call b%cp_to_coo(a,info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_cp_coo_from_fmt_impl


subroutine d_fix_coo_impl(a,info,idir) 
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_fix_coo_impl
  use psb_string_mod
  use psb_ip_reord_mod
  implicit none 

  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  integer, intent(out)                :: info
  integer, intent(in), optional :: idir
  integer, allocatable :: iaux(:)
  !locals
  Integer              :: nza, nzl,iret,idir_, dupl_
  integer              :: i,j, irw, icl, err_act
  integer              :: debug_level, debug_unit
  character(len=20)    :: name = 'psb_fixcoo'

  info  = 0

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if(debug_level >= psb_debug_serial_) &
       & write(debug_unit,*)  trim(name),': start ',&
       & size(a%ia),size(a%ja)
  if (present(idir)) then 
    idir_ = idir
  else
    idir_ = 0
  endif

  nza = a%get_nzeros()
  if (nza < 2) return

  dupl_ = a%get_dupl()

  call d_fix_coo_inner(nza,dupl_,a%ia,a%ja,a%val,i,info,idir_)

  call a%set_sorted()
  call a%set_nzeros(i)
  call a%set_asb()
  

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine d_fix_coo_impl



subroutine d_fix_coo_inner(nzin,dupl,ia,ja,val,nzout,info,idir) 
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_fix_coo_inner
  use psb_string_mod
  use psb_ip_reord_mod
  implicit none 
  
  integer, intent(in)           :: nzin, dupl
  integer, intent(inout)        :: ia(:), ja(:)
  real(psb_dpk_), intent(inout) :: val(:)
  integer, intent(out)          :: nzout, info
  integer, intent(in), optional :: idir
  !locals
  integer, allocatable :: iaux(:)
  Integer              :: nza, nzl,iret,idir_, dupl_
  integer              :: i,j, irw, icl, err_act
  integer              :: debug_level, debug_unit
  character(len=20)    :: name = 'psb_fixcoo'

  info  = 0

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if(debug_level >= psb_debug_serial_) &
       & write(debug_unit,*)  trim(name),': start ',&
       & size(ia),size(ja)
  if (present(idir)) then 
    idir_ = idir
  else
    idir_ = 0
  endif


  if (nzin < 2) return

  dupl_ = dupl
  
  allocate(iaux(nzin+2),stat=info) 
  if (info /= 0) return


  select case(idir_) 

  case(0) !  Row major order

    call msort_up(nzin,ia(1),iaux(1),iret)
    if (iret == 0) &
         & call psb_ip_reord(nzin,val,ia,ja,iaux)
    i    = 1
    j    = i
    do while (i <= nzin)
      do while ((ia(j) == ia(i)))
        j = j+1
        if (j > nzin) exit
      enddo
      nzl = j - i
      call msort_up(nzl,ja(i),iaux(1),iret)
      if (iret == 0) &
           & call psb_ip_reord(nzl,val(i:i+nzl-1),&
           & ia(i:i+nzl-1),ja(i:i+nzl-1),iaux)
      i = j
    enddo

    i = 1
    irw = ia(i)
    icl = ja(i)
    j = 1

    select case(dupl_)
    case(psbn_dupl_ovwrt_)

      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          val(i) = val(j)
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo

    case(psbn_dupl_add_)

      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          val(i) = val(i) + val(j)
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo

    case(psbn_dupl_err_)
      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          call psb_errpush(130,name)          
          goto 9999
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo
    case default
      write(0,*) 'Error in fix_coo: unsafe dupl',dupl_

    end select


    if(debug_level >= psb_debug_serial_)&
         & write(debug_unit,*)  trim(name),': end second loop'

  case(1) !  Col major order

    call msort_up(nzin,ja(1),iaux(1),iret)
    if (iret == 0) &
         & call psb_ip_reord(nzin,val,ia,ja,iaux)
    i    = 1
    j    = i
    do while (i <= nzin)
      do while ((ja(j) == ja(i)))
        j = j+1
        if (j > nzin) exit
      enddo
      nzl = j - i
      call msort_up(nzl,ia(i),iaux(1),iret)
      if (iret == 0) &
           & call psb_ip_reord(nzl,val(i:i+nzl-1),&
           & ia(i:i+nzl-1),ja(i:i+nzl-1),iaux)
      i = j
    enddo

    i = 1
    irw = ia(i)
    icl = ja(i)
    j = 1


    select case(dupl_)
    case(psbn_dupl_ovwrt_)
      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          val(i) = val(j)
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo

    case(psbn_dupl_add_)
      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          val(i) = val(i) + val(j)
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo

    case(psbn_dupl_err_)
      do 
        j = j + 1
        if (j > nzin) exit
        if ((ia(j) == irw).and.(ja(j) == icl)) then 
          call psb_errpush(130,name)
          goto 9999
        else
          i = i+1
          val(i) = val(j)
          ia(i) = ia(j)
          ja(i) = ja(j)
          irw = ia(i) 
          icl = ja(i) 
        endif
      enddo
    case default
      write(0,*) 'Error in fix_coo: unsafe dupl',dupl_
    end select
    if (debug_level >= psb_debug_serial_)&
         & write(debug_unit,*)  trim(name),': end second loop'
  case default
    write(debug_unit,*) trim(name),': unknown direction ',idir_
  end select

  nzout = i 
  
  deallocate(iaux)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return



end subroutine d_fix_coo_inner




subroutine d_mv_coo_to_coo_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_mv_coo_to_coo_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(out) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  info = 0
  call b%set_nzeros(a%get_nzeros())
  call b%set_nrows(a%get_nrows())
  call b%set_ncols(a%get_ncols())
  call b%set_dupl(a%get_dupl())
  call b%set_state(a%get_state())
  call b%set_triangle(a%is_triangle())
  call b%set_upper(a%is_upper()) 

  call b%reallocate(a%get_nzeros())

  call move_alloc(a%ia, b%ia)
  call move_alloc(a%ja, b%ja)
  call move_alloc(a%val, b%val)
  call a%free()

  call b%fix(info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_mv_coo_to_coo_impl
  
subroutine d_mv_coo_from_coo_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_mv_coo_from_coo_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.
  integer :: m,n,nz


  call psb_erractionsave(err_act)
  info = 0
  call a%set_nzeros(b%get_nzeros())
  call a%set_nrows(b%get_nrows())
  call a%set_ncols(b%get_ncols())
  call a%set_dupl(b%get_dupl())
  call a%set_state(b%get_state())
  call a%set_triangle(b%is_triangle())
  call a%set_upper(b%is_upper())

  call a%reallocate(b%get_nzeros())

  call move_alloc(b%ia , a%ia   )
  call move_alloc(b%ja , a%ja   )
  call move_alloc(b%val, a%val )
  call b%free()

  a%ia(:)  = b%ia(:)
  a%ja(:)  = b%ja(:)
  a%val(:) = b%val(:)

  call a%fix(info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_mv_coo_from_coo_impl


subroutine d_mv_coo_to_fmt_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_mv_coo_to_fmt_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_base_sparse_mat), intent(out) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  info = 0

  call b%mv_from_coo(a,info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_mv_coo_to_fmt_impl
  
subroutine d_mv_coo_from_fmt_impl(a,b,info) 
  use psb_error_mod
  use psb_realloc_mod
  use psbn_d_base_mat_mod, psb_protect_name => d_mv_coo_from_fmt_impl
  implicit none 
  class(psbn_d_coo_sparse_mat), intent(inout) :: a
  class(psbn_d_base_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.
  integer :: m,n,nz


  call psb_erractionsave(err_act)
  info = 0

  call b%mv_to_coo(a,info)

  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  call psb_errpush(info,name)

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine d_mv_coo_from_fmt_impl
