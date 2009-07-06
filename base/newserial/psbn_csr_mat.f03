module psbn_d_csr_sparse_mat_mod

  use psbn_d_base_mat_mod

  type, extends(psbn_d_base_sparse_mat) :: psbn_d_csr_sparse_mat
    logical              :: sorted
    integer, allocatable :: irp(:), ja(:)
    real(kind(1.d0)), allocatable :: val(:)
  contains
    procedure, pass(a)  :: d_csr_get_nzeros
    procedure, pass(a)  :: d_base_csmm => d_csr_csmm
    procedure, pass(a)  :: d_base_csmv => d_csr_csmv
    generic, public     :: base_get_nzeros => d_csr_get_nzeros
    procedure, pass(a)  :: d_base_cssm => d_csr_cssm
    procedure, pass(a)  :: d_base_cssv => d_csr_cssv


  end type psbn_d_csr_sparse_mat

contains 

  function d_csr_get_nzeros(a) result(res)
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    integer :: res
    res = a%irp(a%m+1)-1
  end function d_csr_get_nzeros
  

  subroutine d_csr_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans

    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_dpk_) :: acc
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_csr_csmv'
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
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = -acc
          end do

        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = alpha*acc
          end do

        end if


      else if (beta == done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = y(i) -acc
          end do

        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = y(i) + alpha*acc
          end do

        end if

      else if (beta == -done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = -y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = -y(i) -acc
          end do

        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = -y(i) + alpha*acc
          end do

        end if

      else 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = beta*y(i) + acc
          end do

        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j))          
            enddo
            y(i) = beta*y(i) - acc
          end do

        else 

          do i=1,m 
            acc  = zero
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

  end subroutine d_csr_csmv

  subroutine d_csr_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_dpk_), allocatable  :: acc(:)
    logical   :: tra
    Integer :: err_act
    character(len=20)  :: name='d_csr_csmm'
    logical, parameter :: debug=.false.

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
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = acc
          end do
        
        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = -acc
          end do
      
        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = alpha*acc
          end do

        end if
        
        
      else if (beta == done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = y(i,:) + acc
          end do
        
        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = y(i,:) -acc
          end do
      
        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = y(i,:) + alpha*acc
          end do

        end if
        
      else if (beta == -done) then 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = -y(i,:) + acc
          end do
        
        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = -y(i,:) -acc
          end do
      
        else 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = -y(i,:) + alpha*acc
          end do

        end if        

      else 

        if (alpha == done) then 
          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = beta*y(i,:) + acc
          end do
        
        else if (alpha == -done) then 

          do i=1,m 
            acc  = zero
            do j=a%irp(i), a%irp(i+1)-1
              acc  = acc + a%val(j) * x(a%ja(j),:)          
            enddo
            y(i,:) = beta*y(i,:) - acc
          end do
      
        else 

          do i=1,m 
            acc  = zero
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

  end subroutine d_csr_csmm


  subroutine d_csr_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
    class(psbn_d_csr_sparse_mat), intent(in) :: a
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
    character(len=20)  :: name='d_csr_cssv'
    logical, parameter :: debug=.false.

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
      call inner_csrsv(tra,a,x,y)
      do  i = 1, m
        y(i) = alpha*y(i)
      end do
    else 
      allocate(tmp(m), stat=info) 
      if (info /= 0) then 
        write(0,*) 'Memory allocation error in CSRSV '
        return
      end if
      tmp(1:m) = x(1:m)
      call inner_csrsv(tra,a,tmp,y)
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

    subroutine inner_csrsv(tra,a,x,y) 
      use psb_const_mod
      logical, intent(in)                 :: tra  
      class(psbn_d_csr_sparse_mat), intent(in) :: a
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
              y(i) = y(i)/a%val(a%irp(i+1)-1)
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
    end subroutine inner_csrsv

  end subroutine d_csr_cssv



  subroutine d_csr_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
    class(psbn_d_csr_sparse_mat), intent(in) :: a
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
      call inner_csrsm(tra,a,x,y,info)
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
      call inner_csrsm(tra,a,tmp,y,info)
      do  i = 1, m
        y(i,:) = alpha*tmp(i,:) + beta*y(i,:)
      end do
    end if

    if(info /= 0) then
      info=4010
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

    subroutine inner_csrsm(tra,a,x,y,info) 
      use psb_const_mod
      logical, intent(in)                 :: tra  
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: x(:,:)
      real(psb_dpk_), intent(out)         :: y(:,:)
      integer, intent(out)                :: info
      integer :: i,j,k,m, ir, jc
      real(psb_dpk_), allocatable  :: acc(:)
      
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
              y(i,:) = y(i,:)/a%val(a%irp(i+1)-1)
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
    end subroutine inner_csrsm

  end subroutine d_csr_cssm

end module psbn_d_csr_sparse_mat_mod

