module psbn_d_coo_sparse_mat_mod
  
  use psbn_d_base_mat_mod

  type, extends(psbn_d_base_sparse_mat) :: psbn_d_coo_sparse_mat
    integer              :: nnz,  state
    logical              :: sorted
    
    integer, allocatable :: ia(:), ja(:)
    real(kind(1.d0)), allocatable :: val(:)
  contains
    procedure, pass(a)  :: d_coo_get_nzeros
    procedure, pass(a)  :: d_base_csmm => d_coo_csmm
    procedure, pass(a)  :: d_base_csmv => d_coo_csmv
    generic, public     :: base_get_nzeros => d_coo_get_nzeros


  end type psbn_d_coo_sparse_mat

contains 

  function d_coo_get_nzeros(a) result(res)
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    integer :: res
    res = a%nnz
  end function d_coo_get_nzeros
  

  subroutine d_coo_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout)       :: y(:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc
    real(psb_dpk_) :: acc
    logical   :: tra

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
        acc  = zero
        do 
          if (i>nnz) then 
            y(ir) = y(ir) + alpha * acc
            exit
          endif
          if (ia(i) /= ir) then 
            y(ir) = y(ir) + alpha * acc
            ir    = ia(i) 
            acc   = zero
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
          ir = ja(i)
          jc = ia(i)
          y(ir) = y(ir) + alpha*a%val(i)*x(jc)
        enddo
        
      end if                  !.....end testing on alpha
      
    endif

  end subroutine d_coo_csmv

  subroutine d_coo_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
    real(psb_dpk_), intent(inout)       :: y(:,:)
    integer, intent(out)                :: info
    character, optional, intent(in)     :: trans
    
    character :: trans_
    integer   :: i,j,k,m,n, nnz, ir, jc, nc
    real(psb_dpk_), allocatable  :: acc(:)
    logical   :: tra

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

    nc = size(x,2) 
    if (nc /= size(y,2)) then 
      write(0,*) 'Mismatch in column sizes!!'
      return
    end if
    allocate(acc(nc))

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
        acc  = zero
        do 
          if (i>nnz) then 
            y(ir,:) = y(ir,:) + alpha * acc
            exit
          endif
          if (ia(i) /= ir) then 
            y(ir,:) = y(ir,:) + alpha * acc
            ir    = ia(i) 
            acc   = zero
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
          ir = ja(i)
          jc = ia(i)
          y(ir,:) = y(ir,:) + alpha*a%val(i)*x(jc,:)
        enddo
        
      end if                  !.....end testing on alpha
      
    endif

  end subroutine d_coo_csmm
  
end module psbn_d_coo_sparse_mat_mod

