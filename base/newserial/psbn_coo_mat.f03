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
    procedure, pass(a)  :: d_base_cssm => d_coo_cssm
    procedure, pass(a)  :: d_base_cssv => d_coo_cssv
    procedure, pass(a)  :: d_base_csins => d_coo_csins
    procedure, pass(a)  :: base_reallocate_nz => d_coo_reallocate_nz

  end type psbn_d_coo_sparse_mat

contains 


  subroutine  d_coo_reallocate_nz(nz,a) 
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    integer, intent(in) :: nz
    class(psbn_d_coo_sparse_mat), intent(inout) :: a
    Integer :: err_act
    character(len=20)  :: name='d_coo_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    call psb_realloc(nx,a%ia,a%ja,a%val,info)
    
    if (info /= 0) then 
      call psb_errpush(4000,name)
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

  end subroutine d_coo_reallocate_nz


  function d_coo_get_nzeros(a) result(res)
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    integer :: res
    res = a%nnz
  end function d_coo_get_nzeros
  

  subroutine d_coo_csins(nz,val,ia,ja,a,info) 
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    class(psbn_d_coo_sparse_mat), intent(in) :: a
    real(psb_dpk_), intent(in)    :: val(:)
    integer, intent(in)           :: nz, ia(:), ja(:)
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='d_coo_csins'
    logical, parameter :: debug=.false.
    integer :: nza, i,j,k, nzl, isza, int_err(5)

    call psb_erractionsave(err_act)
    info = 0
    
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

    if (a%is_bld()) then 
      ! Build phase. Must handle reallocations in a sensible way.
      nza  = a%get_nzeros()
      isza = a%get_size()
      if (isza < (nza+nz)) then 
        
        
      endif
      
    else  if (a%is_upd()) then 

    else 
      ! State is wrong.
      info = 1121
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

    subroutine psb_inner_upd(nz,ia,ja,val,nza,aspk,maxsz,&
         & imin,imax,jmin,jmax,nzl,info,gtl,ng)
      implicit none 

      integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,maxsz
      integer, intent(in) :: ia(:),ja(:)
      integer, intent(inout) :: nza
      real(psb_dpk_), intent(in) :: val(:)
      real(psb_dpk_), intent(inout) :: aspk(:)
      integer, intent(out) :: info
      integer, intent(in), optional  :: ng,gtl(:)
      integer  :: i,ir,ic
      character(len=20)    :: name, ch_err


      name='psb_inner_upd'

      if (present(gtl)) then 
        if (.not.present(ng)) then 
          info = -1
          return
        endif
        if ((nza > nzl)) then 
          do i=1, nz 
            nza = nza + 1 
            if (nza>maxsz) then 
              call psb_errpush(50,name,i_err=(/7,maxsz,5,0,nza /))
              info = -71
              return
            endif
            aspk(nza) = val(i)
          end do
        else
          do i=1, nz 
            ir = ia(i)
            ic = ja(i) 
            if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
              ir = gtl(ir)
              ic = gtl(ic) 
              if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
                nza = nza + 1 
                if (nza>maxsz) then 
                  info = -72
                  return
                endif
                aspk(nza) = val(i)
              end if
            end if
          end do
        end if
      else
        if ((nza >= nzl)) then 
          do i=1, nz 
            nza = nza + 1 
            if (nza>maxsz) then 
              info = -73
              return
            endif
            aspk(nza) = val(i)
          end do
        else
          do i=1, nz 
            ir = ia(i)
            ic = ja(i) 
            if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
              nza = nza + 1 
              if (nza>maxsz) then 
                info = -74
                return
              endif
              aspk(nza) = val(i)
            end if
          end do
        end if
      end if
    end subroutine psb_inner_upd

    subroutine psb_inner_ins(nz,ia,ja,val,nza,ia1,ia2,aspk,maxsz,&
         & imin,imax,jmin,jmax,info,gtl,ng)
      implicit none 

      integer, intent(in) :: nz, imin,imax,jmin,jmax,maxsz
      integer, intent(in) :: ia(:),ja(:)
      integer, intent(inout) :: nza,ia1(:),ia2(:)
      real(psb_dpk_), intent(in) :: val(:)
      real(psb_dpk_), intent(inout) :: aspk(:)
      integer, intent(out) :: info
      integer, intent(in), optional  :: ng,gtl(:)
      integer :: i,ir,ic

      info = 0
      if (present(gtl)) then 
        if (.not.present(ng)) then 
          info = -1
          return
        endif
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



  end subroutine d_coo_csins


  subroutine d_coo_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
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
    character(len=20)  :: name='d_co_csmv'
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

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_coo_csmv

  subroutine d_coo_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
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
    character(len=20)  :: name='d_coo_csmm'
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

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_coo_csmm

  subroutine d_coo_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
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
    character(len=20)  :: name='d_coo_cssv'
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
      use psb_const_mod
      logical, intent(in)                 :: tra  
      class(psbn_d_coo_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: x(:)
      real(psb_dpk_), intent(out)         :: y(:)

      integer :: i,j,k,m, ir, jc, nnz
      real(psb_dpk_) :: acc

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
                if (ia(j) > i) exit
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
                if (ia(j) > i) exit
                if (ja(j) == i) then 
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
                if (ia(j) < i) exit
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
                if (ia(j) < i) exit
                if (ja(j) == i) then 
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
                if (ia(j) < i) exit
                jc    = a%ja(j)
                y(jc) = y(jc) - a%val(j)*acc 
                j     = j - 1 
              end do
            end do
          else if (.not.a%is_unit()) then 
            j = nnz
            do i=a%get_nrows(), 1, -1
              if (ja(j) == i) then 
                y(i) = y(i) /a%val(j)
                j    = j - 1
              end if
              acc  = y(i) 
              do 
                if (j < 1) exit
                if (ia(j) < i) exit
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
                  if (ia(j) > i) exit
                  jc    = a%ja(j)
                  y(jc) = y(jc) - a%val(j)*acc 
                  j   = j + 1
                end do
              end do
            else if (.not.a%is_unit()) then 
              j = 1
              do i=1, a%get_nrows()
                if (ja(j) == i) then 
                  y(i) = y(i) /a%val(j)
                  j    = j + 1
                end if
                acc = y(i)
                do 
                  if (j > nnz) exit
                  if (ia(j) > i) exit
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

  end subroutine d_coo_cssv



  subroutine d_coo_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_const_mod
    use psb_error_mod
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
      use psb_const_mod
      logical, intent(in)                 :: tra  
      class(psbn_d_coo_sparse_mat), intent(in) :: a
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
                if (ia(j) > i) exit
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
                if (ia(j) > i) exit
                if (ja(j) == i) then 
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
                if (ia(j) < i) exit
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
                if (ia(j) < i) exit
                if (ja(j) == i) then 
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
                if (ia(j) < i) exit
                jc    = a%ja(j)
                y(jc,:) = y(jc,:) - a%val(j)*acc 
                j     = j - 1 
              end do
            end do
          else if (.not.a%is_unit()) then 
            j = nnz
            do i=a%get_nrows(), 1, -1
              if (ja(j) == i) then 
                y(i,:) = y(i,:) /a%val(j)
                j    = j - 1
              end if
              acc  = y(i,:) 
              do 
                if (j < 1) exit
                if (ia(j) < i) exit
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
                  if (ia(j) > i) exit
                  jc    = a%ja(j)
                  y(jc,:) = y(jc,:) - a%val(j)*acc 
                  j   = j + 1
                end do
              end do
            else if (.not.a%is_unit()) then 
              j = 1
              do i=1, a%get_nrows()
                if (ja(j) == i) then 
                  y(i,:) = y(i,:) /a%val(j)
                  j    = j + 1
                end if
                acc = y(i,:)
                do 
                  if (j > nnz) exit
                  if (ia(j) > i) exit
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

  end subroutine d_coo_cssm


  
end module psbn_d_coo_sparse_mat_mod

