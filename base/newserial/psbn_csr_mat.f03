module psbn_d_csr_sparse_mat_mod

  use psbn_d_base_mat_mod

  type, extends(psbn_d_base_sparse_mat) :: psbn_d_csr_sparse_mat

    logical              :: sorted
    integer, allocatable :: irp(:), ja(:)
    real(psb_dpk_), allocatable :: val(:)

  contains
    procedure, pass(a)  :: get_nzeros => d_csr_get_nzeros
    procedure, pass(a)  :: d_base_csmm => d_csr_csmm
    procedure, pass(a)  :: d_base_csmv => d_csr_csmv
    procedure, pass(a)  :: d_base_cssm => d_csr_cssm
    procedure, pass(a)  :: d_base_cssv => d_csr_cssv
    procedure, pass(a)  :: reallocate_nz => d_csr_reallocate_nz
    procedure, pass(a)  :: csins => d_csr_csins


  end type psbn_d_csr_sparse_mat

contains 

  subroutine  d_csr_reallocate_nz(nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    integer, intent(in) :: nz
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='d_csr_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    call psb_realloc(nx,a%ja,info)
    if (info == 0) call psb_realloc(nx,a%val,info)
    if (info == 0) call psb_realloc(max(nx,a%m+1,a%n+1),a%irp,info)
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

  end subroutine d_csr_reallocate_nz

  function d_csr_get_nzeros(a) result(res)
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    integer :: res
    res = a%irp(a%m+1)-1
  end function d_csr_get_nzeros
  

  subroutine d_csr_csins(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_const_mod
    use psb_error_mod
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)


    Integer            :: err_act
    character(len=20)  :: name='d_csr_csins'
    logical, parameter :: debug=.false.
    integer            :: nza, i,j,k, nzl, isza, int_err(5)

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


    nza  = a%get_nzeros()

    if (a%is_bld()) then 
      ! Build phase should only ever be in COO
      info = 1121

    else  if (a%is_upd()) then 
      if (a%is_sorted()) then 


!!$#ifdef FIXED_NAG_SEGV
!!$        call  d_csr_srch_upd(nz,ia,ja,val,a,&
!!$             & imin,imax,jmin,jmax,info,gtl)
!!$#else 
        call  d_csr_srch_upd(nz,ia,ja,val,&
             & a%irp,a%ja,a%val,&
             & a%get_dupl(),a%get_nzeros(),a%get_nrows(),&
             & info,gtl)
!!$#endif

      else
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


!!$#ifdef FIXED_NAG_SEGV
!!$    subroutine d_csr_srch_upd(nz,ia,ja,val,a,&
!!$         & imin,imax,jmin,jmax,info,gtl)
!!$      
!!$      use psb_const_mod
!!$      use psb_realloc_mod
!!$      use psb_string_mod
!!$      use psb_serial_mod
!!$      implicit none 
!!$      
!!$      class(psbn_d_csr_sparse_mat), intent(inout) :: a
!!$      integer, intent(in) :: nz, imin,imax,jmin,jmax
!!$      integer, intent(in) :: ia(:),ja(:)
!!$      real(psb_dpk_), intent(in) :: val(:)
!!$      integer, intent(out) :: info
!!$      integer, intent(in), optional  :: gtl(:)
!!$      integer  :: i,ir,ic, ilr, ilc, ip, &
!!$           & i1,i2,nc,nnz,dupl,ng
!!$      integer              :: debug_level, debug_unit
!!$      character(len=20)    :: name='d_csr_srch_upd'
!!$      
!!$      info = 0
!!$      debug_unit  = psb_get_debug_unit()
!!$      debug_level = psb_get_debug_level()
!!$
!!$      dupl = a%get_dupl()
!!$      
!!$      if (.not.a%is_sorted()) then 
!!$        info = -4
!!$        return
!!$      end if
!!$      
!!$      ilr = -1 
!!$      ilc = -1 
!!$      nnz = a%get_nzeros()
!!$      
!!$      
!!$      if (present(gtl)) then
!!$        ng = size(gtl)
!!$        
!!$        select case(dupl)
!!$        case(psbn_dupl_ovwrt_,psbn_dupl_err_)
!!$          ! Overwrite.
!!$          ! Cannot test for error, should have been caught earlier.
!!$          do i=1, nz
!!$            ir = ia(i)
!!$            ic = ja(i) 
!!$            if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$              ir = gtl(ir)
!!$              if ((ir > 0).and.(ir <= a%m)) then 
!!$                ic = gtl(ic) 
!!$                if (ir /= ilr) then 
!!$                  i1 = psb_ibsrch(ir,nnz,a%ia)
!!$                  i2 = i1
!!$                  do 
!!$                    if (i2+1 > nnz) exit
!!$                    if (a%ia(i2+1) /= a%ia(i2)) exit
!!$                    i2 = i2 + 1
!!$                  end do
!!$                  do 
!!$                    if (i1-1 < 1) exit
!!$                    if (a%ia(i1-1) /= a%ia(i1)) exit
!!$                    i1 = i1 - 1
!!$                  end do
!!$                  ilr = ir
!!$                else
!!$                  i1 = 1
!!$                  i2 = 1
!!$                end if
!!$                nc = i2-i1+1
!!$                ip = psb_issrch(ic,nc,a%ja(i1:i2))
!!$                if (ip>0) then 
!!$                  a%val(i1+ip-1) = val(i)
!!$                else
!!$                  info = i 
!!$                  return
!!$                end if
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Discarding row that does not belong to us.'
!!$              endif
!!$            end if
!!$          end do
!!$        case(psbn_dupl_add_)
!!$          ! Add
!!$          do i=1, nz
!!$            ir = ia(i)
!!$            ic = ja(i) 
!!$            if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
!!$              ir = gtl(ir)
!!$              ic = gtl(ic) 
!!$              if ((ir > 0).and.(ir <= a%m)) then 
!!$                
!!$                if (ir /= ilr) then 
!!$                  i1 = psb_ibsrch(ir,nnz,a%ia)
!!$                  i2 = i1
!!$                  do 
!!$                    if (i2+1 > nnz) exit
!!$                    if (a%ia(i2+1) /= a%ia(i2)) exit
!!$                    i2 = i2 + 1
!!$                  end do
!!$                  do 
!!$                    if (i1-1 < 1) exit
!!$                    if (a%ia(i1-1) /= a%ia(i1)) exit
!!$                    i1 = i1 - 1
!!$                  end do
!!$                  ilr = ir
!!$                else
!!$                  i1 = 1
!!$                  i2 = 1
!!$                end if
!!$                nc = i2-i1+1
!!$                ip = psb_issrch(ic,nc,a%ja(i1:i2))
!!$                if (ip>0) then 
!!$                  a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$                else
!!$                  info = i 
!!$                  return
!!$                end if
!!$              else
!!$                if (debug_level >= psb_debug_serial_) &
!!$                     & write(debug_unit,*) trim(name),&
!!$                     & ': Discarding row that does not belong to us.'              
!!$              end if
!!$            end if
!!$          end do
!!$          
!!$        case default
!!$          info = -3
!!$          if (debug_level >= psb_debug_serial_) &
!!$               & write(debug_unit,*) trim(name),&
!!$               & ': Duplicate handling: ',dupl
!!$        end select
!!$        
!!$      else
!!$        
!!$        select case(dupl)
!!$        case(psbn_dupl_ovwrt_,psbn_dupl_err_)
!!$          ! Overwrite.
!!$          ! Cannot test for error, should have been caught earlier.
!!$          do i=1, nz
!!$            ir = ia(i)
!!$            ic = ja(i) 
!!$            if ((ir > 0).and.(ir <= a%m)) then 
!!$
!!$              if (ir /= ilr) then 
!!$                i1 = psb_ibsrch(ir,nnz,a%ia)
!!$                i2 = i1
!!$                do 
!!$                  if (i2+1 > nnz) exit
!!$                  if (a%ia(i2+1) /= a%ia(i2)) exit
!!$                  i2 = i2 + 1
!!$                end do
!!$                do 
!!$                  if (i1-1 < 1) exit
!!$                  if (a%ia(i1-1) /= a%ia(i1)) exit
!!$                  i1 = i1 - 1
!!$                end do
!!$                ilr = ir
!!$              else
!!$                i1 = 1
!!$                i2 = 1
!!$              end if
!!$              nc = i2-i1+1
!!$              ip = psb_issrch(ic,nc,a%ja(i1:i2))
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = val(i)
!!$              else
!!$                info = i 
!!$                return
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        case(psbn_dupl_add_)
!!$          ! Add
!!$          do i=1, nz
!!$            ir = ia(i)
!!$            ic = ja(i) 
!!$            if ((ir > 0).and.(ir <= a%m)) then 
!!$
!!$              if (ir /= ilr) then 
!!$                i1 = psb_ibsrch(ir,nnz,a%ia)
!!$                i2 = i1
!!$                do 
!!$                  if (i2+1 > nnz) exit
!!$                  if (a%ia(i2+1) /= a%ia(i2)) exit
!!$                  i2 = i2 + 1
!!$                end do
!!$                do 
!!$                  if (i1-1 < 1) exit
!!$                  if (a%ia(i1-1) /= a%ia(i1)) exit
!!$                  i1 = i1 - 1
!!$                end do
!!$                ilr = ir
!!$              else
!!$                i1 = 1
!!$                i2 = 1
!!$              end if
!!$              nc = i2-i1+1
!!$              ip = psb_issrch(ic,nc,a%ja(i1:i2))
!!$              if (ip>0) then 
!!$                a%val(i1+ip-1) = a%val(i1+ip-1) + val(i)
!!$              else
!!$                info = i 
!!$                return
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        case default
!!$          info = -3
!!$          if (debug_level >= psb_debug_serial_) &
!!$               & write(debug_unit,*) trim(name),&
!!$               & ': Duplicate handling: ',dupl
!!$        end select
!!$
!!$      end if
!!$
!!$    end subroutine d_csr_srch_upd
!!$
!!$#else
    subroutine d_csr_srch_upd(nz,ia,ja,val,&
         & airp,aja,aval,dupl,nza,nra,&
         & info,gtl)

      use psb_error_mod
      use psb_sort_mod
      implicit none 

      integer, intent(inout) :: airp(:),aja(:)
      real(psb_dpk_), intent(inout) :: aval(:)
      integer, intent(in) :: nz, dupl,nza, nra
      integer, intent(in) :: ia(:),ja(:)
      real(psb_dpk_), intent(in) :: val(:)
      integer, intent(out) :: info
      integer, intent(in), optional  :: gtl(:)

      integer  :: debug_level, debug_unit
      character(len=20)    :: name='d_csr_srch_upd'
      integer  :: i,ir,ic, ilr, ilc, ip, &
           & i1,i2,nc,lb,ub,m, ng

      info = 0
      debug_unit  = psb_get_debug_unit()
      debug_level = psb_get_debug_level()


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
                i1 = airp(ir)
                i2 = airp(ir+1)
                nc=i2-i1

                ip = psb_ibsrch(ic,nc,aja(i1:i2-1))    
                if (ip>0) then 
                  aval(i1+ip-1) = val(i)
                else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',aja(i1:i2-1)
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
                i1 = airp(ir)
                i2 = airp(ir+1)
                nc = i2-i1
                ip = psb_ibsrch(ic,nc,aja(i1:i2-1))
                if (ip>0) then 
                  aval(i1+ip-1) = aval(i1+ip-1) + val(i)
                else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',aja(i1:i2-1)
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

              i1 = airp(ir)
              i2 = airp(ir+1)
              nc=i2-i1

              ip = psb_ibsrch(ic,nc,aja(i1:i2-1))    
              if (ip>0) then 
                aval(i1+ip-1) = val(i)
              else
                if (debug_level >= psb_debug_serial_) &
                     & write(debug_unit,*) trim(name),&
                     & ': Was searching ',ic,' in: ',i1,i2,&
                     & ' : ',aja(i1:i2-1)
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
              i1 = airp(ir)
              i2 = airp(ir+1)
              nc = i2-i1
              ip = psb_ibsrch(ic,nc,aja(i1:i2-1))
              if (ip>0) then 
                aval(i1+ip-1) = aval(i1+ip-1) + val(i)
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
    end subroutine d_csr_srch_upd
!!$#endif
  end subroutine d_csr_csins






  subroutine d_csr_csmv(alpha,a,x,beta,y,info,trans) 
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

