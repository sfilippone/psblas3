module psbn_d_csr_mat_mod

  use psbn_d_base_mat_mod

  type, extends(psbn_d_base_sparse_mat) :: psbn_d_csr_sparse_mat

    integer, allocatable :: irp(:), ja(:)
    real(psb_dpk_), allocatable :: val(:)

  contains
    procedure, pass(a) :: get_nzeros => d_csr_get_nzeros
    procedure, pass(a) :: get_fmt  => d_csr_get_fmt
    procedure, pass(a) :: d_base_csmm => d_csr_csmm
    procedure, pass(a) :: d_base_csmv => d_csr_csmv
    procedure, pass(a) :: d_base_cssm => d_csr_cssm
    procedure, pass(a) :: d_base_cssv => d_csr_cssv
    procedure, pass(a) :: csnmi => d_csr_csnmi
    procedure, pass(a) :: reallocate_nz => d_csr_reallocate_nz
    procedure, pass(a) :: csput => d_csr_csput
    procedure, pass(a) :: allocate_mnnz => d_csr_allocate_mnnz
    procedure, pass(a) :: cp_to_coo => d_cp_csr_to_coo
    procedure, pass(a) :: cp_from_coo => d_cp_csr_from_coo
    procedure, pass(a) :: cp_to_fmt => d_cp_csr_to_fmt
    procedure, pass(a) :: cp_from_fmt => d_cp_csr_from_fmt
    procedure, pass(a) :: mv_to_coo => d_mv_csr_to_coo
    procedure, pass(a) :: mv_from_coo => d_mv_csr_from_coo
    procedure, pass(a) :: mv_to_fmt => d_mv_csr_to_fmt
    procedure, pass(a) :: mv_from_fmt => d_mv_csr_from_fmt
    procedure, pass(a) :: free => d_csr_free
    procedure, pass(a) :: print => d_csr_print
  end type psbn_d_csr_sparse_mat
  private :: d_csr_get_nzeros, d_csr_csmm, d_csr_csmv, d_csr_cssm, d_csr_cssv, &
       & d_csr_csput, d_csr_reallocate_nz, d_csr_allocate_mnnz, &
       & d_csr_free,  d_csr_print, d_csr_get_fmt, d_csr_csnmi, &
       & d_cp_csr_to_coo, d_cp_csr_from_coo, &
       & d_mv_csr_to_coo, d_mv_csr_from_coo, &
       & d_cp_csr_to_fmt, d_cp_csr_from_fmt, &
       & d_mv_csr_to_fmt, d_mv_csr_from_fmt


  interface 
    subroutine d_cp_csr_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      class(psbn_d_base_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine d_cp_csr_to_fmt_impl
  end interface

  interface 
    subroutine d_cp_csr_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      class(psbn_d_base_sparse_mat), intent(in) :: b
      integer, intent(out)                        :: info
    end subroutine d_cp_csr_from_fmt_impl
  end interface


  interface 
    subroutine d_cp_csr_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      class(psbn_d_coo_sparse_mat), intent(out) :: b
      integer, intent(out)            :: info
    end subroutine d_cp_csr_to_coo_impl
  end interface

  interface 
    subroutine d_cp_csr_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(in)    :: b
      integer, intent(out)                        :: info
    end subroutine d_cp_csr_from_coo_impl
  end interface

  interface 
    subroutine d_mv_csr_to_fmt_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      class(psbn_d_base_sparse_mat), intent(out)  :: b
      integer, intent(out)            :: info
    end subroutine d_mv_csr_to_fmt_impl
  end interface

  interface 
    subroutine d_mv_csr_from_fmt_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout)  :: a
      class(psbn_d_base_sparse_mat), intent(inout) :: b
      integer, intent(out)                         :: info
    end subroutine d_mv_csr_from_fmt_impl
  end interface


  interface 
    subroutine d_mv_csr_to_coo_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(out)   :: b
      integer, intent(out)            :: info
    end subroutine d_mv_csr_to_coo_impl
  end interface

  interface 
    subroutine d_mv_csr_from_coo_impl(a,b,info) 
      use psb_const_mod
      use psbn_d_base_mat_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      class(psbn_d_coo_sparse_mat), intent(inout) :: b
      integer, intent(out)                        :: info
    end subroutine d_mv_csr_from_coo_impl
  end interface

  interface 
    subroutine d_csr_csput_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(inout) :: a
      real(psb_dpk_), intent(in)      :: val(:)
      integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
      integer, intent(out)            :: info
      integer, intent(in), optional   :: gtl(:)
    end subroutine d_csr_csput_impl
  end interface

  interface d_csr_cssm_impl
    subroutine d_csr_cssv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_csr_cssv_impl
    subroutine d_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_csr_cssm_impl
  end interface

  interface d_csr_csmm_impl
    subroutine d_csr_csmv_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
      real(psb_dpk_), intent(inout)       :: y(:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_csr_csmv_impl
    subroutine d_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
      real(psb_dpk_), intent(inout)       :: y(:,:)
      integer, intent(out)                :: info
      character, optional, intent(in)     :: trans
    end subroutine d_csr_csmm_impl
  end interface

  interface d_csr_csnmi_impl
    function d_csr_csnmi_impl(a) result(res)
      use psb_const_mod
      import psbn_d_csr_sparse_mat
      class(psbn_d_csr_sparse_mat), intent(in) :: a
      real(psb_dpk_)         :: res
    end function d_csr_csnmi_impl
  end interface
  


contains 

  !=====================================
  !
  !
  !
  ! Getters 
  !
  !
  !
  !
  !
  !=====================================

  function d_csr_get_fmt(a) result(res)
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    character(len=5) :: res
    res = 'CSR'
  end function d_csr_get_fmt
  
  function d_csr_get_nzeros(a) result(res)
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    integer :: res
    res = a%irp(a%m+1)-1
  end function d_csr_get_nzeros


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


  subroutine  d_csr_reallocate_nz(nz,a) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: nz
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    Integer :: err_act, info
    character(len=20)  :: name='d_csr_reallocate_nz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)

    call psb_realloc(nz,a%ja,info)
    if (info == 0) call psb_realloc(nz,a%val,info)
    if (info == 0) call psb_realloc(max(nz,a%m+1,a%n+1),a%irp,info)
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

  subroutine d_csr_csput(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
    use psb_const_mod
    use psb_error_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    real(psb_dpk_), intent(in)      :: val(:)
    integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
    integer, intent(out)            :: info
    integer, intent(in), optional   :: gtl(:)


    Integer            :: err_act
    character(len=20)  :: name='d_csr_csput'
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

    call d_csr_csput_impl(nz,val,ia,ja,a,imin,imax,jmin,jmax,info,gtl) 
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
  end subroutine d_csr_csput


  subroutine  d_csr_free(a) 
    implicit none 

    class(psbn_d_csr_sparse_mat), intent(inout) :: a

    if (allocated(a%irp)) deallocate(a%irp)
    if (allocated(a%ja)) deallocate(a%ja)
    if (allocated(a%val)) deallocate(a%val)
    call a%set_null()
    call a%set_nrows(0)
    call a%set_ncols(0)
    
    return

  end subroutine d_csr_free


  subroutine d_cp_csr_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    class(psbn_d_coo_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_cp_csr_to_coo_impl(a,b,info)
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

  end subroutine d_cp_csr_to_coo
  
  subroutine d_cp_csr_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_cp_csr_from_coo_impl(a,b,info)
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

  end subroutine d_cp_csr_from_coo


  subroutine d_cp_csr_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    class(psbn_d_base_sparse_mat), intent(out) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_cp_csr_to_fmt_impl(a,b,info)
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

  end subroutine d_cp_csr_to_fmt
  
  subroutine d_cp_csr_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(in) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_cp_csr_from_fmt_impl(a,b,info)
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

  end subroutine d_cp_csr_from_fmt


  subroutine d_mv_csr_to_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(out)   :: b
    integer, intent(out)            :: info 

    Integer :: err_act
    character(len=20)  :: name='to_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_mv_csr_to_coo_impl(a,b,info)
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

  end subroutine d_mv_csr_to_coo
  
  subroutine d_mv_csr_from_coo(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    class(psbn_d_coo_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_coo'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_mv_csr_from_coo_impl(a,b,info)
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

  end subroutine d_mv_csr_from_coo


  subroutine d_mv_csr_to_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    class(psbn_d_base_sparse_mat), intent(out)  :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='to_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_mv_csr_to_fmt_impl(a,b,info)
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

  end subroutine d_mv_csr_to_fmt
  
  subroutine d_mv_csr_from_fmt(a,b,info) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(inout)  :: a
    class(psbn_d_base_sparse_mat), intent(inout) :: b
    integer, intent(out)            :: info

    Integer :: err_act
    character(len=20)  :: name='from_fmt'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    call d_mv_csr_from_fmt_impl(a,b,info)
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

  end subroutine d_mv_csr_from_fmt


  subroutine  d_csr_allocate_mnnz(m,n,a,nz) 
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    integer, intent(in) :: m,n
    class(psbn_d_csr_sparse_mat), intent(inout) :: a
    integer, intent(in), optional :: nz
    Integer :: err_act, info, nz_
    character(len=20)  :: name='allocate_mnz'
    logical, parameter :: debug=.false.

    call psb_erractionsave(err_act)
    info = 0
    if (m < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/1,0,0,0,0/))
      goto 9999
    endif
    if (n < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/2,0,0,0,0/))
      goto 9999
    endif
        if (present(nz)) then 
      nz_ = nz
    else
      nz_ = max(7*m,7*n,1)
    end if
    if (nz_ < 0) then 
      info = 10
      call psb_errpush(info,name,i_err=(/3,0,0,0,0/))
      goto 9999
    endif
      
    if (info == 0) call psb_realloc(m+1,a%irp,info)
    if (info == 0) call psb_realloc(nz_,a%ja,info)
    if (info == 0) call psb_realloc(nz_,a%val,info)
    if (info == 0) then 
      a%irp=0
      call a%set_nrows(m)
      call a%set_ncols(n)
      call a%set_bld()
      call a%set_triangle(.false.)
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

  end subroutine d_csr_allocate_mnnz


  subroutine d_csr_print(iout,a,iv,eirs,eics,head,ivr,ivc)
    use psb_spmat_type
    use psb_string_mod
    implicit none 

    integer, intent(in)               :: iout
    class(psbn_d_csr_sparse_mat), intent(in) :: a   
    integer, intent(in), optional     :: iv(:)
    integer, intent(in), optional     :: eirs,eics
    character(len=*), optional        :: head
    integer, intent(in), optional     :: ivr(:), ivc(:)

    Integer :: err_act
    character(len=20)  :: name='d_csr_print'
    logical, parameter :: debug=.false.

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
    
    write(frmtv,'(a,i3.3,a,i3.3,a)') '(2(i',ni,',1x),es26.18,1x,2(i',ni,',1x))'
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

  end subroutine d_csr_print


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


  subroutine d_csr_csmv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
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

    if (.not.a%is_asb()) then 
      info = 1121
      call psb_errpush(info,name)
      goto 9999
    endif


    call d_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine d_csr_csmv

  subroutine d_csr_csmm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
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



    call d_csr_csmm_impl(alpha,a,x,beta,y,info,trans) 
    
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

  end subroutine d_csr_csmm


  subroutine d_csr_cssv(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
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
    
    call d_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 

    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return


  end subroutine d_csr_cssv



  subroutine d_csr_cssm(alpha,a,x,beta,y,info,trans) 
    use psb_error_mod
    implicit none 
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
    character(len=20)  :: name='d_csr_csmm'
    logical, parameter :: debug=.false.

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
    
    call d_csr_cssm_impl(alpha,a,x,beta,y,info,trans) 
    call psb_erractionrestore(err_act)
    return


9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_csr_cssm
 
  function d_csr_csnmi(a) result(res)
    use psb_error_mod
    use psb_const_mod
    implicit none 
    class(psbn_d_csr_sparse_mat), intent(in) :: a
    real(psb_dpk_)         :: res
    
    Integer :: err_act
    character(len=20)  :: name='csnmi'
    logical, parameter :: debug=.false.
    
    
    res = d_csr_csnmi_impl(a)
    
    return

  end function d_csr_csnmi



end module psbn_d_csr_mat_mod
