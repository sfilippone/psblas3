!
program d_matgen
  use psb_base_mod
  use psb_util_mod
  use psb_d_cxx_mat_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim

  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec 

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a
!!$  type(psb_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  type(psb_d_vect_type) :: b, x
  ! blacs parameters
  integer(psb_ipk_) :: ictxt, iam, np

  ! solver parameters
  integer(psb_ipk_) :: iter, itmax,itrace, istopc, irst
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps
  type(psb_d_csr_sparse_mat) :: acsr
  type(psb_d_cxx_sparse_mat) :: acxx
  
  ! other variables
  integer(psb_ipk_) :: info, err_act
  character(len=20)  :: name,ch_err

  info=psb_success_

  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999

  call psb_set_errverbosity(itwo)
  
  !
  !  get parameters
  !
  call get_parms(ictxt,idim)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  if (.true.) then 
    call psb_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
         & a1,a2,a3,b1,b2,b3,c,g,info,amold=acsr)  
  else if (.false.) then 
    call psb_gen_pde3d(ictxt,idim,a,b,x,desc_a,afmt,&
         & a1,a2,a3,b1,b2,b3,c,g,info,amold=acxx)  
  end if
    
    
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1

  call psb_exit(ictxt)
  stop

9999 call psb_error(ictxt)

  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,idim)
    integer(psb_ipk_) :: ictxt
    integer(psb_ipk_) :: idim
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: intbuf(10), ip

    call psb_info(ictxt, iam, np)

    read(psb_inp_unit,*) idim


    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer(psb_ipk_) :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  pde90 methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage
  
    !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=1.d0/sqrt(3.d0)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=1.d0/sqrt(3.d0)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=1.d0/sqrt(3.d0)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=0.d0
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=1.d0/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=1.d0/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=1.d0/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g

end program d_matgen
