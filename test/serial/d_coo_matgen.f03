!
program d_coo_matgen
  use psb_sparse_mod
!!$  use psb_prec_mod
!!$  use psb_krylov_mod
  use psb_d_base_mat_mod
  use psb_d_csr_mat_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer   :: idim

  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec 

  ! sparse matrix and preconditioner
  type(psb_d_sparse_mat) :: a
!!$  type(psb_dprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  real(psb_dpk_), allocatable :: b(:), x(:)
  ! blacs parameters
  integer            :: ictxt, iam, np

  ! solver parameters
  integer            :: iter, itmax,itrace, istopc, irst
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_dpk_)   :: err, eps

  ! other variables
  integer            :: info, err_act
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

  call psb_set_errverbosity(2)

  !
  !  get parameters
  !
  call get_parms(ictxt,idim)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call create_matrix(idim,a,b,x,desc_a,ictxt,afmt,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    call psb_error(ictxt)
  end if

  call psb_exit(ictxt)
  stop

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
  end if
  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,idim)
    integer      :: ictxt
    integer      :: idim
    integer      :: np, iam
    integer      :: intbuf(10), ip

    call psb_info(ictxt, iam, np)

    read(psb_inp_unit,*) idim


    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer :: iout
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
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs. 
  !
  subroutine create_matrix(idim,a,b,xv,desc_a,ictxt,afmt,info)
    !
    !   discretize the partial diferential equation
    ! 
    !   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
    ! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u 
    !      dxdx     dydy       dzdz        dx       dy         dz   
    !
    ! with Dirichlet boundary conditions, on the unit cube  0<=x,y,z<=1.
    !
    ! Boundary conditions are set in a very simple way, by adding 
    ! equations of the form
    !
    !   u(x,y) = exp(-x^2-y^2-z^2)
    !
    ! Note that if a1=a2=a3=a4=0., the PDE is the well-known Laplace equation.
    !
    use psb_sparse_mod
!!$    use psb_d_base_mat_mod
    use psb_d_csr_mat_mod  
    implicit none
    integer                        :: idim
    integer, parameter             :: nb=20
    real(psb_dpk_), allocatable    :: b(:),xv(:)
    type(psb_desc_type)            :: desc_a
    integer                        :: ictxt, info
    character                      :: afmt*5
    type(psb_d_sparse_mat)   :: a
    real(psb_dpk_)           :: zt(nb),glob_x,glob_y,glob_z
    integer                  :: m,n,nnz,glob_row,nlr,i,ii,ib,k
    integer                  :: x,y,z,ia,indx_owner
    integer                  :: np, iam, nr, nt,nz,isz
    integer                  :: element
    integer, allocatable     :: irow(:),icol(:),myidx(:)
    real(psb_dpk_), allocatable :: val(:)
    type(psb_d_coo_sparse_mat) :: acoo
    type(psb_d_csr_sparse_mat) :: acsr
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_dpk_)         :: deltah
    real(psb_dpk_),parameter   :: rhs=0.d0,one=1.d0,zero=0.d0
    real(psb_dpk_)   :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcpy, tmov
    real(psb_dpk_)   :: a1, a2, a3, a4, b1, b2, b3 
    external         :: a1, a2, a3, a4, b1, b2, b3
    integer          :: err_act

    character(len=20)  :: name, ch_err, asbfmt

    info = psb_success_
    name = 'creatae_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    deltah = 1.d0/(idim-1)

    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes 

    m   = idim*idim*idim
    n   = m
    nnz = ((n*9)/(np))
    if(iam == psb_root_) write(psb_err_unit,'("Generating Matrix (size=",i0,")...")')n

    !
    ! Using a simple BLOCK distribution.
    !
    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))

    nt = nr
    call psb_sum(ictxt,nt) 
    if (nt /= m) write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
    write(psb_err_unit,*) iam, 'Initialization ',nr,nt,m
    nlr = nt
    call psb_barrier(ictxt)

    call acoo%set_null()
    t0 = psb_wtime()

    call acoo%allocate(nr,nr)

    talc = psb_wtime()-t0

!!$    write(psb_out_unit,*) 'Test get size:',d_coo_get_size(acoo)
!!$    write(psb_out_unit,*) 'Test 2 get size:',acoo%get_size(),acoo%get_nzeros()

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    ! we build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    allocate(val(20*nb),irow(20*nb),&
         &icol(20*nb),myidx(nlr),stat=info)
    if (info /= psb_success_ ) then 
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    ! loop over rows belonging to current process in a block
    ! distribution.

    call psb_barrier(ictxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1) 
!!$      write(psb_err_unit,*) 'Row ',ii,ib
      element = 1
      do k=1,ib
        i=ii+k-1
        ! local matrix pointer 
        glob_row=i
        ! compute gridpoint coordinates
        if (mod(glob_row,(idim*idim)) == 0) then
          x = glob_row/(idim*idim)
        else
          x = glob_row/(idim*idim)+1
        endif
        if (mod((glob_row-(x-1)*idim*idim),idim) == 0) then
          y = (glob_row-(x-1)*idim*idim)/idim
        else
          y = (glob_row-(x-1)*idim*idim)/idim+1
        endif
        z = glob_row-(x-1)*idim*idim-(y-1)*idim
        ! glob_x, glob_y, glob_x coordinates
        glob_x=x*deltah
        glob_y=y*deltah
        glob_z=z*deltah

        ! check on boundary points 
        zt(k) = 0.d0
        ! internal point: build discretization
        !   
        !  term depending on   (x-1,y,z)
        !
        if (x == 1) then 
          val(element)=-b1(glob_x,glob_y,glob_z)&
               & -a1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*(-val(element))
        else
          val(element)=-b1(glob_x,glob_y,glob_z)&
               & -a1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-2)*idim*idim+(y-1)*idim+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y-1,z)
        if (y == 1) then 
          val(element)=-b2(glob_x,glob_y,glob_z)&
               & -a2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
        else
          val(element)=-b2(glob_x,glob_y,glob_z)&
               & -a2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*idim*idim+(y-2)*idim+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z-1)
        if (z == 1) then 
          val(element)=-b3(glob_x,glob_y,glob_z)&
               & -a3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
        else
          val(element)=-b3(glob_x,glob_y,glob_z)&
               & -a3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*idim*idim+(y-1)*idim+(z-1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y,z)
        val(element)=2*b1(glob_x,glob_y,glob_z)&
             & +2*b2(glob_x,glob_y,glob_z)&
             & +2*b3(glob_x,glob_y,glob_z)&
             & +a1(glob_x,glob_y,glob_z)&
             & +a2(glob_x,glob_y,glob_z)&
             & +a3(glob_x,glob_y,glob_z)
        val(element) = val(element)/(deltah*&
             & deltah)
        icol(element) = (x-1)*idim*idim+(y-1)*idim+(z)
        irow(element) = glob_row
        element       = element+1                  
        !  term depending on     (x,y,z+1)
        if (z == idim) then 
          val(element)=-b1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
        else
          val(element)=-b1(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*idim*idim+(y-1)*idim+(z+1)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x,y+1,z)
        if (y == idim) then 
          val(element)=-b2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          zt(k) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
        else
          val(element)=-b2(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x-1)*idim*idim+(y)*idim+(z)
          irow(element) = glob_row
          element       = element+1
        endif
        !  term depending on     (x+1,y,z)
        if (x<idim) then 
          val(element)=-b3(glob_x,glob_y,glob_z)
          val(element) = val(element)/(deltah*&
               & deltah)
          icol(element) = (x)*idim*idim+(y-1)*idim+(z)
          irow(element) = glob_row
          element       = element+1
        endif

      end do
      call acoo%csput(element-1,irow,icol,val,1,nr,1,nr,info)

    end do

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    
!!$    call acoo%print(19)
    t1 = psb_wtime()
!!$    write(psb_err_unit,*) 'out of loop ',acoo%get_nzeros()
    call acoo%fix(info)
!!$    write(psb_err_unit,*) '2 out of loop ',acoo%get_nzeros()

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    tasb = psb_wtime()-t1
!!$    call acoo%print(20)
    t1 = psb_wtime()    
    call acsr%cp_from_coo(acoo,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='cp rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    tcpy = psb_wtime()-t1
!!$    call acsr%print(21)
    t1 = psb_wtime()    
    call acsr%mv_from_coo(acoo,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='mv rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    tmov = psb_wtime()-t1
!!$    call acsr%print(22)
    if(iam == psb_root_) then
      asbfmt = acsr%get_fmt()
      write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
           &   asbfmt
      write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
      write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
      write(psb_out_unit,'("-assembly    time : ",es12.5)') tasb
      write(psb_out_unit,'("-copy        time : ",es12.5)') tcpy
      write(psb_out_unit,'("-move        time : ",es12.5)') tmov
!!$      write(psb_out_unit,'("-total       time : ",es12.5)') ttot

    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine create_matrix
end program d_coo_matgen
!
! functions parametrizing the differential equation 
!  
function a1(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) :: a1
  real(psb_dpk_) :: x,y,z
  a1=1.d0
end function a1
function a2(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  a2
  real(psb_dpk_) :: x,y,z
  a2=2.d1*y
end function a2
function a3(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  a3
  real(psb_dpk_) :: x,y,z      
  a3=1.d0
end function a3
function a4(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  a4
  real(psb_dpk_) :: x,y,z      
  a4=1.d0
end function a4
function b1(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  b1   
  real(psb_dpk_) :: x,y,z
  b1=1.d0
end function b1
function b2(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  b2
  real(psb_dpk_) :: x,y,z
  b2=1.d0
end function b2
function b3(x,y,z)
  use psb_sparse_mod, only : psb_dpk_
  real(psb_dpk_) ::  b3
  real(psb_dpk_) :: x,y,z
  b3=1.d0
end function b3


