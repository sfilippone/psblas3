!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
! File: ppde90.f90
!
! Program: ppde90
! This sample program shows how to build and solve a sparse linear
!
! The program  solves a linear system based on the partial differential
! equation 
!
! 
!
! The equation generated is
!
!   b1 d d (u)  b2  d d (u)    a1 d (u))  a2 d (u)))   
! -   ------   -    ------  +    ----- +  ------     + a3 u = 0
!      dx dx         dy dy         dx        dy        
!
! 
! with  Dirichlet boundary conditions on the unit cube 
!
!    0<=x,y,z<=1
! 
! The equation is discretized with finite differences and uniform stepsize;
! the resulting  discrete  equation is
!
! ( u(x,y,z)(2b1+2b2+a1+a2)+u(x-1,y)(-b1-a1)+u(x,y-1)(-b2-a2)+
!  -u(x+1,y)b1-u(x,y+1)b2)*(1/h**2)
!
! Example taken from: C.T.Kelley
!    Iterative Methods for Linear and Nonlinear Equations
!    SIAM 1995
!
!
! In this sample program the index space of the discretized
! computational domain is first numbered sequentially in a standard way, 
! then the corresponding vector is distributed according to an HPF BLOCK
! distribution directive.
!
! Boundary conditions are set in a very simple way, by adding 
! equations of the form
!
!   u(x,y) = rhs(x,y)
!
program pde90
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  implicit none

  ! input parameters
  character :: cmethd*10, prec*10, afmt*5
  integer      :: idim, iret

  ! miscellaneous 
  integer              :: iargc,convert_descr,dim, check_descr
  real(kind(1.d0)), parameter :: one = 1.d0
  real(kind(1.d0)) :: t1, t2, tprec, tsolve, t3, t4 

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a,  l, u, h
  type(psb_dprec_type)  :: pre
  ! descriptor
  type(psb_desc_type)   :: desc_a, desc_a_out
  ! dense matrices
  real(kind(1.d0)), allocatable :: b(:), x(:), d(:),ld(:)
  integer, allocatable :: work(:)
  ! blacs parameters
  integer            :: ictxt, iam, np

  ! solver parameters
  integer            :: iter, itmax,ierr,itrace, methd,iprec, istopc,&
       & iparm(20), ml, novr
  real(kind(1.d0))   :: err, eps, rparm(20)

  ! other variables
  integer            :: i,info
  integer            :: internal, m,ii
  character(len=10)  :: ptype
  character(len=20)  :: name,ch_err

  if(psb_get_errstatus().ne.0) goto 9999
  info=0
  name='pde90'
  call psb_set_errverbosity(2)
  call psb_set_erraction(0)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif

  !
  !  get parameters
  !
  call get_parms(ictxt,cmethd,iprec,novr,afmt,idim,istopc,itmax,itrace,ml)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call create_matrix(idim,a,b,x,desc_a,part_block,ictxt,afmt,info)  
  t2 = psb_wtime() - t1
  if(info.ne.0) then
    info=4010
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_amx(ictxt,t2)
  if (iam == 0) write(*,'("Overall matrix creation time : ",es10.4)')t2
  if (iam == 0) write(*,'(" ")')
  !
  !  prepare the preconditioner.
  !  

  if(iam == psb_root_) write(0,'("Setting preconditioner to : ",a)')pr_to_str(iprec)
  select case(iprec)
  case(noprec_)
    call psb_precset(pre,'noprec',info)
  case(diagsc_)             
    call psb_precset(pre,'diagsc',info)
  case(bja_)             
    call psb_precset(pre,'ilu',info)
  case default
    call psb_precset(pre,'ilu',info)
  end select

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_precbld(a,desc_a,pre,info)
  if(info.ne.0) then
    info=4010
    ch_err='psb_precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  tprec = psb_wtime()-t1

  call psb_amx(ictxt,tprec)

  if (iam == 0) write(*,'("Preconditioner time : ",es10.4)')tprec
  if (iam == 0) write(*,'(" ")')

  !
  ! iterative method parameters 
  !
  if(iam == psb_root_) write(*,'("Calling iterative method ",a)')cmethd
  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  eps   = 1.d-9
  call psb_krylov(cmethd,a,pre,b,x,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=ml)     

  if(info.ne.0) then
    info=4010
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2)

  if (iam == 0) then
    write(*,'(" ")')
    write(*,'("Time to solve matrix : ",es10.4)')t2
    write(*,'("Time per iteration   : ",es10.4)')t2/iter
    write(*,'("Number of iterations : ",i0)')iter
    write(*,'("Error on exit        : ",es10.4)')err
    write(*,'("Info  on exit        : ",i0)')info
  end if

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call psb_precfree(pre,info)
  call psb_cdfree(desc_a,info)
  if(info.ne.0) then
    info=4010
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  !
  ! get iteration parameters from the command line
  !
  subroutine  get_parms(ictxt,cmethd,iprec,novr,afmt,idim,istopc,itmax,itrace,ml)
    integer      :: ictxt
    character    :: cmethd*10, afmt*5
    integer      :: idim, iret, istopc,itmax,itrace,ml, iprec, novr
    character*40 :: charbuf
    integer      :: iargc, np, iam
    external     iargc
    integer      :: intbuf(10), ip

    call psb_info(ictxt, iam, np)

    if (iam==0) then
      read(*,*) ip
      if (ip.ge.3) then
        read(*,*) cmethd
        read(*,*) iprec
        read(*,*) novr
        read(*,*) afmt

        ! convert strings in array
        do i = 1, len(cmethd)
          intbuf(i) = iachar(cmethd(i:i))
        end do
        ! broadcast parameters to all processors
        call psb_bcast(ictxt,intbuf(1:10))

        ! broadcast parameters to all processors
        call psb_bcast(ictxt,iprec)
        call psb_bcast(ictxt,novr)

        do i = 1, len(afmt)
          intbuf(i) = iachar(afmt(i:i))
        end do
        call psb_bcast(ictxt,intbuf(1:10))

        read(*,*) idim
        if (ip.ge.4) then
          read(*,*) istopc
        else
          istopc=1        
        endif
        if (ip.ge.5) then
          read(*,*) itmax
        else
          itmax=500
        endif
        if (ip.ge.6) then
          read(*,*) itrace
        else
          itrace=-1
        endif
        if (ip.ge.7) then
          read(*,*) ml
        else
          ml=1
        endif
        ! broadcast parameters to all processors    

        intbuf(1) = idim
        intbuf(2) = istopc
        intbuf(3) = itmax
        intbuf(4) = itrace
        intbuf(5) = ml
        call psb_bcast(ictxt,intbuf(1:5))

        write(*,'("Solving matrix       : ell1")')      
        write(*,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
        write(*,'("Number of processors : ",i0)')np
        write(*,'("Data distribution    : BLOCK")')
        write(*,'("Preconditioner       : ",a)')pr_to_str(iprec)
        if(iprec.gt.2) write(*,'("Overlapping levels   : ",i0)')novr
        write(*,'("Iterative method     : ",a)')cmethd
        write(*,'(" ")')
      else
        ! wrong number of parameter, print an error message and exit
        call pr_usage(0)      
        call psb_abort(ictxt)
        stop 1
      endif
    else
      call psb_bcast(ictxt,intbuf(1:10))

      do i = 1, 10
        cmethd(i:i) = achar(intbuf(i))
      end do

      call psb_bcast(ictxt,iprec)
      call psb_bcast(ictxt,novr)
      call psb_bcast(ictxt,intbuf(1:10))
      do i = 1, 5
        afmt(i:i) = achar(intbuf(i))
      end do
      call psb_bcast(ictxt,intbuf(1:5))
      idim    = intbuf(1)
      istopc  = intbuf(2)
      itmax   = intbuf(3)
      itrace  = intbuf(4)
      ml      = intbuf(5)
    end if
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
    write(iout,*)'     methd:    cgstab tfqmr cgs' 
    write(iout,*)'     prec :    ilu diagsc none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2 or 3 [1]  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    0  (no tracing, default) or '  
    write(iout,*)'               >= 0 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage

  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs. 
  !
  subroutine create_matrix(idim,a,b,t,desc_a,parts,ictxt,afmt,info)
    !
    !   discretize the partial diferential equation
    ! 
    !   b1 dd(u)  b2 dd(u)    b3 dd(u)    a1 d(u)   a2 d(u)  a3 d(u)  
    ! -   ------ -  ------ -  ------ -  -----  -  ------  -  ------ + a4 u 
    !      dxdx     dydy       dzdz        dx       dy         dz   
    !
    !  = 0 
    ! 
    ! boundary condition: dirichlet
    !    0< x,y,z<1
    !  
    !  u(x,y,z)(2b1+2b2+2b3+a1+a2+a3)+u(x-1,y,z)(-b1-a1)+u(x,y-1,z)(-b2-a2)+
    !  + u(x,y,z-1)(-b3-a3)-u(x+1,y,z)b1-u(x,y+1,z)b2-u(x,y,z+1)b3

    use psb_base_mod
    implicit none
    integer                        :: idim
    integer, parameter             :: nbmax=10
    real(kind(1.d0)), allocatable  :: b(:),t(:)
    type(psb_desc_type)            :: desc_a
    integer                        :: ictxt, info
    character                      :: afmt*5
    interface 
      !   .....user passed subroutine.....
      subroutine parts(global_indx,n,np,pv,nv)
        implicit none
        integer, intent(in)  :: global_indx, n, np
        integer, intent(out) :: nv
        integer, intent(out) :: pv(*) 
      end subroutine parts
    end interface   ! local variables
    type(psb_dspmat_type)    :: a
    real(kind(1.d0))         :: zt(nbmax),glob_x,glob_y,glob_z
    integer                  :: m,n,nnz,glob_row,j
    integer                  :: x,y,z,counter,ia,i,indx_owner
    integer                  :: np, iam
    integer                  :: element
    integer                  :: nv, inv
    integer, allocatable     :: irow(:),icol(:)
    real(kind(1.d0)), allocatable :: val(:)
    integer, allocatable     :: prv(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(kind(1.d0))         :: deltah
    real(kind(1.d0)),parameter   :: rhs=0.d0,one=1.d0,zero=0.d0
    real(kind(1.d0))   :: t1, t2, t3, tins, tasb
    real(kind(1.d0))   :: a1, a2, a3, a4, b1, b2, b3 
    external           :: a1, a2, a3, a4, b1, b2, b3
    integer            :: nb, ir1, ir2, ipr, err_act
    logical            :: own
    ! common area

    character(len=20)  :: name, ch_err

    info = 0
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    deltah = 1.d0/(idim-1)

    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes 

    m   = idim*idim*idim
    n   = m
    nnz = ((n*9)/(np))
    if(iam == psb_root_) write(0,'("Generating Matrix (size=",i0x,")...")')n

    call psb_cdall(ictxt,desc_a,info,mg=n,parts=parts)
    call psb_spall(a,desc_a,info,nnz=nnz)
    ! define  rhs from boundary conditions; also build initial guess 
    call psb_geall(b,desc_a,info)
    call psb_geall(t,desc_a,info)
    if(info.ne.0) then
      info=4010
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    ! we build an auxiliary matrix consisting of one row at a
    ! time; just a small matrix. might be extended to generate 
    ! a bunch of rows per call. 
    ! 
    allocate(val(20*nbmax),irow(20*nbmax),&
         &icol(20*nbmax),prv(np),stat=info)
    if (info.ne.0 ) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    endif

    tins = 0.d0
    call psb_barrier(ictxt)
    t1 = psb_wtime()

    ! loop over rows belonging to current process in a block
    ! distribution.

    !    icol(1)=1    
    do glob_row = 1, n
      call parts(glob_row,n,np,prv,nv)
      do inv = 1, nv
        indx_owner = prv(inv)
        if (indx_owner == iam) then
          ! local matrix pointer 
          element=1
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
          zt(1) = 0.d0
          ! internal point: build discretization
          !   
          !  term depending on   (x-1,y,z)
          !
          if (x==1) then 
            val(element)=-b1(glob_x,glob_y,glob_z)&
                 & -a1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*(-val(element))
          else
            val(element)=-b1(glob_x,glob_y,glob_z)&
                 & -a1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-2)*idim*idim+(y-1)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x,y-1,z)
          if (y==1) then 
            val(element)=-b2(glob_x,glob_y,glob_z)&
                 & -a2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b2(glob_x,glob_y,glob_z)&
                 & -a2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-2)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x,y,z-1)
          if (z==1) then 
            val(element)=-b3(glob_x,glob_y,glob_z)&
                 & -a3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b3(glob_x,glob_y,glob_z)&
                 & -a3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-1)*idim+(z-1)
            element=element+1
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
          icol(element)=(x-1)*idim*idim+(y-1)*idim+(z)
          element=element+1                  
          !  term depending on     (x,y,z+1)
          if (z==idim) then 
            val(element)=-b1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b1(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y-1)*idim+(z+1)
            element=element+1
          endif
          !  term depending on     (x,y+1,z)
          if (y==idim) then 
            val(element)=-b2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            zt(1) = exp(-glob_y**2-glob_z**2)*exp(-glob_x)*(-val(element))  
          else
            val(element)=-b2(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x-1)*idim*idim+(y)*idim+(z)
            element=element+1
          endif
          !  term depending on     (x+1,y,z)
          if (x<idim) then 
            val(element)=-b3(glob_x,glob_y,glob_z)
            val(element) = val(element)/(deltah*&
                 & deltah)
            icol(element)=(x)*idim*idim+(y-1)*idim+(z)
            element=element+1
          endif
          irow(1:element-1)=glob_row
          ia=glob_row

          t3 = psb_wtime()
          call psb_spins(element-1,irow,icol,val,a,desc_a,info)
          if(info.ne.0) exit
          tins = tins + (psb_wtime()-t3)
          call psb_geins(1,(/ia/),zt(1:1),b,desc_a,info)
          if(info.ne.0) exit
          zt(1)=0.d0
          call psb_geins(1,(/ia/),zt(1:1),t,desc_a,info)
          if(info.ne.0) exit
        end if
      end do
    end do

    call psb_barrier(ictxt)    
    t2 = psb_wtime()-t1

    if(info.ne.0) then
      info=4010
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(val,irow,icol)

    t1 = psb_wtime()
    call psb_cdasb(desc_a,info)
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    call psb_barrier(ictxt)
    tasb = psb_wtime()-t1
    if(info.ne.0) then
      info=4010
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_amx(ictxt,t2)
    call psb_amx(ictxt,tins)
    call psb_amx(ictxt,tasb)

    if(iam == psb_root_) then
      write(*,'("The matrix has been generated and assembeld in ",a3," format.")')a%fida(1:3)
      write(*,'("-pspins time   : ",es10.4)')tins
      write(*,'("-insert time   : ",es10.4)')t2
      write(*,'("-assembly time : ",es10.4)')tasb
    end if

    call psb_geasb(b,desc_a,info)
    call psb_geasb(t,desc_a,info)
    if(info.ne.0) then
      info=4010
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
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
end program pde90
!
! functions parametrizing the differential equation 
!  
function a1(x,y,z)
  real(kind(1.d0)) :: a1
  real(kind(1.d0)) :: x,y,z
  a1=1.d0
end function a1
function a2(x,y,z)
  real(kind(1.d0)) ::  a2
  real(kind(1.d0)) :: x,y,z
  a2=2.d1*y
end function a2
function a3(x,y,z)
  real(kind(1.d0)) ::  a3
  real(kind(1.d0)) :: x,y,z      
  a3=1.d0
end function a3
function a4(x,y,z)
  real(kind(1.d0)) ::  a4
  real(kind(1.d0)) :: x,y,z      
  a4=1.d0
end function a4
function b1(x,y,z)
  real(kind(1.d0)) ::  b1   
  real(kind(1.d0)) :: x,y,z
  b1=1.d0
end function b1
function b2(x,y,z)
  real(kind(1.d0)) ::  b2
  real(kind(1.d0)) :: x,y,z
  b2=1.d0
end function b2
function b3(x,y,z)
  real(kind(1.d0)) ::  b3
  real(kind(1.d0)) :: x,y,z
  b3=1.d0
end function b3


