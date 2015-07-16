!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
!
!  subroutine to allocate and fill in the coefficient matrix and
!  the rhs. 
!
subroutine psb_s_gen_pde3d(ictxt,idim,a,bv,xv,desc_a,afmt,&
     & a1,a2,a3,b1,b2,b3,c,g,info,f,amold,vmold,imold,nrl)
  use psb_base_mod
  use psb_s_genpde_mod, psb_protect_name => psb_s_gen_pde3d
  !
  !   Discretizes the partial differential equation
  ! 
  !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
  ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
  !      dxdx     dydy       dzdz        dx       dy         dz   
  !
  ! with Dirichlet boundary conditions
  !   u = g 
  !
  !  on the unit cube  0<=x,y,z<=1.
  !
  !
  ! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
  !
  implicit none
  procedure(s_func_3d)  :: b1,b2,b3,c,a1,a2,a3,g
  integer(psb_ipk_)     :: idim
  type(psb_sspmat_type) :: a
  type(psb_s_vect_type) :: xv,bv
  type(psb_desc_type)   :: desc_a
  integer(psb_ipk_)     :: ictxt, info
  character(len=*)      :: afmt
  procedure(s_func_3d), optional :: f
  class(psb_s_base_sparse_mat), optional :: amold
  class(psb_s_base_vect_type), optional :: vmold
  class(psb_i_base_vect_type), optional :: imold
  integer(psb_ipk_), optional :: nrl

  ! Local variables.

  integer(psb_ipk_), parameter :: nb=20
  type(psb_s_csc_sparse_mat)  :: acsc
  type(psb_s_coo_sparse_mat)  :: acoo
  type(psb_s_csr_sparse_mat)  :: acsr
  real(psb_spk_)           :: zt(nb),x,y,z
  integer(psb_ipk_) :: m,n,nnz,glob_row,nlr,i,ii,ib,k
  integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
  integer(psb_ipk_) :: np, iam, nr, nt
  integer(psb_ipk_) :: icoeff
  integer(psb_ipk_), allocatable     :: irow(:),icol(:),myidx(:)
  real(psb_spk_), allocatable :: val(:)
  ! deltah dimension of each grid cell
  ! deltat discretization time
  real(psb_spk_)            :: deltah, sqdeltah, deltah2
  real(psb_spk_), parameter :: rhs=0.e0,one=1.e0,zero=0.e0
  real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
  integer(psb_ipk_) :: err_act
  procedure(s_func_3d), pointer :: f_
  character(len=20)  :: name, ch_err,tmpfmt

  info = psb_success_
  name = 'create_matrix'
  call psb_erractionsave(err_act)

  call psb_info(ictxt, iam, np)

  
  if (present(f)) then 
    f_ => f
  else
    f_ => s_null_func_3d
  end if

  deltah   = 1.e0/(idim+2)
  sqdeltah = deltah*deltah
  deltah2  = 2.e0* deltah

  ! initialize array descriptor and sparse matrix storage. provide an
  ! estimate of the number of non zeroes 

  m   = idim*idim*idim
  n   = m
  nnz = ((n*9)/(np))
  if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n

  if (present(nrl)) then 
    nr = nrl
  else
    !
    ! Using a simple BLOCK distribution.
    !
    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))
  end if

  nt = nr
  call psb_sum(ictxt,nt) 
  if (nt /= m) then 
    write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
    info = -1
    call psb_barrier(ictxt)
    call psb_abort(ictxt)
    return    
  end if
    
  call psb_barrier(ictxt)
  t0 = psb_wtime()
  call psb_cdall(ictxt,desc_a,info,nl=nr)
  if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
  ! define  rhs from boundary conditions; also build initial guess 
  if (info == psb_success_) call psb_geall(xv,desc_a,info)
  if (info == psb_success_) call psb_geall(bv,desc_a,info)

  call psb_barrier(ictxt)
  talc = psb_wtime()-t0

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
       &icol(20*nb),stat=info)
  if (info /= psb_success_ ) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  myidx = desc_a%get_global_indices()
  nlr = size(myidx)

  ! loop over rows belonging to current process in a block
  ! distribution.

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  do ii=1, nlr,nb
    ib = min(nb,nlr-ii+1) 
    icoeff = 1
    do k=1,ib
      i=ii+k-1
      ! local matrix pointer 
      glob_row=myidx(i)
      ! compute gridpoint coordinates
      if (mod(glob_row,(idim*idim)) == 0) then
        ix = glob_row/(idim*idim)
      else
        ix = glob_row/(idim*idim)+1
      endif
      if (mod((glob_row-(ix-1)*idim*idim),idim) == 0) then
        iy = (glob_row-(ix-1)*idim*idim)/idim
      else
        iy = (glob_row-(ix-1)*idim*idim)/idim+1
      endif
      iz = glob_row-(ix-1)*idim*idim-(iy-1)*idim
      ! x, y, x coordinates
      x = ix*deltah
      y = iy*deltah
      z = iz*deltah
      zt(k) = f_(x,y,z)
      ! internal point: build discretization
      !   
      !  term depending on   (x-1,y,z)
      !
      val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
      if (ix == 1) then 
        zt(k) = g(szero,y,z)*(-val(icoeff)) + zt(k)
      else
        icol(icoeff) = (ix-2)*idim*idim+(iy-1)*idim+(iz)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x,y-1,z)
      val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
      if (iy == 1) then 
        zt(k) = g(x,szero,z)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim*idim+(iy-2)*idim+(iz)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x,y,z-1)
      val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
      if (iz == 1) then 
        zt(k) = g(x,y,szero)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz-1)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif

      !  term depending on     (x,y,z)
      val(icoeff)=2.e0*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah &
           & + c(x,y,z)
      icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz)
      irow(icoeff) = glob_row
      icoeff       = icoeff+1                  
      !  term depending on     (x,y,z+1)
      val(icoeff)=-a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
      if (iz == idim) then 
        zt(k) = g(x,y,sone)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz+1)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x,y+1,z)
      val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
      if (iy == idim) then 
        zt(k) = g(x,sone,z)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim*idim+(iy)*idim+(iz)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x+1,y,z)
      val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
      if (ix==idim) then 
        zt(k) = g(sone,y,z)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix)*idim*idim+(iy-1)*idim+(iz)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif

    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
    if(info /= psb_success_) exit
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
    if(info /= psb_success_) exit
    zt(:)=0.e0
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
    if(info /= psb_success_) exit
  end do

  tgen = psb_wtime()-t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='insert rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  deallocate(val,irow,icol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_cdasb(desc_a,info,mold=imold)
  tcdasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  if (info == psb_success_) then 
    if (present(amold)) then 
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,mold=amold)
    else
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    end if
  end if
  call psb_barrier(ictxt)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (info == psb_success_) call psb_geasb(xv,desc_a,info,mold=vmold)
  if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  tasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  ttot = psb_wtime() - t0 

  call psb_amx(ictxt,talc)
  call psb_amx(ictxt,tgen)
  call psb_amx(ictxt,tasb)
  call psb_amx(ictxt,ttot)
  if(iam == psb_root_) then
    tmpfmt = a%get_fmt()
    write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
         &   tmpfmt
    write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
    write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
    write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
    write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
    write(psb_out_unit,'("-total       time : ",es12.5)') ttot

  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_s_gen_pde3d



!
!  subroutine to allocate and fill in the coefficient matrix and
!  the rhs. 
!
subroutine psb_s_gen_pde2d(ictxt,idim,a,bv,xv,desc_a,afmt,&
     & a1,a2,b1,b2,c,g,info,f,amold,vmold,imold,nrl)
  use psb_base_mod
  use psb_s_genpde_mod, psb_protect_name => psb_s_gen_pde2d
  !
  !   Discretizes the partial differential equation
  ! 
  !   a1 dd(u)  a2 dd(u)    b1 d(u)  b2 d(u) 
  ! -   ------ -  ------ +  -----  +  ------  + c u = f
  !      dxdx     dydy         dx       dy     
  !
  ! with Dirichlet boundary conditions
  !   u = g 
  !
  !  on the unit square  0<=x,y<=1.
  !
  !
  ! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
  !
  implicit none
  procedure(s_func_2d)  :: b1,b2,c,a1,a2,g
  integer(psb_ipk_)     :: idim
  type(psb_sspmat_type) :: a
  type(psb_s_vect_type) :: xv,bv
  type(psb_desc_type)   :: desc_a
  integer(psb_ipk_)     :: ictxt, info
  character(len=*)      :: afmt
  procedure(s_func_2d), optional :: f
  class(psb_s_base_sparse_mat), optional :: amold
  class(psb_s_base_vect_type), optional :: vmold
  class(psb_i_base_vect_type), optional :: imold
  integer(psb_ipk_), optional :: nrl

  ! Local variables.

  integer(psb_ipk_), parameter :: nb=20
  type(psb_s_csc_sparse_mat)  :: acsc
  type(psb_s_coo_sparse_mat)  :: acoo
  type(psb_s_csr_sparse_mat)  :: acsr
  real(psb_spk_)           :: zt(nb),x,y,z
  integer(psb_ipk_) :: m,n,nnz,glob_row,nlr,i,ii,ib,k
  integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
  integer(psb_ipk_) :: np, iam, nr, nt
  integer(psb_ipk_) :: icoeff
  integer(psb_ipk_), allocatable     :: irow(:),icol(:),myidx(:)
  real(psb_spk_), allocatable :: val(:)
  ! deltah dimension of each grid cell
  ! deltat discretization time
  real(psb_spk_)            :: deltah, sqdeltah, deltah2
  real(psb_spk_), parameter :: rhs=0.e0,one=1.e0,zero=0.e0
  real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
  integer(psb_ipk_) :: err_act
  procedure(s_func_2d), pointer :: f_
  character(len=20)  :: name, ch_err,tmpfmt

  info = psb_success_
  name = 'create_matrix'
  call psb_erractionsave(err_act)

  call psb_info(ictxt, iam, np)

  
  if (present(f)) then 
    f_ => f
  else
    f_ => s_null_func_2d
  end if

  deltah   = 1.e0/(idim+2)
  sqdeltah = deltah*deltah
  deltah2  = 2.e0* deltah

  ! initialize array descriptor and sparse matrix storage. provide an
  ! estimate of the number of non zeroes 

  m   = idim*idim
  n   = m
  nnz = ((n*7)/(np))
  if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n

  if (present(nrl)) then 
    nr = nrl
  else
    !
    ! Using a simple BLOCK distribution.
    !
    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))
  end if
  
  nt = nr
  call psb_sum(ictxt,nt) 
  if (nt /= m) then 
    write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
    info = -1
    call psb_barrier(ictxt)
    call psb_abort(ictxt)
    return    
  end if
  call psb_barrier(ictxt)
  t0 = psb_wtime()
  call psb_cdall(ictxt,desc_a,info,nl=nr)
  if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
  ! define  rhs from boundary conditions; also build initial guess 
  if (info == psb_success_) call psb_geall(xv,desc_a,info)
  if (info == psb_success_) call psb_geall(bv,desc_a,info)
  call psb_barrier(ictxt)
  talc = psb_wtime()-t0

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
       &icol(20*nb),stat=info)
  if (info /= psb_success_ ) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  
  myidx = desc_a%get_global_indices()
  nlr = size(myidx)

  ! loop over rows belonging to current process in a block
  ! distribution.

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  do ii=1, nlr,nb
    ib = min(nb,nlr-ii+1) 
    icoeff = 1
    do k=1,ib
      i=ii+k-1
      ! local matrix pointer 
      glob_row=myidx(i)
      ! compute gridpoint coordinates
      if (mod(glob_row,(idim)) == 0) then
        ix = glob_row/(idim)
      else
        ix = glob_row/(idim)+1
      endif
      iy = (glob_row-(ix-1)*idim)
      ! x, y
      x = ix*deltah
      y = iy*deltah

      zt(k) = f_(x,y)
      ! internal point: build discretization
      !   
      !  term depending on   (x-1,y)
      !
      val(icoeff) = -a1(x,y)/sqdeltah-b1(x,y)/deltah2
      if (ix == 1) then 
        zt(k) = g(szero,y)*(-val(icoeff)) + zt(k)
      else
        icol(icoeff) = (ix-2)*idim+iy
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x,y-1)
      val(icoeff)  = -a2(x,y)/sqdeltah-b2(x,y)/deltah2
      if (iy == 1) then 
        zt(k) = g(x,szero)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim+(iy-1)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif

      !  term depending on     (x,y)
      val(icoeff)=2.e0*(a1(x,y) + a2(x,y))/sqdeltah + c(x,y)
      icol(icoeff) = (ix-1)*idim+iy
      irow(icoeff) = glob_row
      icoeff       = icoeff+1                  
      !  term depending on     (x,y+1)
      val(icoeff)=-a2(x,y)/sqdeltah+b2(x,y)/deltah2
      if (iy == idim) then 
        zt(k) = g(x,sone)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix-1)*idim+(iy+1)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif
      !  term depending on     (x+1,y)
      val(icoeff)=-a1(x,y)/sqdeltah+b1(x,y)/deltah2
      if (ix==idim) then 
        zt(k) = g(sone,y)*(-val(icoeff))   + zt(k)
      else
        icol(icoeff) = (ix)*idim+(iy)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
      endif

    end do
    call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
    if(info /= psb_success_) exit
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
    if(info /= psb_success_) exit
    zt(:)=0.e0
    call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
    if(info /= psb_success_) exit
  end do

  tgen = psb_wtime()-t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='insert rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  deallocate(val,irow,icol)

  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_cdasb(desc_a,info,mold=imold)
  tcdasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  if (info == psb_success_) then 
    if (present(amold)) then 
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,mold=amold)
    else
      call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    end if
  end if
  call psb_barrier(ictxt)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (info == psb_success_) call psb_geasb(xv,desc_a,info,mold=vmold)
  if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  tasb = psb_wtime()-t1
  call psb_barrier(ictxt)
  ttot = psb_wtime() - t0 

  call psb_amx(ictxt,talc)
  call psb_amx(ictxt,tgen)
  call psb_amx(ictxt,tasb)
  call psb_amx(ictxt,ttot)
  if(iam == psb_root_) then
    tmpfmt = a%get_fmt()
    write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")')&
         &   tmpfmt
    write(psb_out_unit,'("-allocation  time : ",es12.5)') talc
    write(psb_out_unit,'("-coeff. gen. time : ",es12.5)') tgen
    write(psb_out_unit,'("-desc asbly  time : ",es12.5)') tcdasb
    write(psb_out_unit,'("- mat asbly  time : ",es12.5)') tasb
    write(psb_out_unit,'("-total       time : ",es12.5)') ttot

  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psb_s_gen_pde2d
