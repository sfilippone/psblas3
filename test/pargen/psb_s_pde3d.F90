!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: psb_s_pde3d.f90
!
! Program: psb_s_pde3d
! This sample program solves a linear system obtained by discretizing a
! PDE with Dirichlet BCs.
!
!
! The PDE is a general second order equation in 3d
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
! There are three choices available for data distribution:
! 1. A simple BLOCK distribution
! 2. A ditribution based on arbitrary assignment of indices to processes,
!    typically from a graph partitioner
! 3. A 3D distribution in which the unit cube is partitioned
!    into subcubes, each one assigned to a process.
!
!
module psb_s_pde3d_mod


  use psb_base_mod, only : psb_spk_, psb_ipk_, psb_lpk_, psb_desc_type,&
       &  psb_sspmat_type, psb_s_vect_type, szero,&
       &  psb_s_base_sparse_mat, psb_s_base_vect_type, &
       &  psb_i_base_vect_type, psb_l_base_vect_type

  interface
    function s_func_3d(x,y,z) result(val)
      import :: psb_spk_
      real(psb_spk_), intent(in) :: x,y,z
      real(psb_spk_) :: val
    end function s_func_3d
  end interface

  interface psb_gen_pde3d
    module procedure  psb_s_gen_pde3d
  end interface psb_gen_pde3d

contains

  function s_null_func_3d(x,y,z) result(val)

    real(psb_spk_), intent(in) :: x,y,z
    real(psb_spk_) :: val

    val = szero

  end function s_null_func_3d
  !
  ! functions parametrizing the differential equation
  !

  !
  ! Note: b1, b2 and b3 are the coefficients of the first
  ! derivative of the unknown function. The default
  ! we apply here is to have them zero, so that the resulting
  ! matrix is symmetric/hermitian and suitable for
  ! testing with CG and FCG.
  ! When testing methods for non-hermitian matrices you can
  ! change the B1/B2/B3 functions to e.g. sone/sqrt((3*sone))
  !
  function b1(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y,z
    b1=szero
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y,z
    b2=szero
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  b3
    real(psb_spk_), intent(in) :: x,y,z
    b3=szero
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  c
    real(psb_spk_), intent(in) :: x,y,z
    c=szero
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  a1
    real(psb_spk_), intent(in) :: x,y,z
    a1=sone/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  a2
    real(psb_spk_), intent(in) :: x,y,z
    a2=sone/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  a3
    real(psb_spk_), intent(in) :: x,y,z
    a3=sone/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none
    real(psb_spk_) ::  g
    real(psb_spk_), intent(in) :: x,y,z
    g = szero
    if (x == sone) then
      g = sone
    else if (x == szero) then
      g = exp(y**2-z**2)
    end if
  end function g


  !
  !  subroutine to allocate and fill in the coefficient matrix and
  !  the rhs.
  !
  subroutine psb_s_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info,&
       & f,amold,vmold,imold,partition,nrl,iv)
    use psb_base_mod
    use psb_util_mod
#if defined(OPENMP)
    use omp_lib
#endif
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
    integer(psb_ipk_)     :: idim
    type(psb_sspmat_type) :: a
    type(psb_s_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type)   :: ctxt
    integer(psb_ipk_)     :: info
    character(len=*)      :: afmt
    procedure(s_func_3d), optional :: f
    class(psb_s_base_sparse_mat), optional :: amold
    class(psb_s_base_vect_type), optional :: vmold
    class(psb_i_base_vect_type), optional :: imold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)

    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_s_csc_sparse_mat)  :: acsc
    type(psb_s_coo_sparse_mat)  :: acoo
    type(psb_s_csr_sparse_mat)  :: acsr
    real(psb_spk_)           :: zt(nb),x,y,z
    integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_, mysz
    integer(psb_lpk_) :: m,n,glob_row,nt
    integer(psb_ipk_) :: ix,iy,iz,ia,indx_owner
    ! For 3D partition
    ! Note: integer control variables going directly into an MPI call
    ! must be 4 bytes, i.e. psb_mpk_
    integer(psb_mpk_) :: npdims(3), npp, minfo
    integer(psb_ipk_) :: npx,npy,npz, iamx,iamy,iamz,mynx,myny,mynz
    integer(psb_ipk_), allocatable :: bndx(:),bndy(:),bndz(:)
    ! Process grid
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable     :: myidx(:)
    ! deltah dimension of each grid cell
    ! deltat discretization time
    real(psb_spk_)            :: deltah, sqdeltah, deltah2
    real(psb_spk_), parameter :: rhs=szero,one=sone,zero=szero
    real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    procedure(s_func_3d), pointer :: f_
    character(len=20)  :: name, ch_err,tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)


    if (present(f)) then
      f_ => f
    else
      f_ => s_null_func_3d
    end if

    deltah   = sone/(idim+1)
    sqdeltah = deltah*deltah
    deltah2  = (2*sone)* deltah

    if (present(partition)) then
      if ((1<= partition).and.(partition <= 3)) then
        partition_ = partition
      else
        write(*,*) 'Invalid partition choice ',partition,' defaulting to 3'
        partition_ = 3
      end if
    else
      partition_ = 3
    end if

    ! initialize array descriptor and sparse matrix storage. provide an
    ! estimate of the number of non zeroes

    m   = (1_psb_lpk_*idim)*idim*idim
    n   = m
    nnz = ((n*7)/(np))
    if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n
    t0 = psb_wtime()
    select case(partition_)
    case(1)
      ! A BLOCK partition
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
      call psb_sum(ctxt,nt)
      if (nt /= m) then
        write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return
      end if

      !
      ! First example  of use of CDALL: specify for each process a number of
      ! contiguous rows
      !
      call psb_cdall(ctxt,desc_a,info,nl=nr)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(2)
      ! A  partition  defined by the user through IV

      if (present(iv)) then
        if (size(iv) /= m) then
          write(psb_err_unit,*) iam, 'Initialization error: wrong IV size',size(iv),m
          info = -1
          call psb_barrier(ctxt)
          call psb_abort(ctxt)
          return
        end if
      else
        write(psb_err_unit,*) iam, 'Initialization error: IV not present'
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
        return
      end if

      !
      ! Second example  of use of CDALL: specify for each row the
      ! process that owns it
      !
      call psb_cdall(ctxt,desc_a,info,vg=iv)
      myidx = desc_a%get_global_indices()
      nlr = size(myidx)

    case(3)
      ! A 3-dimensional partition

      ! A nifty MPI function will split the process list
      npdims = 0
#if defined(SERIAL_MPI)
      npdims = 1
#else
      call mpi_dims_create(np,3,npdims,info)
#endif
      npx = npdims(1)
      npy = npdims(2)
      npz = npdims(3)

      allocate(bndx(0:npx),bndy(0:npy),bndz(0:npz))
      ! We can reuse idx2ijk for process indices as well.
      call idx2ijk(iamx,iamy,iamz,iam,npx,npy,npz,base=0)
      ! Now let's split the 3D cube in hexahedra
      call dist1Didx(bndx,idim,npx)
      mynx = bndx(iamx+1)-bndx(iamx)
      call dist1Didx(bndy,idim,npy)
      myny = bndy(iamy+1)-bndy(iamy)
      call dist1Didx(bndz,idim,npz)
      mynz = bndz(iamz+1)-bndz(iamz)

      ! How many indices do I own?
      nlr = mynx*myny*mynz
      allocate(myidx(nlr))
      ! Now, let's generate the list of indices I own
      nr = 0
      do i=bndx(iamx),bndx(iamx+1)-1
        do j=bndy(iamy),bndy(iamy+1)-1
          do k=bndz(iamz),bndz(iamz+1)-1
            nr = nr + 1
            call ijk2idx(myidx(nr),i,j,k,idim,idim,idim)
          end do
        end do
      end do
      if (nr /= nlr) then
        write(psb_err_unit,*) iam,iamx,iamy,iamz, 'Initialization error: NR vs NLR ',&
             & nr,nlr,mynx,myny,mynz
        info = -1
        call psb_barrier(ctxt)
        call psb_abort(ctxt)
      end if

      !
      ! Third example  of use of CDALL: specify for each process
      ! the set of global indices it owns.
      !
      call psb_cdall(ctxt,desc_a,info,vl=myidx)

      !
      ! Specify process topology
      !
      block
        !
        ! Use adjcncy methods 
        ! 
        integer(psb_mpk_), allocatable :: neighbours(:)
        integer(psb_mpk_) :: cnt
        logical, parameter :: debug_adj=.true.
        if (debug_adj.and.(np > 1)) then 
          cnt = 0
          allocate(neighbours(np))
          if (iamx < npx-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx+1,iamy,iamz,npx,npy,npz,base=0)
          end if
          if (iamy < npy-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy+1,iamz,npx,npy,npz,base=0)
          end if
          if (iamz < npz-1) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz+1,npx,npy,npz,base=0)
          end if
          if (iamx >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx-1,iamy,iamz,npx,npy,npz,base=0)
          end if
          if (iamy >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy-1,iamz,npx,npy,npz,base=0)
          end if
          if (iamz >0) then
            cnt = cnt + 1 
            call ijk2idx(neighbours(cnt),iamx,iamy,iamz-1,npx,npy,npz,base=0)
          end if
          call psb_realloc(cnt, neighbours,info)
          call desc_a%set_p_adjcncy(neighbours)
          !write(0,*) iam,' Check on neighbours: ',desc_a%get_p_adjcncy()
        end if
      end block

    case default
      write(psb_err_unit,*) iam, 'Initialization error: should not get here'
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end select


    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz, &
         & bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)
    ! define  rhs from boundary conditions; also build initial guess
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    if (info == psb_success_) call psb_geall(bv,desc_a,info,&
         & bldmode=psb_matbld_remote_,dupl=psb_dupl_add_)

    call psb_barrier(ctxt)
    talc = psb_wtime()-t0

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    !$omp parallel shared(deltah,myidx,a,desc_a)
    !
    block 
      integer(psb_ipk_) :: i,j,k,ii,ib,icoeff, ix,iy,iz, ith,nth
      integer(psb_lpk_) :: glob_row
      integer(psb_lpk_), allocatable     :: irow(:),icol(:)
      real(psb_spk_), allocatable :: val(:)
      real(psb_spk_)    :: x,y,z, zt(nb)
#if defined(OPENMP)
      nth = omp_get_num_threads()
      ith = omp_get_thread_num()
#else
      nth = 1
      ith = 0
#endif
      allocate(val(20*nb),irow(20*nb),&
           &icol(20*nb),stat=info)
      if (info /= psb_success_ ) then
        info=psb_err_alloc_dealloc_
        call psb_errpush(info,name)
        !goto 9999
      endif

      !$omp  do schedule(dynamic)
      !     
      do ii=1, nlr, nb
        if(info /= psb_success_) cycle
        ib = min(nb,nlr-ii+1)
        !ib = min(nb,mysz-ii+1)
        icoeff = 1
        do k=1,ib
          i=ii+k-1
          ! local matrix pointer
          glob_row=myidx(i)
          ! compute gridpoint coordinates
          call idx2ijk(ix,iy,iz,glob_row,idim,idim,idim)
          ! x, y, z coordinates
          x = (ix-1)*deltah
          y = (iy-1)*deltah
          z = (iz-1)*deltah
          zt(k) = f_(x,y,z)
          ! internal point: build discretization
          !
          !  term depending on   (x-1,y,z)
          !
          val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
          if (ix == 1) then
            zt(k) = g(szero,y,z)*(-val(icoeff)) + zt(k)
          else
            call ijk2idx(icol(icoeff),ix-1,iy,iz,idim,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y-1,z)
          val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
          if (iy == 1) then
            zt(k) = g(x,szero,z)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy-1,iz,idim,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif
          !  term depending on     (x,y,z-1)
          val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
          if (iz == 1) then
            zt(k) = g(x,y,szero)*(-val(icoeff))   + zt(k)
          else
            call ijk2idx(icol(icoeff),ix,iy,iz-1,idim,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
          endif

        !  term depending on     (x,y,z)
        val(icoeff)=(2*sone)*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah &
             & + c(x,y,z)
        call ijk2idx(icol(icoeff),ix,iy,iz,idim,idim,idim)
        irow(icoeff) = glob_row
        icoeff       = icoeff+1
        !  term depending on     (x,y,z+1)
        val(icoeff)=-a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
        if (iz == idim) then
          zt(k) = g(x,y,sone)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy,iz+1,idim,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y+1,z)
        val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
        if (iy == idim) then
          zt(k) = g(x,sone,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy+1,iz,idim,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x+1,y,z)
        val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
        if (ix==idim) then
          zt(k) = g(sone,y,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix+1,iy,iz,idim,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif

      end do
#if defined(OPENMP)
!!$        write(0,*) omp_get_thread_num(),' Check insertion ',&
!!$             & irow(1:icoeff-1),':',icol(1:icoeff-1)
#endif
        call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
        if(info /= psb_success_) cycle
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
        if(info /= psb_success_) cycle
        zt(:)=szero
        call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),xv,desc_a,info)
        if(info /= psb_success_) cycle
      end do
      !$omp end do
      deallocate(val,irow,icol)
    end block
    !$omp end parallel

    tgen = psb_wtime()-t1
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_cdasb(desc_a,info,mold=imold)
    tcdasb = psb_wtime()-t1

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) then
      if (present(amold)) then
        call psb_spasb(a,desc_a,info,mold=amold)
      else
        call psb_spasb(a,desc_a,info,afmt=afmt)
      end if
    end if
    call psb_barrier(ctxt)
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
    call psb_barrier(ctxt)
    ttot = psb_wtime() - t0

    call psb_amx(ctxt,talc)
    call psb_amx(ctxt,tgen)
    call psb_amx(ctxt,tasb)
    call psb_amx(ctxt,ttot)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_s_gen_pde3d
  function outside(i,j,k,bndx,bndy,bndz,iamx,iamy,iamz) result(res)
    logical :: res
    integer(psb_ipk_), intent(in) :: i,j,k,iamx,iamy,iamz
    integer(psb_ipk_), intent(in) :: bndx(0:),bndy(0:),bndz(0:)

    res = (i<bndx(iamx)).or.(i>=bndx(iamx+1)) &
         & .or.(j<bndy(iamy)).or.(j>=bndy(iamy+1)) &
         & .or.(k<bndz(iamz)).or.(k>=bndz(iamz+1))
  end function outside
end module psb_s_pde3d_mod

program psb_s_pde3d
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use psb_s_pde3d_mod
#if defined(OPENMP)
  use omp_lib
#endif
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim
  integer(psb_epk_) :: system_size

  ! miscellaneous
  real(psb_spk_), parameter :: one = sone
  real(psb_dpk_) :: t1, t2, tprec

  ! sparse matrix and preconditioner
  type(psb_sspmat_type) :: a
  type(psb_sprec_type)  :: prec
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense vectors
  type(psb_s_vect_type) :: xxv,bv
  ! parallel environment
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: iam, np, nth

  ! solver parameters
  integer(psb_ipk_) :: iter, itmax,itrace, istopc, irst, ipart
  integer(psb_epk_) :: amatsize, precsize, descsize, d2size
  real(psb_spk_)   :: err, eps

  ! Parameters for solvers in Block-Jacobi preconditioner
  type ainvparms
    character(len=12) :: alg, orth_alg, ilu_alg, ilut_scale
    integer(psb_ipk_) :: fill, inv_fill
    real(psb_spk_)    :: thresh, inv_thresh
  end type ainvparms
  type(ainvparms)     :: parms

  ! other variables
  integer(psb_ipk_) :: info, i
  character(len=20) :: name,ch_err
  character(len=40) :: fname

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#if defined(OPENMP)
  !$OMP parallel shared(nth)
  !$OMP master
  nth = omp_get_num_threads()
  !$OMP end master
  !$OMP end parallel
#else
  nth = 1
#endif
  
  if (iam < 0) then
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif
  if(psb_errstatus_fatal()) goto 9999
  name='pde3d90'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  !
  !  get parameters
  !
  call get_parms(ctxt,kmethd,ptype,afmt,idim,istopc,itmax,itrace,irst,ipart,parms)
  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_gen_pde3d(ctxt,idim,a,bv,xxv,desc_a,afmt,info,partition=ipart)
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_gen_pde3d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  !
  !  prepare the preconditioner.
  !
  if(iam == psb_root_) write(psb_out_unit,'("Setting preconditioner to : ",a)')ptype
  call prec%init(ctxt,ptype,info)
  !
  ! Set the options for the BJAC preconditioner
  !
  if (psb_toupper(ptype) == "BJAC") then
      call prec%set('sub_solve',       parms%alg,   info)
      select case (psb_toupper(parms%alg))
      case ("ILU")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('ilu_alg',         parms%ilu_alg,    info)
      case ("ILUT")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('sub_iluthrs',     parms%thresh,     info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case ("AINV")
        call prec%set('inv_thresh',      parms%inv_thresh, info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
        call prec%set('ainv_alg',        parms%orth_alg,   info)
      case ("INVK")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case ("INVT")
        call prec%set('sub_fillin',      parms%fill,       info)
        call prec%set('inv_fillin',      parms%inv_fill,   info)
        call prec%set('sub_iluthrs',     parms%thresh,     info)
        call prec%set('inv_thresh',      parms%inv_thresh, info)
        call prec%set('ilut_scale',      parms%ilut_scale, info)
      case default
        ! Do nothing, use default setting in the init routine
      end select
  else
    ! nothing to set for NONE or DIAG preconditioner
  end if

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call prec%build(a,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_precbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  tprec = psb_wtime()-t1

  call psb_amx(ctxt,tprec)

  if (iam == psb_root_) write(psb_out_unit,'("Preconditioner time : ",es12.5)')tprec
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  call prec%descr(info)
  !
  ! iterative method parameters
  !
  if(iam == psb_root_) write(psb_out_unit,'("Calling iterative method ",a)')kmethd
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  eps   = 1.d-6
  call psb_krylov(kmethd,a,prec,bv,xxv,eps,desc_a,info,&
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,irst=irst)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='solver routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)
  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  precsize = prec%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  call psb_sum(ctxt,precsize)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Number of processes           : ",i12)')np
    write(psb_out_unit,'("Number of threads             : ",i12)')nth
    write(psb_out_unit,'("Total number of tasks         : ",i12)')nth*np
    write(psb_out_unit,'("Linear system size            : ",i12)') system_size
    write(psb_out_unit,'("Time to solve system          : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration            : ",es12.5)')t2/iter
    write(psb_out_unit,'("Number of iterations          : ",i12)')iter
    write(psb_out_unit,'("Convergence indicator on exit : ",es12.5)')err
    write(psb_out_unit,'("Info  on exit                 : ",i12)')info
    write(psb_out_unit,'("Total memory occupation for      A: ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for   PREC: ",i12)')precsize
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(psb_out_unit,'("Storage format for               A: ",a)') a%get_fmt()
    write(psb_out_unit,'("Storage format for          DESC_A: ",a)') desc_a%get_fmt()
  end if


  !
  !  cleanup storage and exit
  !
  call psb_gefree(bv,desc_a,info)
  call psb_gefree(xxv,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ctxt,kmethd,ptype,afmt,idim,istopc,&
       & itmax,itrace,irst,ipart,parms)
    type(psb_ctxt_type) :: ctxt
    character(len=*) :: kmethd, ptype, afmt
    integer(psb_ipk_) :: idim, istopc,itmax,itrace,irst,ipart
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: ip, inp_unit
    character(len=1024)   :: filename
    type(ainvparms)   :: parms

    call psb_info(ctxt, iam, np)

    if (iam == 0) then
      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ctxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      read(inp_unit,*) ip
      if (ip >= 3) then
        read(inp_unit,*) kmethd
        read(inp_unit,*) ptype
        read(inp_unit,*) afmt

        read(inp_unit,*) idim
        if (ip >= 4) then
          read(inp_unit,*) ipart
        else
          ipart = 3
        endif
        if (ip >= 5) then
          read(inp_unit,*) istopc
        else
          istopc=1
        endif
        if (ip >= 6) then
          read(inp_unit,*) itmax
        else
          itmax=500
        endif
        if (ip >= 7) then
          read(inp_unit,*) itrace
        else
          itrace=-1
        endif
        if (ip >= 8) then
          read(inp_unit,*) irst
        else
          irst=1
        endif
        if (ip >= 9) then
          read(inp_unit,*) parms%alg
          read(inp_unit,*) parms%ilu_alg
          read(inp_unit,*) parms%ilut_scale
          read(inp_unit,*) parms%fill
          read(inp_unit,*) parms%inv_fill
          read(inp_unit,*) parms%thresh
          read(inp_unit,*) parms%inv_thresh
          read(inp_unit,*) parms%orth_alg
        else
          parms%alg =  'ILU' ! Block Solver ILU,ILUT,INVK,AINVT,AORTH
          parms%ilu_alg = 'NONE' ! If ILU : MILU or NONE othewise ignored
          parms%ilut_scale = 'NONE' ! If ILUT: NONE, MAXVAL, DIAG, ARWSUM, ACLSUM, ARCSUM
          parms%fill = 0     ! Level of fill for forward factorization
          parms%inv_fill = 1 ! Level of fill for inverse factorization (only INVK)
          parms%thresh = 1E-1_psb_spk_ ! Threshold for forward factorization
          parms%inv_thresh = 1E-1_psb_spk_ ! Threshold for inverse factorization
          parms%orth_alg = 'LLK'  ! What orthogonalization algorithm?
        endif

        write(psb_out_unit,'("Solving matrix       : ell1")')
        write(psb_out_unit,&
             & '("Grid dimensions      : ",i4," x ",i4," x ",i4)') &
             & idim,idim,idim
        write(psb_out_unit,'("Number of processors : ",i0)')np
        select case(ipart)
        case(1)
          write(psb_out_unit,'("Data distribution    : BLOCK")')
        case(3)
          write(psb_out_unit,'("Data distribution    : 3D")')
        case default
          ipart = 3
          write(psb_out_unit,'("Unknown data distrbution, defaulting to 3D")')
        end select
        write(psb_out_unit,'("Preconditioner       : ",a)') ptype
        if( psb_toupper(ptype) == "BJAC" ) then
          write(psb_out_unit,'("Block subsolver      : ",a)') parms%alg
          select case (psb_toupper(parms%alg))
          case ('ILU')
            write(psb_out_unit,'("Fill in              : ",i0)') parms%fill
            write(psb_out_unit,'("MILU                 : ",a)') parms%ilu_alg
          case ('ILUT')
            write(psb_out_unit,'("Fill in              : ",i0)') parms%fill
            write(psb_out_unit,'("Threshold            : ",es12.5)') parms%thresh
            write(psb_out_unit,'("Scaling              : ",a)') parms%ilut_scale
          case ('INVK')
            write(psb_out_unit,'("Fill in              : ",i0)') parms%fill
            write(psb_out_unit,'("Invese Fill in       : ",i0)') parms%inv_fill
            write(psb_out_unit,'("Scaling              : ",a)') parms%ilut_scale
          case ('INVT')
            write(psb_out_unit,'("Fill in              : ",i0)') parms%fill
            write(psb_out_unit,'("Threshold            : ",es12.5)') parms%thresh
            write(psb_out_unit,'("Invese Fill in       : ",i0)') parms%inv_fill
            write(psb_out_unit,'("Inverse Threshold    : ",es12.5)') parms%inv_thresh
            write(psb_out_unit,'("Scaling              : ",a)') parms%ilut_scale
          case ('AINV','AORTH')
            write(psb_out_unit,'("Inverse Threshold    : ",es12.5)') parms%inv_thresh
            write(psb_out_unit,'("Invese Fill in       : ",i0)') parms%inv_fill
            write(psb_out_unit,'("Orthogonalization    : ",a)') parms%orth_alg
            write(psb_out_unit,'("Scaling              : ",a)') parms%ilut_scale
          case default
              write(psb_out_unit,'("Unknown diagonal solver")')
          end select
        end if
        write(psb_out_unit,'("Iterative method     : ",a)') kmethd
        write(psb_out_unit,'(" ")')
      else
        ! wrong number of parameter, print an error message and exit
        call pr_usage(izero)
        call psb_abort(ctxt)
        stop 1
      endif
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if

    end if
    ! broadcast parameters to all processors
    call psb_bcast(ctxt,kmethd)
    call psb_bcast(ctxt,afmt)
    call psb_bcast(ctxt,ptype)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,ipart)
    call psb_bcast(ctxt,istopc)
    call psb_bcast(ctxt,itmax)
    call psb_bcast(ctxt,itrace)
    call psb_bcast(ctxt,irst)
    call psb_bcast(ctxt,parms%alg)
    call psb_bcast(ctxt,parms%fill)
    call psb_bcast(ctxt,parms%inv_fill)
    call psb_bcast(ctxt,parms%thresh)
    call psb_bcast(ctxt,parms%inv_thresh)
    call psb_bcast(ctxt,parms%orth_alg)
    call psb_bcast(ctxt,parms%ilut_scale)

    return

  end subroutine get_parms
  !
  !  print an error message
  !
  subroutine pr_usage(iout)
    integer(psb_ipk_) :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  pde3d90 methd prec dim &
         &[istop itmax itrace]'
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl'
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     ipart     data partition  1  3      '
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations '
  end subroutine pr_usage

end program psb_s_pde3d
