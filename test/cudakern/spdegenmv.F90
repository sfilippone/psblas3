!    
!                Parallel Sparse BLAS  GPU plugin
!      (C) Copyright 2013
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
! File: spdegenmv.f90
!
! Program: pdegenmv
! This sample program measures the performance of the matrix-vector product.
! The matrix is generated in the same way as for the pdegen test case of
! the main PSBLAS library.
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
  function b1(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none 
    real(psb_spk_) :: b1
    real(psb_spk_), intent(in) :: x,y,z
    b1=sone/sqrt((3*sone))
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none 
    real(psb_spk_) ::  b2
    real(psb_spk_), intent(in) :: x,y,z
    b2=sone/sqrt((3*sone))
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_spk_, sone, szero
    implicit none 
    real(psb_spk_) ::  b3
    real(psb_spk_), intent(in) :: x,y,z      
    b3=sone/sqrt((3*sone))
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
    type(psb_ctxt_type) :: ctxt
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
    integer(psb_ipk_) :: nnz,nr,nlr,i,j,ii,ib,k, partition_
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
    integer(psb_lpk_), allocatable     :: irow(:),icol(:),myidx(:)
    real(psb_spk_), allocatable :: val(:)
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

    deltah   = sone/(idim+2)
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
      call mpi_dims_create(np,3,npdims,info)
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
      
    case default
      write(psb_err_unit,*) iam, 'Initialization error: should not get here'
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end select

    
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz,&
         & dupl=psb_dupl_err_)
    ! define  rhs from boundary conditions; also build initial guess 
    if (info == psb_success_) call psb_geall(xv,desc_a,info)
    if (info == psb_success_) call psb_geall(bv,desc_a,info)

    call psb_barrier(ctxt)
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


    ! loop over rows belonging to current process in a block
    ! distribution.

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1) 
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
      call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
      if(info /= psb_success_) exit
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),bv,desc_a,info)
      if(info /= psb_success_) exit
      zt(:)=szero
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


end module psb_s_pde3d_mod


program pdgenmv
  use psb_base_mod
  use psb_util_mod 
  use psb_ext_mod
#ifdef HAVE_CUDA
  use psb_cuda_mod
#endif
  use psb_s_pde3d_mod
  implicit none

  ! input parameters
  character(len=5)  :: acfmt, agfmt
  integer   :: idim

  ! miscellaneous 
  real(psb_spk_), parameter :: one = 1.e0
  real(psb_dpk_) :: t1, t2, tprec, flops, tflops,&
       & tt1, tt2, gt1, gt2, gflops, bdwdth,&
       & tcnvcsr, tcnvc1, tcnvgpu, tcnvg1

  ! sparse matrix and preconditioner
  type(psb_sspmat_type) :: a, agpu, aux_a
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  type(psb_s_vect_type), target :: xv,bv, xg, bg 
#ifdef HAVE_CUDA
  type(psb_s_vect_cuda)  :: vmold
  type(psb_i_vect_cuda)  :: imold 
#endif
  real(psb_spk_), allocatable :: x1(:), x2(:), x0(:)
  ! blacs parameters
  type(psb_ctxt_type) :: ctxt
  integer            :: iam, np

  ! solver parameters
  integer(psb_epk_) :: amatsize, precsize, descsize, annz, nbytes
  real(psb_spk_)   :: err, eps
  integer, parameter :: ntests=200, ngpu=50, ncnv=20
  type(psb_s_coo_sparse_mat), target   :: acoo
  type(psb_s_csr_sparse_mat), target   :: acsr
  type(psb_s_ell_sparse_mat), target   :: aell
  type(psb_s_hll_sparse_mat), target   :: ahll
  type(psb_s_dia_sparse_mat), target   :: adia
  type(psb_s_hdia_sparse_mat), target   :: ahdia
#ifdef HAVE_CUDA
  type(psb_s_cuda_elg_sparse_mat), target   :: aelg
  type(psb_s_cuda_csrg_sparse_mat), target  :: acsrg
#if CUDA_SHORT_VERSION <= 10
  type(psb_s_cuda_hybg_sparse_mat), target  :: ahybg
#endif
  type(psb_s_cuda_hlg_sparse_mat), target   :: ahlg
  type(psb_s_cuda_hdiag_sparse_mat), target   :: ahdiag
  type(psb_s_cuda_dnsg_sparse_mat), target   :: adnsg
#endif
  class(psb_s_base_sparse_mat), pointer :: agmold, acmold
  ! other variables
  logical, parameter :: dump=.false.
  integer(psb_ipk_)  :: info, i, j, nr, nrg
  integer(psb_lpk_)  :: ig
  character(len=20)  :: name,ch_err
  character(len=40)  :: fname

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

#ifdef HAVE_CUDA
  call psb_cuda_init(ctxt)
#endif
#ifdef HAVE_RSB
  call psb_rsb_init()
#endif

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='pdegenmv-cuda'
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
#ifdef HAVE_CUDA
  write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
  write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())  
#endif
  !
  !  get parameters
  !
  call get_parms(ctxt,acfmt,agfmt,idim)

  !
  !  allocate and fill in the coefficient matrix and initial vectors
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,'CSR  ',info,partition=3)  
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')

  if (dump) then 
    write(fname,'(a,i3.3,a,i3.3,a,i3.3,a)')  'pde',idim,'-',iam,'-',np,'.mtx'
    call a%print(fname,head='PDEGEN test matrix')
  end if

  select case(psb_toupper(acfmt))
  case('ELL')
    acmold => aell
  case('HLL')
    acmold => ahll
  case('DIA')
    acmold => adia
  case('HDIA')
    acmold => ahdia
  case('CSR')
    acmold => acsr
  case('COO')
    acmold => acoo
#ifdef HAVE_RSB
  case('RSB')
    acmold => arsb
#endif
  case default
    write(*,*) 'Unknown format defaulting to HLL'
    acmold => ahll
  end select
  call a%cscnv(info,mold=acmold)
  if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
    write(0,*) 'From cscnv ',info
    call psb_error()
    stop
  end if

#ifdef HAVE_CUDA
  select case(psb_toupper(agfmt))
  case('ELG')
    agmold => aelg
  case('HLG')
    agmold => ahlg
  case('HDIAG')
    agmold => ahdiag
  case('CSRG')
    agmold => acsrg
  case('DNSG')
    agmold => adnsg
#if CUDA_SHORT_VERSION <= 10
  case('HYBG')
    agmold => ahybg
#endif
  case default
    write(*,*) 'Unknown format defaulting to HLG'
    agmold => ahlg
  end select
  call a%cscnv(agpu,info,mold=agmold)
  if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
    write(0,*) 'From cscnv ',info
    call psb_error()
    stop
  end if
  call desc_a%cnv(mold=imold)

  call psb_geasb(bg,desc_a,info,scratch=.true.,mold=vmold)
  call psb_geasb(xg,desc_a,info,scratch=.true.,mold=vmold)
#endif
  nr       = desc_a%get_local_rows()
  nrg      = desc_a%get_global_rows() 
  call psb_geall(x0,desc_a,info)
  do i=1, nr
    call desc_a%l2g(i,ig,info)
    x0(i) = 1.0 + (1.0*ig)/nrg
  end do
  call a%cscnv(aux_a,info,mold=acoo)
  tcnvcsr = 0
  tcnvgpu = 0
  call psb_geall(x1,desc_a,info)
  do j=1, ncnv
    call aux_a%cscnv(a,info,mold=acoo)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call a%cscnv(info,mold=acmold)
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    tcnvcsr = tcnvcsr + t2
    if (j==1) tcnvc1 = t2
    call psb_geasb(x1,desc_a,info)
    call xv%bld(x0)
    call psb_geasb(bv,desc_a,info,scratch=.true.)
    
#ifdef HAVE_CUDA
    
    call aux_a%cscnv(agpu,info,mold=acoo)
    call xg%bld(x0,mold=vmold)
    call psb_geasb(bg,desc_a,info,scratch=.true.,mold=vmold)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call agpu%cscnv(info,mold=agmold)
    call psb_cuda_DeviceSync()
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    if (j==1) tcnvg1 = t2
    tcnvgpu = tcnvgpu + t2
#endif
  end do


  call xv%set(x0)
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(sone,a,xv,szero,bv,desc_a,info)
  end do
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)

#ifdef HAVE_CUDA
  call xg%set(x0)

  ! FIXME: cache flush needed here
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  
  call psb_barrier(ctxt)
  tt1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(sone,agpu,xv,szero,bg,desc_a,info)
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 1 spmm',info,i,ntests
      call psb_error()
      stop
    end if

  end do
  call psb_cuda_DeviceSync()
  call psb_barrier(ctxt)
  tt2 = psb_wtime() - tt1
  call psb_amx(ctxt,tt2)
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  nr       = desc_a%get_local_rows() 
  eps = maxval(abs(x1(1:nr)-x2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on xGPU',eps


  ! FIXME: cache flush needed here
  call xg%set(x0)
  call xg%sync()
  call psb_barrier(ctxt)
  gt1 = psb_wtime()
  do i=1,ntests*ngpu
    call psb_spmm(sone,agpu,xg,szero,bg,desc_a,info)
    ! For timing purposes we need to make sure all threads
    ! in the device are done. 
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 2 spmm',info,i,ntests
      call psb_error()
      stop
    end if
    
  end do
  call psb_cuda_DeviceSync()
  call psb_barrier(ctxt)
  gt2 = psb_wtime() - gt1
  call psb_amx(ctxt,gt2)
  call bg%sync()
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  call psb_geaxpby(-sone,bg,+sone,bv,desc_a,info)
  eps = psb_geamax(bv,desc_a,info)

  call psb_amx(ctxt,t2)
  eps = maxval(abs(x1(1:nr)-x2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on GPU',eps
  if (dump) then 
    write(fname,'(a,i3.3,a,i3.3,a)')'XCPU-out-',iam,'-',np,'.mtx'
    call mm_array_write(x1(1:nr),'Local part CPU',info,filename=fname)
    write(fname,'(a,i3.3,a,i3.3,a)')'XGPU-out-',iam,'-',np,'.mtx'
    call mm_array_write(x2(1:nr),'Local part GPU',info,filename=fname)
  end if
#endif
  annz     = a%get_nzeros()
  amatsize = a%sizeof()
  descsize = psb_sizeof(desc_a)
  call psb_sum(ctxt,nr)
  call psb_sum(ctxt,annz)
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  
  if (iam == psb_root_) then
    write(psb_out_unit,&
         & '("Matrix: ell1 ",i0)') idim
    write(psb_out_unit,&
         &'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,&
         &'("Size of matrix                   : ",i20,"           ")') nr
    write(psb_out_unit,&
         &'("Number of nonzeros               : ",i20,"           ")') annz
    write(psb_out_unit,&
         &'("Memory occupation                : ",i20,"           ")') amatsize
    flops  = ntests*(2.d0*annz)
    tflops = flops
    gflops = flops * ngpu
    write(psb_out_unit,'("Storage type for    A: ",a)') a%get_fmt()
#ifdef HAVE_CUDA
    write(psb_out_unit,'("Storage type for AGPU: ",a)') agpu%get_fmt()
    write(psb_out_unit,'("Time to convert A from COO to CPU (1): ",F20.9)')&
         & tcnvc1
    write(psb_out_unit,'("Time to convert A from COO to CPU (t): ",F20.9)')&
         & tcnvcsr
    write(psb_out_unit,'("Time to convert A from COO to CPU (a): ",F20.9)')&
         & tcnvcsr/ncnv
    write(psb_out_unit,'("Time to convert A from COO to GPU (1): ",F20.9)')&
         & tcnvg1
    write(psb_out_unit,'("Time to convert A from COO to GPU (t): ",F20.9)')&
         & tcnvgpu
    write(psb_out_unit,'("Time to convert A from COO to GPU (a): ",F20.9)')&
         & tcnvgpu/ncnv

#endif
    write(psb_out_unit,&
         & '("Number of flops (",i0," prod)        : ",F20.0,"           ")') &
         &  ntests,flops

    flops  = flops / (t2)
    tflops = tflops / (tt2)
    gflops = gflops / (gt2)

    write(psb_out_unit,'("Time for ",i6," products (s) (CPU)   : ",F20.3)')&
         &  ntests,t2
    write(psb_out_unit,'("Time per product    (ms)     (CPU)   : ",F20.3)')&
         & t2*1.d3/(1.d0*ntests)
    write(psb_out_unit,'("MFLOPS                       (CPU)   : ",F20.3)')&
         & flops/1.d6
#ifdef HAVE_CUDA
    write(psb_out_unit,'("Time for ",i6," products (s) (xGPU)  : ",F20.3)')&
         & ntests, tt2
    write(psb_out_unit,'("Time per product    (ms)     (xGPU)  : ",F20.3)')&
         & tt2*1.d3/(1.d0*ntests)
    write(psb_out_unit,'("MFLOPS                       (xGPU)  : ",F20.3)')&
         & tflops/1.d6

    write(psb_out_unit,'("Time for ",i6," products (s) (GPU.)  : ",F20.3)')&
         & ngpu*ntests,gt2
    write(psb_out_unit,'("Time per product    (ms)     (GPU.)  : ",F20.3)')&
         & gt2*1.d3/(1.d0*ntests*ngpu)
    write(psb_out_unit,'("MFLOPS                       (GPU.)  : ",F20.3)')&
         & gflops/1.d6
#endif
    !
    ! This computation assumes the data movement associated with CSR:
    ! it is minimal in terms of coefficients. Other formats may either move
    ! more data (padding etc.) or less data (if they can save on the indices). 
    !
    nbytes = nr*(2*psb_sizeof_sp + psb_sizeof_ip)+&
         & annz*(psb_sizeof_sp + psb_sizeof_ip)
    bdwdth = ntests*nbytes/(t2*1.d6)
    write(psb_out_unit,*)
    write(psb_out_unit,'("MBYTES/S sust. effective bandwidth  (CPU)  : ",F20.3)') bdwdth
#ifdef HAVE_CUDA
    bdwdth = ngpu*ntests*nbytes/(gt2*1.d6)
    write(psb_out_unit,'("MBYTES/S sust. effective bandwidth  (GPU)  : ",F20.3)') bdwdth
    bdwdth = psb_cuda_MemoryPeakBandwidth()
    write(psb_out_unit,'("MBYTES/S peak bandwidth             (GPU)  : ",F20.3)') bdwdth
#endif
    write(psb_out_unit,'("Storage type for DESC_A: ",a)') desc_a%indxmap%get_fmt()
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize

  end if

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(bv,desc_a,info)
  call psb_gefree(xv,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
#ifdef HAVE_CUDA
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ctxt,acfmt,agfmt,idim)
    type(psb_ctxt_type) :: ctxt
    character(len=*) :: agfmt, acfmt
    integer      :: idim
    integer      :: np, iam
    integer      :: intbuf(10), ip

    call psb_info(ctxt, iam, np)

    if (iam == 0) then
      write(*,*) 'CPU side format?'
      read(psb_inp_unit,*) acfmt
      write(*,*) 'CUDA side format?'
      read(psb_inp_unit,*) agfmt
      write(*,*) 'Size of discretization cube?'
      read(psb_inp_unit,*) idim
    endif
    call psb_bcast(ctxt,acfmt)
    call psb_bcast(ctxt,agfmt)
    call psb_bcast(ctxt,idim)

    if (iam == 0) then
      write(psb_out_unit,'("Testing matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)')np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'(" ")')
    end if
    return

  end subroutine get_parms


end program pdgenmv
