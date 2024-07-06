module psb_d_pde3d_mod

   use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_lpk_, psb_desc_type,&
   &  psb_dspmat_type, psb_d_multivect_type, dzero,&
   &  psb_d_base_sparse_mat, psb_d_base_multivect_type, &
   &  psb_i_base_vect_type, psb_l_base_vect_type

   interface
      function d_func_3d(x,y,z) result(val)
         import :: psb_dpk_
         real(psb_dpk_), intent(in) :: x,y,z
         real(psb_dpk_) :: val
      end function d_func_3d
   end interface

   interface psb_gen_pde3d
      module procedure  psb_d_gen_pde3d
   end interface psb_gen_pde3d

contains

   function d_null_func_3d(x,y,z) result(val)

      real(psb_dpk_), intent(in) :: x,y,z
      real(psb_dpk_) :: val

      val = dzero

   end function d_null_func_3d
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
   ! change the B1/B2/B3 functions to e.g. done/sqrt((3*done))
   !
   function b1(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) :: b1
      real(psb_dpk_), intent(in) :: x,y,z
      b1=done/sqrt((3*done))
   end function b1
   function b2(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  b2
      real(psb_dpk_), intent(in) :: x,y,z
      b2=done/sqrt((3*done))
   end function b2
   function b3(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  b3
      real(psb_dpk_), intent(in) :: x,y,z
      b3=done/sqrt((3*done))
   end function b3
   function c(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  c
      real(psb_dpk_), intent(in) :: x,y,z
      c=dzero
   end function c
   function a1(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  a1
      real(psb_dpk_), intent(in) :: x,y,z
      a1=done/80
   end function a1
   function a2(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  a2
      real(psb_dpk_), intent(in) :: x,y,z
      a2=done/80
   end function a2
   function a3(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  a3
      real(psb_dpk_), intent(in) :: x,y,z
      a3=done/80
   end function a3
   function g(x,y,z)
      use psb_base_mod, only : psb_dpk_, done, dzero
      implicit none
      real(psb_dpk_) ::  g
      real(psb_dpk_), intent(in) :: x,y,z
      g = dzero
      if (x == done) then
         g = done
      else if (x == dzero) then
         g = exp(y**2-z**2)
      end if
   end function g


   !
   !  subroutine to allocate and fill in the coefficient matrix and
   !  the rhs.
   !
   subroutine psb_d_gen_pde3d(ctxt,idim,a,bmv,xmv,nrhs,desc_a,afmt,info,&
   & f,amold,vmold,imold,partition,nrl,iv,tnd)
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
      type(psb_dspmat_type) :: a
      type(psb_d_multivect_type) :: xmv,bmv
      integer(psb_ipk_)     :: nrhs
      type(psb_desc_type)   :: desc_a
      type(psb_ctxt_type) :: ctxt
      integer(psb_ipk_)     :: info
      character(len=*)      :: afmt
      procedure(d_func_3d), optional :: f
      class(psb_d_base_sparse_mat), optional :: amold
      class(psb_d_base_multivect_type), optional :: vmold
      class(psb_i_base_vect_type), optional :: imold
      integer(psb_ipk_), optional :: partition, nrl,iv(:)
      logical, optional :: tnd
      ! Local variables.

      integer(psb_ipk_), parameter :: nb=20
      type(psb_d_csc_sparse_mat)  :: acsc
      type(psb_d_coo_sparse_mat)  :: acoo
      type(psb_d_csr_sparse_mat)  :: acsr
      real(psb_dpk_)           :: zt(nb,nrhs),x,y,z
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
      real(psb_dpk_), allocatable :: val(:)
      ! deltah dimension of each grid cell
      ! deltat discretization time
      real(psb_dpk_)            :: deltah, sqdeltah, deltah2
      real(psb_dpk_), parameter :: rhs=dzero,one=done,zero=dzero
      real(psb_dpk_)    :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
      integer(psb_ipk_) :: err_act
      procedure(d_func_3d), pointer :: f_
      logical :: tnd_
      character(len=20)  :: name, ch_err,tmpfmt

      info = psb_success_
      name = 'create_matrix'
      call psb_erractionsave(err_act)

      call psb_info(ctxt, iam, np)


      if (present(f)) then
         f_ => f
      else
         f_ => d_null_func_3d
      end if

      deltah   = done/(idim+2)
      sqdeltah = deltah*deltah
      deltah2  = (2*done)* deltah

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
      if (info == psb_success_) call psb_geall(xmv,desc_a,info,n=nrhs)
      if (info == psb_success_) call psb_geall(bmv,desc_a,info,n=nrhs)

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
            zt(k,:) = f_(x,y,z)
            ! internal point: build discretization
            !
            !  term depending on   (x-1,y,z)
            !
            val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
            if (ix == 1) then
               zt(k,:) = g(dzero,y,z)*(-val(icoeff)) + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix-1,iy,iz,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif
            !  term depending on     (x,y-1,z)
            val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
            if (iy == 1) then
               zt(k,:) = g(x,dzero,z)*(-val(icoeff))   + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix,iy-1,iz,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif
            !  term depending on     (x,y,z-1)
            val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
            if (iz == 1) then
               zt(k,:) = g(x,y,dzero)*(-val(icoeff))   + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix,iy,iz-1,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif

            !  term depending on     (x,y,z)
            val(icoeff)=(2*done)*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah &
            & + c(x,y,z)
            call ijk2idx(icol(icoeff),ix,iy,iz,idim,idim,idim)
            irow(icoeff) = glob_row
            icoeff       = icoeff+1
            !  term depending on     (x,y,z+1)
            val(icoeff)=-a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
            if (iz == idim) then
               zt(k,:) = g(x,y,done)*(-val(icoeff))   + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix,iy,iz+1,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif
            !  term depending on     (x,y+1,z)
            val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
            if (iy == idim) then
               zt(k,:) = g(x,done,z)*(-val(icoeff))   + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix,iy+1,iz,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif
            !  term depending on     (x+1,y,z)
            val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
            if (ix==idim) then
               zt(k,:) = g(done,y,z)*(-val(icoeff))   + zt(k,:)
            else
               call ijk2idx(icol(icoeff),ix+1,iy,iz,idim,idim,idim)
               irow(icoeff) = glob_row
               icoeff       = icoeff+1
            endif

         end do
         call psb_spins(icoeff-1,irow,icol,val,a,desc_a,info)
         if(info /= psb_success_) exit
         call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,:),bmv,desc_a,info)
         if(info /= psb_success_) exit
         zt(:,:)=dzero
         call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,:),xmv,desc_a,info)
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
            call psb_spasb(a,desc_a,info,mold=amold,bld_and=tnd)
         else
            call psb_spasb(a,desc_a,info,afmt=afmt,bld_and=tnd)
         end if
      end if
      call psb_barrier(ctxt)
      if(info /= psb_success_) then
         info=psb_err_from_subroutine_
         ch_err='asb rout.'
         call psb_errpush(info,name,a_err=ch_err)
         goto 9999
      end if
      if (info == psb_success_) call psb_geasb(xmv,desc_a,info,mold=vmold)
      if (info == psb_success_) call psb_geasb(bmv,desc_a,info,mold=vmold)
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

9999  call psb_error_handler(ctxt,err_act)

      return
   end subroutine psb_d_gen_pde3d

end module psb_d_pde3d_mod

program dpdegen
   use psb_base_mod
   use psb_util_mod
   use psb_prec_mod
   use psb_krylov_mod
   use psb_ext_mod
   use psb_cuda_mod
   use psb_d_pde3d_mod

   ! input parameters
   character(len=40) :: kmethd, ptype
   character(len=5)  :: afmt
   integer(psb_ipk_) :: idim, nrhs, istopc, itmax, itrace
   real(psb_dpk_)    :: eps

   ! sparse matrix
   type(psb_dspmat_type) :: a

   ! preconditioner data
   type(psb_dprec_type) :: prec

   ! dense matrices
   real(psb_dpk_), allocatable :: b_mv_glob(:,:), x_mv_glob(:,:), r_mv_glob(:,:)
   type(psb_d_multivect_type)  :: b_mv, x_mv, r_mv
   type(psb_d_multivect_cuda)  :: gpumold
   type(psb_i_vect_cuda)       :: imold
   integer(psb_ipk_)           :: m

   ! blacs parameters
   type(psb_desc_type) :: desc_a
   type(psb_ctxt_type) :: ctxt
   integer             :: iam, np

   ! solver parameters
   real(psb_dpk_)     :: err, cond
   integer(psb_epk_)  :: amatsize, precsize, descsize, annz
   integer(psb_ipk_)  :: iter, ierr, ircode

   ! molds
   type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
   type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
   type(psb_d_cuda_elg_sparse_mat), target   :: aelg
   type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
   class(psb_d_base_sparse_mat), pointer     :: agmold

   ! other variables
   integer(psb_ipk_)           :: info, i, j, rep
   real(psb_dpk_)              :: t1, t2, tprec
   character(len=20)           :: name, ch_err
   real(psb_dpk_)              :: random_value
   real(psb_dpk_)              :: resmxp
   real(psb_dpk_), allocatable :: resmx(:)
   logical                     :: tnd = .false.
   logical                     :: print_matrix = .false.

   ! Init environment
   info=psb_success_
   call psb_init(ctxt)
   call psb_info(ctxt,iam,np)
   call psb_cuda_init(ctxt)
   if (iam < 0) then
      ! This should not happen, but just in case
      call psb_exit(ctxt)
      stop
   endif
   if(psb_get_errstatus() /= 0) goto 9999
   name='pdegenmm_cuda'
   !
   ! Hello world
   !
   if (iam == psb_root_) then
      write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
      write(*,*) 'This is the ',trim(name),' sample program'
      write(psb_out_unit,'("Number of processors: ",i8)')np
   end if
   write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice()+1,' out of', psb_cuda_getDeviceCount()
   write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice()+1,' is a: ', trim(psb_cuda_DeviceName())
   !
   !  get parameters
   !
   call get_parms(ctxt,kmethd,ptype,idim,afmt,nrhs,istopc,itmax,itrace,eps)
   !
   !  allocate and fill in the coefficient matrix and initial vectors
   !
   call psb_barrier(ctxt)
   t1 = psb_wtime()
   call psb_gen_pde3d(ctxt,idim,a,b_mv,x_mv,nrhs,desc_a,'CSR  ',info,vmold=gpumold,partition=3,tnd=tnd)
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

   ! building the preconditioner
   call prec%init(ctxt,ptype,info)
   call prec%build(a,desc_a,info)
   if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
   end if

   select case(psb_toupper(afmt))
    case('ELG')
      agmold => aelg
    case('HLG')
      agmold => ahlg
    case('HDIAG')
      agmold => ahdiag
    case('CSRG')
      agmold => acsrg
    case default
      write(*,*) 'Unknown format defaulting to HLG'
      agmold => ahlg
   end select
   call a%cscnv(info,mold=agmold)
   call desc_a%cnv(mold=imold)
   if ((info /= 0).or.(psb_get_errstatus()/=0)) then
      write(0,*) 'From cscnv ',info
      call psb_error()
      stop
   end if

   ! set random RHS
   call b_mv%zero()
   call random_number(b_mv%v%v(1:desc_a%get_local_rows(),:))
   b_mv%v%v(1:desc_a%get_local_rows(),:) = -10 + (20)*b_mv%v%v
   call b_mv%v%set_host()
   call b_mv%sync()

   call psb_geall(r_mv,desc_a,info,nrhs)
   call r_mv%zero()
   call psb_geasb(r_mv,desc_a,info)

   t2 = psb_wtime() - t1
   call psb_amx(ctxt, t2)

   if (iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
      write(psb_out_unit,'(" ")')
   end if

   ! building the preconditioner
   call prec%init(ctxt,ptype,info)
   t1 = psb_wtime()
   call prec%build(a,desc_a,info)
   tprec = psb_wtime()-t1
   if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_precbld')
      goto 9999
   end if

   call psb_amx(ctxt,tprec)

   if(iam == psb_root_) then
      write(psb_out_unit,'("Preconditioner time: ",es12.5)')tprec
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Starting algorithm")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_barrier(ctxt)
   t1 = psb_wtime()

   call psb_krylov(kmethd,a,prec,b_mv,x_mv,eps,desc_a,info,&
   & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc)

   call psb_barrier(ctxt)
   t2 = psb_wtime() - t1
   call psb_amx(ctxt,t2)

   if(iam == psb_root_) then
      write(psb_out_unit,'("Finished algorithm")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_geaxpby(done,b_mv,dzero,r_mv,desc_a,info)
   call psb_spmm(-done,a,x_mv,done,r_mv,desc_a,info)

   resmx  = psb_genrm2(r_mv,desc_a,info)
   resmxp = psb_geamax(r_mv,desc_a,info)

   amatsize = a%sizeof()
   annz     = a%get_nzeros()
   descsize = desc_a%sizeof()
   precsize = prec%sizeof()

   call psb_sum(ctxt,amatsize)
   call psb_sum(ctxt,annz)
   call psb_sum(ctxt,descsize)
   call psb_sum(ctxt,precsize)

   call psb_gather(x_mv_glob,x_mv,desc_a,info,root=psb_root_)
   if (info == psb_success_) call psb_gather(b_mv_glob,b_mv,desc_a,info,root=psb_root_)
   if (info == psb_success_) call psb_gather(r_mv_glob,r_mv,desc_a,info,root=psb_root_)
   if (info /= psb_success_) goto 9999

   if (iam == psb_root_) then
      call prec%descr(info)
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Computed solution on:         ",i8," processors")')np
      write(psb_out_unit,'("Dimesion of A:                      ",i12)')desc_a%get_global_rows()
      write(psb_out_unit,'("Number of nonzeros:                 ",i12)')annz
      write(psb_out_unit,'("Number of RHS:                      ",i12)')nrhs
      write(psb_out_unit,'("Storage format for A:                ",a)')a%get_fmt()
      write(psb_out_unit,'("Storage format for DESC_A:           ",a)')desc_a%get_fmt()
      write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
      write(psb_out_unit,'("Total memory occupation for PREC:   ",i12)')precsize
      write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
      write(psb_out_unit,'("Iterations to convergence:          ",i12)')iter
      write(psb_out_unit,'("Error estimate on exit:             ",es12.5)')err
      write(psb_out_unit,'("Time to buil prec.:                 ",es12.5)')tprec
      write(psb_out_unit,'("Time to solve system:               ",es12.5)')t2
      write(psb_out_unit,'("Time per iteration:                 ",es12.5)')t2/(iter)
      write(psb_out_unit,'("Total time:                         ",es12.5)')t2+tprec
      write(psb_out_unit,'("Residual norm 2:                    ",es12.5)')sum(resmx)/size(resmx)
      write(psb_out_unit,'("Residual norm inf:                  ",es12.5)')resmxp
      if (print_matrix) then
         write(psb_out_unit,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
         do i=1,m
            write(psb_out_unit,993) i, x_mv_glob(i,:), r_mv_glob(i,:), b_mv_glob(i,:)
         end do
      end if
   end if

998 format(i8,20(2x,g20.14))
993 format(i6,20(1x,e12.6))

   !
   !  cleanup storage and exit
   !
   call psb_gefree(b_mv,desc_a,info)
   call psb_gefree(x_mv,desc_a,info)
   call psb_spfree(a,desc_a,info)
   call prec%free(info)
   call psb_cdfree(desc_a,info)
   call psb_cuda_exit()
   call psb_exit(ctxt)

   return

9999 call psb_error(ctxt)

   return

contains
   !
   ! get iteration parameters from standard input
   !
   subroutine get_parms(ctxt,kmethd,ptype,idim,afmt,nrhs,istopc,itmax,itrace,eps)
      type(psb_ctxt_type) :: ctxt
      character(len=40)   :: kmethd, ptype
      character(len=5)    :: afmt
      integer(psb_ipk_)   :: idim, nrhs, istopc, itmax, itrace
      real(psb_dpk_)      :: eps

      integer(psb_ipk_)   :: np, iam
      integer(psb_ipk_)   :: ip, inp_unit
      character(len=1024) :: filename

      call psb_info(ctxt, iam, np)

      if (iam == 0) then
         if (command_argument_count()>0) then
            call get_command_argument(1, filename)
            inp_unit = 30
            open(inp_unit, file=filename, action='read', iostat=info)
            if (info /= 0) then
               write(psb_err_unit,*) 'Could not open file ',filename,' for input'
               call psb_abort(ctxt)
               stop
            else
               write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
            end if
         else
            inp_unit = psb_inp_unit
         end if

         ! Read Input Parameters
         read(inp_unit,*) kmethd
         read(inp_unit,*) ptype
         read(inp_unit,*) idim
         read(inp_unit,*) afmt
         read(inp_unit,*) nrhs
         read(inp_unit,*) istopc
         read(inp_unit,*) itmax
         read(inp_unit,*) itrace
         read(inp_unit,*) eps

         call psb_bcast(ctxt,kmethd)
         call psb_bcast(ctxt,ptype)
         call psb_bcast(ctxt,idim)
         call psb_bcast(ctxt,afmt)
         call psb_bcast(ctxt,nrhs)
         call psb_bcast(ctxt,istopc)
         call psb_bcast(ctxt,itmax)
         call psb_bcast(ctxt,itrace)
         call psb_bcast(ctxt,eps)

         if (inp_unit /= psb_inp_unit) then
            close(inp_unit)
         end if
      else
         ! Receive Parameters
         call psb_bcast(ctxt,kmethd)
         call psb_bcast(ctxt,ptype)
         call psb_bcast(ctxt,idim)
         call psb_bcast(ctxt,afmt)
         call psb_bcast(ctxt,nrhs)
         call psb_bcast(ctxt,istopc)
         call psb_bcast(ctxt,itmax)
         call psb_bcast(ctxt,itrace)
         call psb_bcast(ctxt,eps)
      end if

      if (iam == 0) then
         write(psb_out_unit,'(" ")')
         write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)') idim,idim,idim
         write(psb_out_unit,'("Iterative method     : ",a)') kmethd
         write(psb_out_unit,'("Number of processors : ",i0)') np
         write(psb_out_unit,'("Number of RHS        : ",i4)') nrhs
         write(psb_out_unit,'("Number of iterations : ",i3)') itmax
         write(psb_out_unit,'("Storage format       : ",a)') afmt
         write(psb_out_unit,'(" ")')
      end if

      return

   end subroutine get_parms

end program dpdegen
