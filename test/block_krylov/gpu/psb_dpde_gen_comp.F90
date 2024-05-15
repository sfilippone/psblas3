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
   character(len=40) :: kmethd 	  = "GMRES"
   character(len=40) :: ptype     = "NONE"
   character(len=5)  :: agfmt     = "CSRG"
   integer(psb_ipk_) :: nrhs      = 5
   integer(psb_ipk_) :: istopbg   = 1
   integer(psb_ipk_) :: istoprg   = 2
   integer(psb_ipk_) :: itmax     = 500
   integer(psb_ipk_) :: itrace	  = -1
   integer(psb_ipk_) :: itrs	  = 100
   real(psb_dpk_)    :: eps	      = 1.d-7
   integer(psb_ipk_) :: idim      = 20
   logical           :: tnd       = .false.


   ! sparse matrix
   type(psb_dspmat_type) :: a, aux_a

   ! preconditioner data
   type(psb_dprec_type) :: prec

   ! miscellaneous
   real(psb_dpk_)               :: tb1, tb2, tr1, tr2
   real(psb_dpk_), allocatable  :: tb(:), tr(:)

   ! descriptor
   type(psb_desc_type)   :: desc_a

   ! dense matrices
   type(psb_d_multivect_type), target :: x_mv, b_mv
   type(psb_d_vect_type), target      :: x_col, b_col
   type(psb_d_multivect_cuda)  :: gpumold_mv
   type(psb_d_vect_cuda)       :: gpumold_col
   type(psb_i_vect_cuda)       :: imold
   real(psb_dpk_), allocatable :: b_mv_glob(:,:)
   real(psb_dpk_), allocatable :: b_col_glob(:)

   ! blacs parameters
   type(psb_ctxt_type) :: ctxt
   integer             :: iam, np

   ! solver parameters
   real(psb_dpk_)     :: err, cond
   integer(psb_ipk_)  :: reps = 2

   ! molds
   type(psb_d_cuda_elg_sparse_mat), target   :: aelg
   type(psb_d_cuda_csrg_sparse_mat), target  :: acsrg
   type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
   class(psb_d_base_sparse_mat), pointer     :: agmold

   ! other variables
   integer(psb_ipk_)  :: info, i, j, m_problem, iter, rep
   character(len=20)  :: name, ch_err
   real(psb_dpk_)     :: random_value

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
   name='pdegenmm-cuda'
   !
   ! Hello world
   !
   if (iam == psb_root_) then
      write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
      write(*,*) 'This is the ',trim(name),' sample program'
      write(psb_out_unit,'("Number of processors: ",i8)')np
   end if
   write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
   write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())

   !  allocate and fill in the coefficient matrix and initial vectors
   !
   call psb_barrier(ctxt)
   t1 = psb_wtime()
   call psb_gen_pde3d(ctxt,idim,a,b_mv,x_mv,nrhs,desc_a,'CSR  ',info,partition=3,tnd=tnd)
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

   select case(psb_toupper(agfmt))
    case('ELG')
      agmold => aelg
    case('HLG')
      agmold => ahlg
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

   ! Set RHS
   call psb_geall(b_mv_glob,desc_a,info,n=nrhs)
   do i=1,b_mv%get_nrows()
      do j=1,b_mv%get_ncols()
         call random_number(random_value)
         b_mv_glob(i,j) = random_value
      end do
   end do

   if (iam == psb_root_) then
      allocate(tb(reps),tr(reps))
      tb = dzero
      tr = dzero
   end if

   do rep=1,reps
      call psb_scatter(b_mv_glob,b_mv,desc_a,info,root=psb_root_,mold=gpumold_mv)
      call psb_geall(x_mv,desc_a,info,nrhs)
      call x_mv%zero()
      call psb_geasb(x_mv,desc_a,info,mold=gpumold_mv)

      call psb_barrier(ctxt)
      tb1 = psb_wtime()

      call psb_krylov(kmethd,a,prec,b_mv,x_mv,eps,desc_a,info,&
      & itmax=itmax,iter=iter,err=err,itrace=itrace,&
      & itrs=itrs,istop=istopbg)

      call psb_barrier(ctxt)
      tb2 = psb_wtime() - tb1
      call psb_amx(ctxt,tb2)

      if (iam == psb_root_) then
         tb(rep) = tb2
         write(*,*) 'Time', rep, tb(rep), iter
      end if
      call psb_barrier(ctxt)

      call psb_gefree(b_mv,desc_a,info)
      call psb_gefree(x_mv,desc_a,info)
   end do

   if(iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Finished BGMRES")')
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Starting sGMRES")')
      write(psb_out_unit,'(" ")')
   end if

   call psb_barrier(ctxt)

   do rep=1,reps
      do i=1,nrhs
         b_col_glob = b_mv_glob(:,i)
         call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_,mold=gpumold_col)
         call psb_geall(x_col,desc_a,info)
         call x_col%zero()
         call psb_geasb(x_col,desc_a,info,mold=gpumold_col)

         cond = dzero

         call psb_barrier(ctxt)
         tr1 = psb_wtime()

         call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,&
         & itmax=itmax,iter=iter,err=err,itrace=itrace,&
         & istop=istoprg,irst=itrs,cond=cond)

         call psb_barrier(ctxt)
         tr2 = psb_wtime() - tr1
         call psb_amx(ctxt,tr2)

         if (iam == psb_root_) then
            tr(rep) = tr(rep) + tr2
         end if
         call psb_barrier(ctxt)

         call psb_gefree(b_col, desc_a,info)
         call psb_gefree(x_col, desc_a,info)
      end do

      if (iam == psb_root_) then
         write(*,*) 'Time', rep, tr(rep), iter
      end if
   end do

   call psb_barrier(ctxt)

   if(iam == psb_root_) then
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Finished sGMRES")')
      write(psb_out_unit,'(" ")')
   end if

   !
   !  cleanup storage and exit
   !
   call psb_spfree(a,desc_a,info)
   call psb_cdfree(desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='free routine'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
   end if

   call psb_cuda_exit()
   call psb_exit(ctxt)
   return

9999 continue

   call psb_error(ctxt)

end program dpdegen
