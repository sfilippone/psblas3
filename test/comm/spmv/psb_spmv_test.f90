!> Test program for overlapping communication and computation with psb_spmm.
!!
!!
module psb_spmv_overlap_test

  use psb_base_mod
  use psb_util_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif

  implicit none

  interface 
    function d_func_3d(x,y,z) result(val)
      import :: psb_dpk_
      real(psb_dpk_), intent(in) :: x,y,z
      real(psb_dpk_) :: val
    end function d_func_3d
  end interface 

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
    b1=dzero
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=dzero
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_, done, dzero
    implicit none 
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=dzero
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
  subroutine psb_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info,&
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
    type(psb_d_vect_type) :: xv,bv
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    character(len=*)      :: afmt
    procedure(d_func_3d), optional :: f
    class(psb_d_base_sparse_mat), optional :: amold
    class(psb_d_base_vect_type), optional :: vmold 
    class(psb_i_base_vect_type), optional :: imold
    integer(psb_ipk_), optional :: partition, nrl,iv(:)
    logical, optional :: tnd
    ! Local variables.

    integer(psb_ipk_), parameter :: nb=20
    type(psb_d_csc_sparse_mat)  :: acsc
    type(psb_d_coo_sparse_mat)  :: acoo
    type(psb_d_csr_sparse_mat)  :: acsr
    real(psb_dpk_)           :: zt(nb),x,y,z
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
    integer(psb_ipk_) :: np, iam, nth
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

    deltah   = done/(idim+1)
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

    tnd_ = .false.
    if (present(tnd)) tnd_ = tnd
    
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

    
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
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
          zt(k) = g(dzero,y,z)*(-val(icoeff)) + zt(k)
        else
          call ijk2idx(icol(icoeff),ix-1,iy,iz,idim,idim,idim)
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y-1,z)
        val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
        if (iy == 1) then 
          zt(k) = g(x,dzero,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy-1,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y,z-1)
        val(icoeff)=-a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
        if (iz == 1) then 
          zt(k) = g(x,y,dzero)*(-val(icoeff))   + zt(k)
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
          zt(k) = g(x,y,done)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy,iz+1,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x,y+1,z)
        val(icoeff)=-a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
        if (iy == idim) then 
          zt(k) = g(x,done,z)*(-val(icoeff))   + zt(k)
        else
          call ijk2idx(icol(icoeff),ix,iy+1,iz,idim,idim,idim)          
          irow(icoeff) = glob_row
          icoeff       = icoeff+1
        endif
        !  term depending on     (x+1,y,z)
        val(icoeff)=-a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
        if (ix==idim) then 
          zt(k) = g(done,y,z)*(-val(icoeff))   + zt(k)
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
      zt(:)=dzero
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
    if (present(imold)) then
      call psb_cdasb(desc_a,info,mold=imold)
    else
      call psb_cdasb(desc_a,info)
    end if
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
    if (info == psb_success_) then
      if (present(vmold)) then
        call psb_geasb(xv,desc_a,info,mold=vmold)
        if (info == psb_success_) call psb_geasb(bv,desc_a,info,mold=vmold)
      else
        call psb_geasb(xv,desc_a,info)
        if (info == psb_success_) call psb_geasb(bv,desc_a,info)
      end if
    end if
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
  end subroutine psb_d_gen_pde3d

  subroutine run_spmv_kernel(ctxt,use_gpu,matrix_file,matrix_fmt,cpu_fmt,gpu_fmt,idim_in,times_in,comm_mode)
    use psb_base_mod
#ifdef PSB_HAVE_CUDA
    use psb_cuda_mod
#endif
    implicit none

    type(psb_ctxt_type), intent(in) :: ctxt
    logical, intent(in)             :: use_gpu
    character(len=*), intent(in)    :: matrix_file
    character(len=*), intent(in)    :: matrix_fmt
    character(len=*), intent(in)    :: cpu_fmt
    character(len=*), intent(in)    :: gpu_fmt
    integer(psb_ipk_), intent(in)   :: idim_in, times_in
    character(len=*), intent(in)    :: comm_mode

    type(psb_dspmat_type)           :: a
    type(psb_d_vect_type)           :: x, y
    type(psb_desc_type)             :: desc_a
    character(len=8)                :: afmt
    character(len=64)               :: env_buf
    integer(psb_ipk_)               :: my_rank, np, info, err_act
    integer(psb_ipk_)               :: idim, times, i
    integer                         :: env_len, env_status, ios
    real(psb_dpk_)                  :: alpha, beta, t0, t1, dt, avg_t
    logical                         :: use_external_matrix

    integer(psb_ipk_)               :: comm_type

#ifdef PSB_HAVE_CUDA
    type(psb_d_vect_cuda)                     :: cuda_vector_mold
    type(psb_i_vect_cuda)                     :: cuda_index_mold
    type(psb_d_cuda_elg_sparse_mat), target   :: cuda_ell_sparse_mold
    type(psb_d_cuda_csrg_sparse_mat), target  :: cuda_csr_sparse_mold
    type(psb_d_cuda_hdiag_sparse_mat), target :: cuda_hdia_sparse_mold
    type(psb_d_cuda_hlg_sparse_mat), target   :: cuda_hll_sparse_mold
    class(psb_d_base_sparse_mat), pointer     :: cuda_sparse_mold
#endif

    select case(psb_toupper(trim(comm_mode)))
    case('P2P','ISEND_IRECV')
      comm_type = psb_comm_isend_irecv_
    case('NEIGHBOR','INEIGHBOR_ALLTOALLV')
      comm_type = psb_comm_ineighbor_alltoallv_
    case('PNEIGHBOR','PERSISTENT','PERSISTENT_INEIGHBOR_A2AV')
      comm_type = psb_comm_persistent_ineighbor_alltoallv_
    case('MPI_GET','RMA_PULL')
      comm_type = psb_comm_rma_pull_
    case('MPI_PUT','RMA_PUSH')
      comm_type = psb_comm_rma_push_
    case default
      comm_type = psb_comm_isend_irecv_
      if (my_rank == psb_root_) then
        write(psb_err_unit,'("Unknown comm backend: ",a,", defaulting to P2P")') trim(comm_mode)
      end if
    end select

    info = psb_success_
    afmt = psb_toupper(trim(cpu_fmt))
    if (len_trim(afmt) == 0) afmt = 'CSR'
    if (idim_in > 0) then
      idim = idim_in
    else
      idim = 10
    end if

    if (times_in > 0) then
      times = times_in
    else
      times = 100
    end if
    alpha = done
    beta  = dzero

    call psb_erractionsave(err_act)
    call psb_info(ctxt, my_rank, np)
    use_external_matrix = (len_trim(matrix_file) > 0)

    if (idim_in <= 0) then
      call get_environment_variable('IDIM', env_buf, length=env_len, status=env_status)
      if ((env_status == 0) .and. (env_len > 0)) then
        read(env_buf(1:env_len), *, iostat=ios) idim
        if ((ios /= 0) .or. (idim < 2)) idim = 10
      end if
    end if

    if (times_in <= 0) then
      call get_environment_variable('TIMES', env_buf, length=env_len, status=env_status)
      if ((env_status == 0) .and. (env_len > 0)) then
        read(env_buf(1:env_len), *, iostat=ios) times
        if ((ios /= 0) .or. (times < 1)) times = 100
      end if
    end if

    call psb_barrier(ctxt)
    if (use_external_matrix) then
      call load_external_matrix(ctxt, matrix_file, matrix_fmt, a, y, x, desc_a, afmt, info)
    else
      call psb_d_gen_pde3d(ctxt,idim,a,y,x,desc_a,afmt,info)
    end if
    if (info /= psb_success_) goto 9999

#ifdef PSB_HAVE_CUDA
    if (use_gpu) then
      select case(psb_toupper(trim(gpu_fmt)))
      case('ELG')
        cuda_sparse_mold => cuda_ell_sparse_mold
      case('CSRG')
        cuda_sparse_mold => cuda_csr_sparse_mold
      case('HDIAG','HDIA')
        cuda_sparse_mold => cuda_hdia_sparse_mold
      case default
        cuda_sparse_mold => cuda_hll_sparse_mold
      end select
      call a%cscnv(info,mold=cuda_sparse_mold)
      if (info /= psb_success_) goto 9999
      call desc_a%cnv(mold=cuda_index_mold)
      if (info /= psb_success_) goto 9999
      call x%cnv(mold=cuda_vector_mold)
      if (info /= psb_success_) goto 9999
      call y%cnv(mold=cuda_vector_mold)
      if (info /= psb_success_) goto 9999
    end if
#endif

    ! Set the comm scheme on the DESCRIPTOR, not the vector. desc_a%comm_type
    ! survives the GPU x%cnv above (which allocates a fresh x%v WITHOUT copying
    ! comm_handle); x then lazy-inits its handle from desc_a%comm_type at the
    ! first halo. Setting psb_comm_set on x%v%comm_handle *before* cnv silently
    ! reverted every backend to baseline on GPU (CPU was unaffected: no cnv).
    call desc_a%set_comm_scheme(comm_type, info)
    if (info /= psb_success_) goto 9999

    ! warm-up
    call psb_spmm(alpha, a, x, beta, y, desc_a, info, doswap=.true.)
    if (info /= psb_success_) goto 9999

    call psb_barrier(ctxt)
    t0 = psb_wtime()
    do i = 1, times
      call psb_spmm(alpha, a, x, beta, y, desc_a, info, doswap=.true.)
      if (info /= psb_success_) exit
    end do
    t1 = psb_wtime()
    if (info /= psb_success_) goto 9999

    dt = t1 - t0
    call psb_amx(ctxt, dt)
    avg_t = dt / real(times, psb_dpk_)

    if (my_rank == psb_root_) then
      write(psb_out_unit,'(/,"SpMV benchmark")')
      write(psb_out_unit,'("  cpu matrix fmt   : ",a)') trim(afmt)
      if (use_gpu) write(psb_out_unit,'("  gpu matrix fmt   : ",a)') trim(psb_toupper(trim(gpu_fmt)))
      if (use_external_matrix) then
        write(psb_out_unit,'("  matrix file     : ",a)') trim(matrix_file)
        write(psb_out_unit,'("  matrix format   : ",a)') trim(matrix_fmt)
      else
        write(psb_out_unit,'("  idim            : ",i0)') idim
      end if
      write(psb_out_unit,'("  global non zeros : ",i0)') a%get_nzeros()
      write(psb_out_unit,'("  global rows     : ",i0)') a%get_nrows()
      write(psb_out_unit,'("  global cols     : ",i0)') a%get_ncols()
      write(psb_out_unit,'("  repetitions     : ",i0)') times
      write(psb_out_unit,'("  comm backend    : ",a)') trim(psb_toupper(trim(comm_mode)))
      write(psb_out_unit,'("  total time [s]  : ",es12.5)') dt
      write(psb_out_unit,'("  avg time   [s]  : ",es12.5)') avg_t
    end if

    call psb_gefree(y, desc_a, info)
    call psb_gefree(x, desc_a, info)
    call psb_spfree(a, desc_a, info)
    call psb_cdfree(desc_a, info)
    call psb_erractionrestore(err_act)
    return

9999 call psb_error(ctxt)
    call psb_error_handler(ctxt, err_act)
  end subroutine run_spmv_kernel

  subroutine load_external_matrix(ctxt, matrix_file, matrix_fmt, a, bv, xv, desc_a, afmt, info)
    type(psb_ctxt_type), intent(in)      :: ctxt
    character(len=*), intent(in)         :: matrix_file
    character(len=*), intent(in)         :: matrix_fmt
    type(psb_dspmat_type), intent(out)   :: a
    type(psb_d_vect_type), intent(out)   :: bv, xv
    type(psb_desc_type), intent(out)     :: desc_a
    character(len=*), intent(in)         :: afmt
    integer(psb_ipk_), intent(out)       :: info

    type(psb_ldspmat_type) :: aux_a
    real(psb_dpk_), allocatable :: rhs_glob(:), x_glob(:)
    integer(psb_lpk_) :: nrows, ncols
    integer(psb_ipk_) :: iam, np

    info = psb_success_
    call psb_info(ctxt, iam, np)

    ! Only the root holds the global matrix; psb_matdist scatters from there.
    ! Reading it on every rank costs one full copy per rank: at 448 ranks over
    ! 4 nodes that is more than 100 GB per node for a 60M nonzero matrix, and
    ! the job is OOM-killed long before it reaches the SpMV.
    if (iam == psb_root_) then
      select case(psb_toupper(trim(matrix_fmt)))
      case('MM')
        call mm_mat_read(aux_a,info,filename=trim(matrix_file))
      case('HB')
        call hb_read(aux_a,info,filename=trim(matrix_file))
      case default
        info = psb_err_internal_error_
      end select
      if (info == psb_success_) then
        nrows = aux_a%get_nrows()
        ncols = aux_a%get_ncols()
        if (nrows /= ncols) then
          write(psb_err_unit,'("Input matrix must be square: ",a)') trim(matrix_file)
          info = psb_err_internal_error_
        end if
      end if
    end if
    call psb_bcast(ctxt, info)
    if (info /= psb_success_) return
    call psb_bcast(ctxt, nrows)
    call psb_bcast(ctxt, ncols)

    call psb_matdist(aux_a, a, ctxt, desc_a, info, fmt=afmt, parts=part_block)
    if (info /= psb_success_) return

    call psb_geall(xv,desc_a,info)
    if (info /= psb_success_) return
    call psb_geall(bv,desc_a,info)
    if (info /= psb_success_) return

    ! Same reasoning: psb_scatter only reads these on the root.
    if (iam == psb_root_) then
      allocate(rhs_glob(nrows), x_glob(ncols), stat=info)
    else
      allocate(rhs_glob(1), x_glob(1), stat=info)
    end if
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      return
    end if
    rhs_glob = done
    x_glob = dzero

    call psb_scatter(rhs_glob,bv,desc_a,info,root=psb_root_)
    if (info /= psb_success_) return
    call psb_scatter(x_glob,xv,desc_a,info,root=psb_root_)
    if (info /= psb_success_) return

    deallocate(rhs_glob, x_glob)
  end subroutine load_external_matrix

end module psb_spmv_overlap_test

program psb_spmv_kernel
  use psb_spmv_overlap_test
  use psb_base_mod
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif
  implicit none

  type(psb_ctxt_type) :: ctxt
  logical :: use_gpu
  integer(psb_ipk_) :: my_rank, np, k
  integer :: ios
  character(len=256) :: arg
  character(len=256) :: matrix_file
  character(len=2)   :: matrix_fmt
  character(len=8)   :: cpu_fmt
  character(len=8)   :: gpu_fmt
  integer(psb_ipk_) :: idim_arg, times_arg
  integer :: kmode
  integer, parameter :: n_comm_modes = 5
  character(len=20), parameter :: comm_modes(n_comm_modes) = [character(len=20) :: &
    & 'P2P', 'NEIGHBOR', 'PNEIGHBOR', 'MPI_GET', 'MPI_PUT']

  idim_arg  = -1
  times_arg = -1

  matrix_file = ''
  matrix_fmt  = 'MM'
  cpu_fmt     = 'CSR'
  gpu_fmt     = 'HLG'

  call psb_init(ctxt)
  call psb_info(ctxt, my_rank, np)

#ifdef PSB_HAVE_CUDA
  use_gpu = .true.
#else
  use_gpu = .false.
#endif

  do k = 1, command_argument_count()
    call get_command_argument(k, arg)
    if (index(psb_toupper(trim(arg)), '--GPU=') == 1) then
      select case (psb_toupper(adjustl(arg(7:len_trim(arg)))))
      case ('TRUE','T','1','YES','Y','ON')
        use_gpu = .true.
      case ('FALSE','F','0','NO','N','OFF')
        use_gpu = .false.
      end select
    else if (index(psb_toupper(trim(arg)), '--MATRIX=') == 1) then
      matrix_file = adjustl(arg(10:len_trim(arg)))
    else if (index(psb_toupper(trim(arg)), '--FMT=') == 1) then
      arg = psb_toupper(adjustl(arg(7:len_trim(arg))))
      if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
    else if (index(psb_toupper(trim(arg)), '--MTX_FMT=') == 1) then
      arg = psb_toupper(adjustl(arg(10:len_trim(arg))))
      if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
    else if (index(psb_toupper(trim(arg)), '--DIM=') == 1) then
      read(arg(7:len_trim(arg)),*,iostat=ios) idim_arg
      if ((ios /= 0) .or. (idim_arg < 2)) idim_arg = -1
    else if (index(psb_toupper(trim(arg)), '--TIMES=') == 1) then
      read(arg(9:len_trim(arg)),*,iostat=ios) times_arg
      if ((ios /= 0) .or. (times_arg < 1)) times_arg = -1
    else if (index(psb_toupper(trim(arg)), '--ITERS=') == 1) then
      read(arg(9:len_trim(arg)),*,iostat=ios) times_arg
      if ((ios /= 0) .or. (times_arg < 1)) times_arg = -1
    else if (index(psb_toupper(trim(arg)), '--CPU_FORMAT=') == 1) then
      cpu_fmt = psb_toupper(adjustl(arg(14:len_trim(arg))))
    else if (index(psb_toupper(trim(arg)), '--CPU_FMT=') == 1) then
      cpu_fmt = psb_toupper(adjustl(arg(11:len_trim(arg))))
    else if (index(psb_toupper(trim(arg)), '--GPU_FORMAT=') == 1) then
      gpu_fmt = psb_toupper(adjustl(arg(14:len_trim(arg))))
    else if (index(psb_toupper(trim(arg)), '--GPU_FMT=') == 1) then
      gpu_fmt = psb_toupper(adjustl(arg(11:len_trim(arg))))
    else if (trim(psb_toupper(arg)) == '--MATRIX') then
      if (k < command_argument_count()) call get_command_argument(k+1,matrix_file)
    else if (trim(psb_toupper(arg)) == '--FMT') then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        arg = psb_toupper(trim(arg))
        if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
      end if
    else if (trim(psb_toupper(arg)) == '--MTX_FMT') then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        arg = psb_toupper(trim(arg))
        if ((trim(arg) == 'MM') .or. (trim(arg) == 'HB')) matrix_fmt = trim(arg)
      end if
    else if (trim(psb_toupper(arg)) == '--DIM') then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        read(arg,*,iostat=ios) idim_arg
        if ((ios /= 0) .or. (idim_arg < 2)) idim_arg = -1
      end if
    else if ((trim(psb_toupper(arg)) == '--TIMES') .or. (trim(psb_toupper(arg)) == '--ITERS')) then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        read(arg,*,iostat=ios) times_arg
        if ((ios /= 0) .or. (times_arg < 1)) times_arg = -1
      end if
    else if ((trim(psb_toupper(arg)) == '--CPU_FORMAT') .or. (trim(psb_toupper(arg)) == '--CPU_FMT')) then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        cpu_fmt = psb_toupper(trim(arg))
      end if
    else if ((trim(psb_toupper(arg)) == '--GPU_FORMAT') .or. (trim(psb_toupper(arg)) == '--GPU_FMT')) then
      if (k < command_argument_count()) then
        call get_command_argument(k+1,arg)
        gpu_fmt = psb_toupper(trim(arg))
      end if
    end if
  end do

#ifdef PSB_HAVE_CUDA
  if (use_gpu) call psb_cuda_init(ctxt)
#else
  use_gpu = .false.
#endif

  if (my_rank == psb_root_) then
    write(psb_out_unit,*) 'Welcome to PSBLAS version: ', psb_version_string_
    write(psb_out_unit,*) 'This is the psb_spmv_kernel sample program'
    write(psb_out_unit,'("GPU enabled          : ",l1)') use_gpu
    write(psb_out_unit,'("Usage: ./psb_spmv_kernel [--gpu=TRUE|FALSE] [--dim=N] [--times=N] ",&
      &"[--cpu_fmt=CSR|COO|CSC|ELL|HLL] [--gpu_fmt=HLL|ELL|CSR|HDIA] [--matrix=<path>] [--fmt=MM|HB] ",&
      &"(runs all comm backends)")')
  end if

  do kmode = 1, n_comm_modes
    if (my_rank == psb_root_) then
      write(psb_out_unit,'(/,"=== Backend sweep: ",a," ===")') trim(comm_modes(kmode))
    end if
    call run_spmv_kernel(ctxt, use_gpu, matrix_file, matrix_fmt, cpu_fmt, gpu_fmt, &
      & idim_arg, times_arg, comm_modes(kmode))
  end do

#ifdef PSB_HAVE_CUDA
  if (use_gpu) call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
end program psb_spmv_kernel
