program psb_comm_cg_test
  use psb_base_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_comm_factory_mod

  implicit none

  type(psb_ctxt_type) :: ctxt
  type(psb_dspmat_type) :: a
  type(psb_desc_type) :: desc_a
  type(psb_d_vect_type) :: b, x
  type(psb_dprec_type)  :: prec

  integer(psb_ipk_) :: info, iam, np
  integer(psb_ipk_) :: idim, itmax, itrace, istop, iter, is
  integer(psb_ipk_) :: iter_arr(3), info_arr(3)
  integer(psb_ipk_) :: scheme_types(3)
  real(psb_dpk_) :: eps, err, t1, t2
  real(psb_dpk_) :: tsolve(3), err_arr(3)
  character(len=25) :: scheme_names(3)
  character(len=5) :: afmt
  character(len=256) :: arg

  info = psb_success_
  afmt = 'CSR'
  idim = 40
  itmax = 500
    itrace = 0
  istop = 2
  eps = 1.d-6
  scheme_types = (/ psb_comm_isend_irecv_, psb_comm_ineighbor_alltoallv_, &
       & psb_comm_persistent_ineighbor_alltoallv_ /)
  scheme_names(1) = 'isend_irecv'
  scheme_names(2) = 'ineighbor_alltoallv'
  scheme_names(3) = 'persistent_ineighbor_a2av'

  call get_command_argument(1,arg)
  if (len_trim(arg) > 0) then
    read(arg,*,iostat=info) idim
    if (info /= 0) then
      idim = 40
      info = psb_success_
    end if
  end if

  call psb_init(ctxt)
  call psb_info(ctxt, iam, np)

  if (iam == psb_root_) then
    write(psb_out_unit,*) 'Welcome to PSBLAS version: ', psb_version_string_
    write(psb_out_unit,*) 'This is the comm/cg test program'
    write(psb_out_unit,'("Grid dimensions      : ",i4," x ",i4," x ",i4)') idim,idim,idim
    write(psb_out_unit,'("Number of processors : ",i0)') np
    write(psb_out_unit,'("Iterative method     : CG")')
    write(psb_out_unit,'("Preconditioner       : NONE")')
    write(psb_out_unit,'(" ")')
  end if

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_d_gen_pde3d(ctxt,idim,a,b,x,desc_a,afmt,info)
  if (info /= psb_success_) goto 9999

  call prec%init(ctxt,'NONE',info)
  if (info /= psb_success_) goto 9999

  call prec%build(a,desc_a,info)
  if (info /= psb_success_) goto 9999

  do is = 1, 3
    call psb_geaxpby(dzero,b,dzero,x,desc_a,info)
    if (info /= psb_success_) goto 9999

    call psb_comm_init(scheme_types(is),x%v%comm_handle,info)
    if (info /= psb_success_) goto 9999

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call psb_krylov('CG',a,prec,b,x,eps,desc_a,info,&
         & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istop)
    t2 = psb_wtime() - t1
    call psb_amx(ctxt,t2)

    tsolve(is) = t2
    iter_arr(is) = iter
    err_arr(is) = err
    info_arr(is) = info

    if (info /= psb_success_) goto 9999
  end do

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("CG solve time by communication scheme")')
    write(psb_out_unit,'("--------------------------------------")')
    do is = 1, 3
      write(psb_out_unit,'(a25,2x,"time=",es12.5,2x,"iter=",i8,2x,"err=",es12.5,2x,"info=",i6)') &
           & trim(scheme_names(is)), tsolve(is), iter_arr(is), err_arr(is), info_arr(is)
    end do
  end if

  call psb_gefree(b,desc_a,info)
  call psb_gefree(x,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call prec%free(info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)
  stop 1

contains

  function b1(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b1

  function b2(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b2

  function b3(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function b3

  function cfun(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
  end function cfun

  function a1(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a1

  function a2(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a2

  function a3(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = done/80
  end function a3

  function gfun(x,y,z) result(val)
    real(psb_dpk_), intent(in) :: x,y,z
    real(psb_dpk_) :: val
    val = dzero
    if (x == done) then
      val = done
    else if (x == dzero) then
      val = exp(y**2-z**2)
    end if
  end function gfun

  subroutine psb_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,afmt,info)
    implicit none
    integer(psb_ipk_), intent(in)     :: idim
    type(psb_dspmat_type), intent(out) :: a
    type(psb_d_vect_type), intent(out) :: xv,bv
    type(psb_desc_type), intent(out)   :: desc_a
    type(psb_ctxt_type), intent(in)    :: ctxt
    integer(psb_ipk_), intent(out)     :: info
    character(len=*), intent(in)       :: afmt

    integer(psb_ipk_), parameter :: nb=20
    real(psb_dpk_) :: zt(nb),x,y,z
    integer(psb_lpk_) :: m,n,glob_row
    integer(psb_ipk_) :: nnz,nlr,i,ii,ib,k
    integer(psb_ipk_) :: ix,iy,iz
    integer(psb_ipk_) :: np, iam, nr, nt
    integer(psb_ipk_) :: icoeff
    integer(psb_lpk_), allocatable :: irow(:),icol(:),myidx(:)
    real(psb_dpk_), allocatable :: val(:)
    real(psb_dpk_) :: deltah, sqdeltah, deltah2
    real(psb_dpk_) :: t0, t1, t2, t3, tasb, talc, ttot, tgen, tcdasb
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name, ch_err, tmpfmt

    info = psb_success_
    name = 'create_matrix'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    deltah   = done/(idim+2)
    sqdeltah = deltah*deltah
    deltah2  = 2.d0*deltah

    m   = idim*idim*idim
    n   = m
    nnz = ((n*9)/(np))
    if(iam == psb_root_) write(psb_out_unit,'("Generating Matrix (size=",i0,")...")')n

    nt = (m+np-1)/np
    nr = max(0,min(nt,m-(iam*nt)))
    nt = nr
    call psb_sum(ctxt,nt)
    if (nt /= m) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,m
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if

    call psb_barrier(ctxt)
    t0 = psb_wtime()
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    if (info == psb_success_) call psb_spall(a,desc_a,info,nnz=nnz)
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

    allocate(val(20*nb),irow(20*nb),icol(20*nb),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif

    myidx = desc_a%get_global_indices()
    nlr = size(myidx)

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    do ii=1, nlr,nb
      ib = min(nb,nlr-ii+1)
      icoeff = 1
      do k=1,ib
        i=ii+k-1
        glob_row=myidx(i)

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

        x = (ix-1)*deltah
        y = (iy-1)*deltah
        z = (iz-1)*deltah
        zt(k) = dzero

        val(icoeff) = -a1(x,y,z)/sqdeltah-b1(x,y,z)/deltah2
        if (ix == 1) then
          zt(k) = gfun(dzero,y,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-2)*idim*idim+(iy-1)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff)  = -a2(x,y,z)/sqdeltah-b2(x,y,z)/deltah2
        if (iy == 1) then
          zt(k) = gfun(x,dzero,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-2)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a3(x,y,z)/sqdeltah-b3(x,y,z)/deltah2
        if (iz == 1) then
          zt(k) = gfun(x,y,dzero)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz-1)
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff)=2.d0*(a1(x,y,z)+a2(x,y,z)+a3(x,y,z))/sqdeltah + cfun(x,y,z)
        icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+iz
        irow(icoeff) = glob_row
        icoeff = icoeff+1

        val(icoeff) = -a3(x,y,z)/sqdeltah+b3(x,y,z)/deltah2
        if (iz == idim) then
          zt(k) = gfun(x,y,done)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+(iy-1)*idim+(iz+1)
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a2(x,y,z)/sqdeltah+b2(x,y,z)/deltah2
        if (iy == idim) then
          zt(k) = gfun(x,done,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = (ix-1)*idim*idim+iy*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
        endif

        val(icoeff) = -a1(x,y,z)/sqdeltah+b1(x,y,z)/deltah2
        if (ix == idim) then
          zt(k) = gfun(done,y,z)*(-val(icoeff)) + zt(k)
        else
          icol(icoeff) = ix*idim*idim+(iy-1)*idim+iz
          irow(icoeff) = glob_row
          icoeff = icoeff+1
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
    call psb_cdasb(desc_a,info)
    tcdasb = psb_wtime()-t1

    call psb_barrier(ctxt)
    t1 = psb_wtime()
    if (info == psb_success_) call psb_spasb(a,desc_a,info,afmt=afmt)
    call psb_barrier(ctxt)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='asb rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (info == psb_success_) call psb_geasb(xv,desc_a,info)
    if (info == psb_success_) call psb_geasb(bv,desc_a,info)
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
      write(psb_out_unit,'("The matrix has been generated and assembled in ",a3," format.")') tmpfmt
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

end program psb_comm_cg_test
