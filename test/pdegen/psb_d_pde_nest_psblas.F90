!
!  Test code for all subroutines in psb_d_nest_psblas_mod.
!  Extends psb_d_pde_nest_first.F90: copies the same 2x2 saddle-point
!  setup (n=10 global DOFs per block), then exercises every operation
!  exported by psb_d_nest_psblas_mod.
!
!  Tested subroutines
!  ------------------
!  T01  psb_d_nest_spmm          y = alpha*A*x + beta*y  (block SpMV)
!  T02  psb_d_nest_geaxpby       y = alpha*x + beta*y
!  T03  psb_d_nest_genrm2        ||x||_2  (function)
!  T04  psb_d_nest_genrm2s       ||x||_2  (subroutine)
!  T05  psb_d_nest_gedot         dot(x, y)
!  T06  psb_d_nest_geamax        ||x||_inf
!  T07  psb_d_nest_geasum        ||x||_1
!  T08  psb_d_nest_gemin         min(x)
!  T09  psb_d_nest_minquotient   min(x/y)
!  T10  psb_d_nest_gemlt         y = y * x  (element-wise)
!  T11  psb_d_nest_gediv         x = x / y  (element-wise; result in x)
!  T12  psb_d_nest_geinv         y = 1/x    (result in y)
!  T13  psb_d_nest_geabs         y = |x|    (result in y)
!  T14  psb_d_nest_geaddconst    z = x + b  (result in z)
!  T15  psb_d_nest_gecmp         z(i)=1 if |x(i)|>=c, 0 otherwise
!  T16  psb_d_nest_mask          mask operation; t=.true. if all satisfied
!  T17  psb_d_nest_upd_xyz       y=alpha*x+beta*y; z=gamma*y_new+delta*z
!
program psb_d_pde_nest_psblas
  use psb_base_mod
  use psb_desc_nest_mod
  use psb_d_nest_mod
  use psb_d_nest_tools_mod, only : psb_geall_nest, psb_geasb_nest, psb_gefree_nest, &
       &                           psb_geins_nest, psb_spall_nest, psb_spins_nest,   &
       &                           psb_spasb_nest
  implicit none

  !------------------------------------------------------------------
  ! Parallel context
  !------------------------------------------------------------------
  type(psb_ctxt_type)  :: ctxt
  integer(psb_ipk_)    :: iam, np, info

  !------------------------------------------------------------------
  ! Problem size
  !------------------------------------------------------------------
  integer(psb_ipk_), parameter :: n = 100
  integer(psb_ipk_)    :: nlr, iloc, i
  integer(psb_lpk_)    :: iglob

  !------------------------------------------------------------------
  ! Per-block descriptors (identical layout as psb_d_pde_nest_first)
  !------------------------------------------------------------------
  type(psb_desc_type)  :: desc1, desc2, desc3, desc4

  !------------------------------------------------------------------
  ! Nested descriptor and nested sparse matrix
  !------------------------------------------------------------------
  type(psb_desc_nest_type)     :: descs
  type(psb_d_nest_sparse_mat)  :: anest

  !------------------------------------------------------------------
  ! Individual sparse matrices (A11 = Laplacian, A12 = I, A21 = I)
  !------------------------------------------------------------------
  type(psb_dspmat_type) :: a11, a12, a21

  !------------------------------------------------------------------
  ! Work nested vectors  (xnest, ynest, znest  reused across tests)
  !------------------------------------------------------------------
  type(psb_d_nest_vect_type) :: xnest, ynest, znest

  !------------------------------------------------------------------
  ! Insert buffers
  !------------------------------------------------------------------
  integer(psb_lpk_) :: grow(1)
  real(psb_dpk_)    :: gval(1)

  !------------------------------------------------------------------
  ! Scalar results and expected values
  !------------------------------------------------------------------
  real(psb_dpk_) :: res, res2, expected
  real(psb_dpk_), parameter :: tol = 1.0e-10_psb_dpk_
  logical :: t_mask

  !------------------------------------------------------------------
  ! Test pass/fail counter
  !------------------------------------------------------------------
  integer(psb_ipk_) :: npass, nfail

  character(len=40) :: name = 'psb_d_pde_nest_psblas'

  !==================================================================
  ! Initialise MPI
  !==================================================================
  call psb_init(ctxt)
  call psb_info(ctxt, iam, np)
  call psb_set_errverbosity(itwo)

  npass = 0
  nfail = 0

  !==================================================================
  ! Build per-block descriptors
  ! Use exact block distribution so total rows = n exactly.
  ! Ceiling division (n+np-1)/np gives nlr*np > n phantom rows
  ! when np does not divide n evenly.
  !==================================================================
  nlr = n / np
  if (iam < mod(n, np)) nlr = nlr + 1_psb_ipk_
  nlr = max(1_psb_ipk_, nlr)
  call psb_cdall(ctxt, desc1, info, nl=nlr)
  call psb_cdall(ctxt, desc2, info, nl=nlr)
  call desc1%clone(desc3, info)

  !==================================================================
  ! A(1,1) = tridiagonal Laplacian
  !==================================================================
  call psb_spall(a11, desc1, info, nnz=3*desc1%get_local_rows())
  do iloc = 1, desc1%get_local_rows()
    call desc1%l2g(iloc, iglob, info)
    grow(1)=iglob; gval(1)=2.0_psb_dpk_
    call psb_spins(1,grow,grow,gval,a11,desc1,info)
    if (iglob>1) then
      grow(1)=iglob; gval(1)=-1.0_psb_dpk_
      call psb_spins(1,grow,[iglob-1_psb_lpk_],gval,a11,desc1,info)
    end if
    if (iglob<n) then
      grow(1)=iglob; gval(1)=-1.0_psb_dpk_
      call psb_spins(1,grow,[iglob+1_psb_lpk_],gval,a11,desc1,info)
    end if
  end do

  !==================================================================
  ! A(1,2) = Identity
  !==================================================================
  call psb_spall(a12, desc2, info, nnz=desc2%get_local_rows())
  do iloc = 1, desc2%get_local_rows()
    call desc2%l2g(iloc, iglob, info)
    grow(1)=iglob; gval(1)=1.0_psb_dpk_
    call psb_spins(1,grow,grow,gval,a12,desc2,info)
  end do

  !==================================================================
  ! A(2,1) = Identity
  !==================================================================
  call psb_spall(a21, desc3, info, nnz=desc3%get_local_rows())
  do iloc = 1, desc3%get_local_rows()
    call desc3%l2g(iloc, iglob, info)
    grow(1)=iglob; gval(1)=1.0_psb_dpk_
    call psb_spins(1,grow,grow,gval,a21,desc3,info)
  end do

  !------------------------------------------------------------------
  ! Finalise descriptors and assemble matrices
  !------------------------------------------------------------------
  call psb_cdasb(desc1, info)
  call psb_cdasb(desc2, info)
  call psb_cdasb(desc3, info)
  call desc2%clone(desc4, info)

  call psb_spasb(a11, desc1, info, dupl=psb_dupl_add_)
  call psb_spasb(a12, desc2, info, dupl=psb_dupl_add_)
  call psb_spasb(a21, desc3, info, dupl=psb_dupl_add_)

  !==================================================================
  ! Assemble nested matrix
  !==================================================================
  anest%nrblocks = 2
  anest%ncblocks = 2
  allocate(anest%mats(2,2), anest%blk_present(2,2), stat=info)
  if (info /= 0) goto 9999
  anest%blk_present = .false.

  call psb_move_alloc(a11, anest%mats(1,1), info)
  call psb_move_alloc(a12, anest%mats(1,2), info)
  call psb_move_alloc(a21, anest%mats(2,1), info)

  anest%blk_present(1,1) = .true.
  anest%blk_present(1,2) = .true.
  anest%blk_present(2,1) = .true.

  !==================================================================
  ! Assemble nested descriptor
  !==================================================================
  descs%nrblocks = 2
  descs%ncblocks = 2
  allocate(descs%descs(2,2), stat=info)
  if (info /= 0) goto 9999

  call desc1%clone(descs%descs(1,1), info)
  call desc2%clone(descs%descs(1,2), info)
  call desc3%clone(descs%descs(2,1), info)
  call desc4%clone(descs%descs(2,2), info)

  !==================================================================
  ! Allocate and assemble work nested vectors
  !==================================================================
  call psb_geall_nest(xnest, descs, info)
  if (info /= psb_success_) goto 9999
  call psb_geasb_nest(xnest, descs, info)
  if (info /= psb_success_) goto 9999

  call psb_geall_nest(ynest, descs, info)
  if (info /= psb_success_) goto 9999
  call psb_geasb_nest(ynest, descs, info)
  if (info /= psb_success_) goto 9999

  call psb_geall_nest(znest, descs, info)
  if (info /= psb_success_) goto 9999
  call psb_geasb_nest(znest, descs, info)
  if (info /= psb_success_) goto 9999

  !==================================================================
  ! T01: psb_d_nest_spmm
  !   x = all 1s;  y = anest * x
  !   Block 1 (Laplacian * ones):  interior rows give 0, boundary rows
  !   give 1 (first) or 1 (last), interior give 0.
  !   Block 2 (A21 * ones) = ones.
  !   Print both result blocks for visual inspection.
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T01: psb_d_nest_spmm  (y = anest * x,  x = all-1s)'

  call set_nest_val(xnest, done)
  call ynest%zero()

  call psb_d_nest_spmm(done, anest, xnest, dzero, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  ! call print_nest_vec(ynest, 'y = anest * [1,1]^T', iam, np, ctxt, descs)

  ! Expected values for y = anest * ones (ng=n ensures exactly n global rows):
  !   Block 1: row 1 and row n give Lap=1, I=1 => y1=2  (2 boundary rows)
  !            rows 2..n-1 give Lap=0, I=1     => y1=1  (n-2 interior rows)
  !   Block 2: I*ones = 1 for all n rows
  !   amax=2, gemin=1, geasum = (n+2) + n = 2n+2
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T01 spmm amax(y)=2', res, 2.0_psb_dpk_, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T01 spmm gemin(y)=1', res, done, tol, npass, nfail, iam)
  res = psb_d_nest_geasum(ynest, descs, info)
  expected = 2.0_psb_dpk_ * real(n, psb_dpk_) + 2.0_psb_dpk_
  call check('T01 spmm geasum(y)=2n+2', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T02: psb_d_nest_geaxpby
  !   x = all 3s,  y = all 2s
  !   y = 2*x + (-1)*y  =>  y = 2*3 - 2 = 4  (all 4s)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T02: psb_d_nest_geaxpby  (y = 2*x - y,  x=3 y=2 => 4)'

  call set_nest_val(xnest, 3.0_psb_dpk_)
  call set_nest_val(ynest, 2.0_psb_dpk_)

  call psb_d_nest_geaxpby(2.0_psb_dpk_, xnest, -done, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 4.0_psb_dpk_
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T02 geaxpby amax(y)=4', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T02 geaxpby amin(y)=4', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T03: psb_d_nest_genrm2
  !   x = all 1s  =>  ||x||_2 = sqrt(2*n) = sqrt(20)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T03: psb_d_nest_genrm2  (x=1 => sqrt(2n))'

  call set_nest_val(xnest, done)

  res      = psb_d_nest_genrm2(xnest, descs, info)
  expected = sqrt(2.0_psb_dpk_ * real(n, psb_dpk_))
  call check('T03 genrm2(ones)', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T04: psb_d_nest_genrm2s  (subroutine form; result must equal T03)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T04: psb_d_nest_genrm2s  (subroutine form of genrm2)'

  call psb_d_nest_genrm2s(res2, xnest, descs, info)
  call check('T04 genrm2s == genrm2', res2, res, tol, npass, nfail, iam)

  !==================================================================
  ! T05: psb_d_nest_gedot
  !   x = all 1s,  y = all 2s  =>  dot = 2 * 2*n = 40
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T05: psb_d_nest_gedot  (x=1 y=2 => 2*2n=40)'

  call set_nest_val(xnest, done)
  call set_nest_val(ynest, 2.0_psb_dpk_)

  res      = psb_d_nest_gedot(xnest, ynest, descs, info)
  expected = 2.0_psb_dpk_ * 2.0_psb_dpk_ * real(n, psb_dpk_)
  call check('T05 gedot', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T06: psb_d_nest_geamax
  !   x = all 5s  =>  ||x||_inf = 5
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T06: psb_d_nest_geamax  (x=5 => 5)'

  call set_nest_val(xnest, 5.0_psb_dpk_)

  res      = psb_d_nest_geamax(xnest, descs, info)
  expected = 5.0_psb_dpk_
  call check('T06 geamax', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T07: psb_d_nest_geasum
  !   x = all 1s  =>  ||x||_1 = 2*n = 20
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T07: psb_d_nest_geasum  (x=1 => 2n=20)'

  call set_nest_val(xnest, done)

  res      = psb_d_nest_geasum(xnest, descs, info)
  expected = 2.0_psb_dpk_ * real(n, psb_dpk_)
  call check('T07 geasum', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T08: psb_d_nest_gemin
  !   x = all 7s  =>  min = 7
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T08: psb_d_nest_gemin  (x=7 => 7)'

  call set_nest_val(xnest, 7.0_psb_dpk_)

  res      = psb_d_nest_gemin(xnest, descs, info)
  expected = 7.0_psb_dpk_
  call check('T08 gemin', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T09: psb_d_nest_minquotient
  !   x = all 3s,  y = all 6s  =>  min(x/y) = 0.5
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T09: psb_d_nest_minquotient  (x=3 y=6 => 0.5)'

  call set_nest_val(xnest, 3.0_psb_dpk_)
  call set_nest_val(ynest, 6.0_psb_dpk_)

  res      = psb_d_nest_minquotient(xnest, ynest, descs, info)
  expected = 0.5_psb_dpk_
  call check('T09 minquotient', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T10: psb_d_nest_gemlt
  !   x = all 2s,  y = all 4s  =>  y = y * x = 8  (result in y)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T10: psb_d_nest_gemlt  (x=2 y=4 => y=8)'

  call set_nest_val(xnest, 2.0_psb_dpk_)
  call set_nest_val(ynest, 4.0_psb_dpk_)

  call psb_d_nest_gemlt(xnest, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 8.0_psb_dpk_
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T10 gemlt amax(y)=8', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T10 gemlt amin(y)=8', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T11: psb_d_nest_gediv
  !   x = all 6s,  y = all 3s  =>  x = x / y = 2  (result in x)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T11: psb_d_nest_gediv  (x=6 y=3 => x=2)'

  call set_nest_val(xnest, 6.0_psb_dpk_)
  call set_nest_val(ynest, 3.0_psb_dpk_)

  call psb_d_nest_gediv(xnest, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 2.0_psb_dpk_
  res = psb_d_nest_geamax(xnest, descs, info)
  call check('T11 gediv amax(x)=2', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(xnest, descs, info)
  call check('T11 gediv amin(x)=2', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T12: psb_d_nest_geinv
  !   x = all 4s  =>  y = 1/x = 0.25  (result in y)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T12: psb_d_nest_geinv  (x=4 => y=0.25)'

  call set_nest_val(xnest, 4.0_psb_dpk_)

  call psb_d_nest_geinv(xnest, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 0.25_psb_dpk_
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T12 geinv amax(y)=0.25', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T12 geinv amin(y)=0.25', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T13: psb_d_nest_geabs
  !   x = all -3s  =>  y = |x| = 3  (result in y)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T13: psb_d_nest_geabs  (x=-3 => y=3)'

  call set_nest_val(xnest, -3.0_psb_dpk_)

  call psb_d_nest_geabs(xnest, ynest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 3.0_psb_dpk_
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T13 geabs amax(y)=3', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T13 geabs amin(y)=3', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T14: psb_d_nest_geaddconst
  !   x = all 2s,  b = 7.0  =>  z = x + 7 = 9  (result in z)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T14: psb_d_nest_geaddconst  (x=2 b=7 => z=9)'

  call set_nest_val(xnest, 2.0_psb_dpk_)

  call psb_d_nest_geaddconst(xnest, 7.0_psb_dpk_, znest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 9.0_psb_dpk_
  res = psb_d_nest_geamax(znest, descs, info)
  call check('T14 geaddconst amax(z)=9', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(znest, descs, info)
  call check('T14 geaddconst amin(z)=9', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T15a: psb_d_nest_gecmp — entries satisfy threshold
  !   x = all 3s,  c = 2.0  =>  z(i)=1  (since |3| >= 2)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T15a: psb_d_nest_gecmp  (x=3 c=2 => z=1)'

  call set_nest_val(xnest, 3.0_psb_dpk_)

  call psb_d_nest_gecmp(xnest, 2.0_psb_dpk_, znest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = done
  res = psb_d_nest_geamax(znest, descs, info)
  call check('T15a gecmp amax(z)=1', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(znest, descs, info)
  call check('T15a gecmp amin(z)=1', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T15b: psb_d_nest_gecmp — entries do not satisfy threshold
  !   x = all 1s,  c = 2.0  =>  z(i)=0  (since |1| < 2)
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T15b: psb_d_nest_gecmp  (x=1 c=2 => z=0)'

  call set_nest_val(xnest, done)

  call psb_d_nest_gecmp(xnest, 2.0_psb_dpk_, znest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = dzero
  res = psb_d_nest_geamax(znest, descs, info)
  call check('T15b gecmp amax(z)=0', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! T16: psb_d_nest_mask
  !   Semantics: mask(c, x, m, t)
  !     c = values to test (first arg)
  !     x = constraint-type indicators (second arg):
  !           2  => satisfied if c(i) > 0
  !           1  => satisfied if c(i) >= 0
  !          -1  => satisfied if c(i) <= 0
  !          -2  => satisfied if c(i) < 0
  !     m = output mask (0=satisfied, 1=violated)
  !     t = .true. iff all entries satisfied
  !
  !   Case: c = all +3 (positive), x = all 2 (check > 0)  => m=0 t=T
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') 'T16: psb_d_nest_mask  (c=3 x=2 => m=0 t=.true.)'

  call set_nest_val(xnest, 3.0_psb_dpk_)   ! values to test
  call set_nest_val(ynest, 2.0_psb_dpk_)   ! constraint indicators (type 2: check > 0)

  call psb_d_nest_mask(xnest, ynest, znest, t_mask, descs, info)
  if (info /= psb_success_) goto 9999

  if (iam == 0) then
    if (t_mask) then
      write(*,'(a)') '  T16 mask: t=.true.  PASS  (all constraints satisfied)'
      npass = npass + 1
    else
      write(*,'(a)') '  T16 mask: t=.false. FAIL  (expected .true.)'
      nfail = nfail + 1
    end if
  end if

  !------------------------------------------------------------------
  ! T16b: c = all -3 (negative), x = all 2 (check > 0)  => m=1 t=F
  !------------------------------------------------------------------
  if (iam == 0) write(*,'(a)') 'T16b: psb_d_nest_mask  (c=-3 x=2 => m=1 t=.false.)'

  call set_nest_val(xnest, -3.0_psb_dpk_)  ! values (negative)
  call set_nest_val(ynest,  2.0_psb_dpk_)  ! indicators (type 2: check > 0)

  call psb_d_nest_mask(xnest, ynest, znest, t_mask, descs, info)
  if (info /= psb_success_) goto 9999

  if (iam == 0) then
    if (.not. t_mask) then
      write(*,'(a)') '  T16b mask: t=.false. PASS  (all constraints violated)'
      npass = npass + 1
    else
      write(*,'(a)') '  T16b mask: t=.true.  FAIL  (expected .false.)'
      nfail = nfail + 1
    end if
  end if

  !==================================================================
  ! T17: psb_d_nest_upd_xyz
  !   Computes:  y_new = alpha*x + beta*y
  !              z_new = gamma*y_new + delta*z
  !
  !   x=1, y=2, z=3, alpha=2, beta=3, gamma=4, delta=5
  !   => y_new = 2*1 + 3*2 = 8
  !   => z_new = 4*8 + 5*3 = 47
  !==================================================================
  if (iam == 0) write(*,'(/,a)') repeat('=',60)
  if (iam == 0) write(*,'(a)') &
       'T17: psb_d_nest_upd_xyz  (x=1 y=2 z=3 a=2 b=3 g=4 d=5 => y=8 z=47)'

  call set_nest_val(xnest, done)
  call set_nest_val(ynest, 2.0_psb_dpk_)
  call set_nest_val(znest, 3.0_psb_dpk_)

  call psb_d_nest_upd_xyz(2.0_psb_dpk_, 3.0_psb_dpk_, &
       &                  4.0_psb_dpk_, 5.0_psb_dpk_, &
       &                  xnest, ynest, znest, descs, info)
  if (info /= psb_success_) goto 9999

  expected = 8.0_psb_dpk_
  res = psb_d_nest_geamax(ynest, descs, info)
  call check('T17 upd_xyz amax(y)=8', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(ynest, descs, info)
  call check('T17 upd_xyz amin(y)=8', res, expected, tol, npass, nfail, iam)

  expected = 47.0_psb_dpk_
  res = psb_d_nest_geamax(znest, descs, info)
  call check('T17 upd_xyz amax(z)=47', res, expected, tol, npass, nfail, iam)
  res = psb_d_nest_gemin(znest, descs, info)
  call check('T17 upd_xyz amin(z)=47', res, expected, tol, npass, nfail, iam)

  !==================================================================
  ! Summary
  !==================================================================
  if (iam == 0) then
    write(*,'(/,a)') repeat('=',60)
    write(*,'(a,i0,a,i0,a)') &
         ' RESULTS: ', npass, ' passed,  ', nfail, ' failed'
    write(*,'(a)') repeat('=',60)
  end if

  !==================================================================
  ! Clean up
  !==================================================================
  call psb_gefree_nest(xnest, descs, info)
  call psb_gefree_nest(ynest, descs, info)
  call psb_gefree_nest(znest, descs, info)

  call psb_cdfree(desc1, info)
  call psb_cdfree(desc2, info)
  call psb_cdfree(desc3, info)
  call psb_cdfree(desc4, info)

  call psb_exit(ctxt)
  stop

9999 continue
  write(psb_err_unit,*) trim(name), ': error info=', info, ' rank=', iam
  call psb_error(ctxt)
  call psb_exit(ctxt)
  stop

contains

  !------------------------------------------------------------------
  ! Set every local entry of every block to val
  !------------------------------------------------------------------
  subroutine set_nest_val(v, val)
    use psb_base_mod
    type(psb_d_nest_vect_type), intent(inout) :: v
    real(psb_dpk_),             intent(in)    :: val
    integer(psb_ipk_) :: k, linfo
    linfo = 0
    do k = 1, v%nblocks
      call v%vects(k)%set(val, linfo)
    end do
  end subroutine set_nest_val

  !------------------------------------------------------------------
  ! Scalar pass/fail check with tolerance
  !------------------------------------------------------------------
  subroutine check(label, got, expected, tol, np_, nf_, myrank)
    use psb_base_mod
    character(len=*),  intent(in)    :: label
    real(psb_dpk_),    intent(in)    :: got, expected, tol
    integer(psb_ipk_), intent(inout) :: np_, nf_
    integer(psb_ipk_), intent(in)    :: myrank

    if (myrank /= 0) return
    if (abs(got - expected) <= tol * max(done, abs(expected))) then
      write(*,'(2x,a,a,f16.10,a,f16.10)') &
           'PASS  ', trim(label)//'  got=', got, '  exp=', expected
      np_ = np_ + 1
    else
      write(*,'(2x,a,a,f16.10,a,f16.10)') &
           'FAIL  ', trim(label)//'  got=', got, '  exp=', expected
      nf_ = nf_ + 1
    end if
  end subroutine check

  !------------------------------------------------------------------
  ! Print every block of a nested vector (one rank at a time).
  ! Each process flushes stdout before the barrier so that buffered
  ! output does not bleed into the next process's print window.
  !------------------------------------------------------------------
  subroutine print_nest_vec(v, label, myrank, nprocs, myctxt, ds)
    use psb_base_mod
    use iso_fortran_env, only: output_unit
    type(psb_d_nest_vect_type), intent(inout) :: v
    character(len=*),           intent(in)    :: label
    integer(psb_ipk_),          intent(in)    :: myrank, nprocs
    type(psb_ctxt_type),        intent(in)    :: myctxt
    type(psb_desc_nest_type),   intent(in)    :: ds

    integer(psb_ipk_) :: blk, ip, k, nr, linfo
    real(psb_dpk_), allocatable :: vals(:)

    do blk = 1, v%nblocks
      nr = ds%descs(blk,blk)%get_local_rows()
      do ip = 0, nprocs-1
        call psb_barrier(myctxt)
        if (myrank == ip) then
          write(*,'(a,a,a,i0,a)') '  [', trim(label), '] block ', blk, ':'
          linfo = 0
          allocate(vals(nr), stat=linfo)
          if (linfo == 0) vals = v%vects(blk)%get_vect()
          do k = 1, nr
            write(*,'(4x,i4,f14.6)') k, vals(k)
          end do
          deallocate(vals)
          flush(output_unit)
        end if
      end do
      call psb_barrier(myctxt)
    end do
  end subroutine print_nest_vec

end program psb_d_pde_nest_psblas
