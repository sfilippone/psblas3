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
!    Moved here from MLD2P4, original copyright below.
!
!
!
!                             MLD2P4  version 2.2
!    MultiLevel Domain Decomposition Parallel Preconditioners Package
!               based on PSBLAS (Parallel Sparse BLAS version 3.5)
!
!    (C) Copyright 2008-2018
!
!        Salvatore Filippone
!        Pasqua D'Ambra
!        Daniela di Serafino
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the MLD2P4 group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
! File: psb_diluk_fact.f90
!
! Subroutine: psb_diluk_fact
! Version:    real
! Contains:   psb_diluk_factint, iluk_copyin, iluk_fact, iluk_copyout.
!
!  This routine computes either the ILU(k) or the MILU(k) factorization of the
!  diagonal blocks of a distributed matrix. These factorizations are used to
!  build the 'base preconditioner' (block-Jacobi preconditioner/solver,
!  Additive Schwarz preconditioner) corresponding to a certain level of a
!  multilevel preconditioner.
!
!  Details on the above factorizations can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix is stored into a and blck, as specified in
!  the description of the arguments below. The storage format for both the L and
!  U factors is CSR. The diagonal of the U factor is stored separately (actually,
!  the inverse of the diagonal entries is stored; this is then managed in the solve
!  stage associated to the ILU(k)/MILU(k) factorization).
!
!
! Arguments:
!    fill_in -  integer, input.
!               The fill-in level k in ILU(k)/MILU(k).
!    ialg    -  integer, input.
!               The type of incomplete factorization to be performed.
!               The ILU(k) factorization is computed if ialg = 1 (= psb_ilu_n_);
!               the MILU(k) one if ialg = 2 (= psb_milu_n_); other values are
!               not allowed.
!    a       -  type(psb_dspmat_type), input.
!               The sparse matrix structure containing the local matrix.
!               Note that if the 'base' Additive Schwarz preconditioner
!               has overlap greater than 0 and the matrix has not been reordered
!               (see psb_fact_bld), then a contains only the 'original' local part
!               of the distributed matrix, i.e. the rows of the matrix held
!               by the calling process according to the initial data distribution.
!    l       -  type(psb_dspmat_type), input/output.
!               The L factor in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    u       -  type(psb_dspmat_type), input/output.
!               The U factor (except its diagonal) in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    d       -  real(psb_dpk_), dimension(:), input/output.
!               The inverse of the diagonal entries of the U factor in the incomplete
!               factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    info    -  integer, output.
!               Error code.
!    blck    -  type(psb_dspmat_type), input, optional, target.
!               The sparse matrix structure containing the remote rows of the
!               distributed matrix, that have been retrieved by psb_as_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 or the matrix has been reordered
!               (see psb_fact_bld), then blck does not contain any row.
!
subroutine psb_diluk_fact(fill_in,ialg,a,l,u,d,info,blck,shft)

  use psb_base_mod
  use psb_d_ilu_fact_mod, psb_protect_name => psb_diluk_fact

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)       :: fill_in, ialg
  integer(psb_ipk_), intent(out)      :: info
  type(psb_dspmat_type),intent(in)    :: a
  type(psb_dspmat_type),intent(inout) :: l,u
  type(psb_dspmat_type),intent(in), optional, target :: blck
  real(psb_dpk_), intent(inout)    ::  d(:)
  real(psb_dpk_), intent(in), optional :: shft
  !     Local Variables
  integer(psb_ipk_)   :: l1, l2, m, err_act

  real(psb_dpk_) :: shft_
  type(psb_dspmat_type), pointer  :: blck_
  type(psb_d_csr_sparse_mat)      :: ll, uu
  character(len=20)   :: name, ch_err

  name='psb_diluk_fact'
  info = psb_success_
  call psb_erractionsave(err_act)

  !
  ! Point to / allocate memory for the incomplete factorization
  !
  if (present(blck)) then
    blck_ => blck
  else
    allocate(blck_,stat=info)
    if (info == psb_success_) call blck_%allocate(izero,izero,info,ione,type='CSR') 
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  endif
  if (present(shft)) then
    shft_ = shft
  else
    shft_ = dzero
  end if

  m = a%get_nrows() + blck_%get_nrows()
  if ((m /= l%get_nrows()).or.(m /= u%get_nrows()).or.&
       & (m > size(d))    ) then
    write(0,*) 'Wrong allocation status for L,D,U? ',&
         & l%get_nrows(),size(d),u%get_nrows()
    info = -1
    return
  end if

  call l%mv_to(ll)
  call u%mv_to(uu)

  !
  ! Compute the ILU(k) or the MILU(k) factorization, depending on ialg
  !
  call psb_diluk_factint(fill_in,ialg,a,blck_,&
       & d,ll%val,ll%ja,ll%irp,uu%val,uu%ja,uu%irp,l1,l2,info,shft_)
  if (info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_diluk_factint'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  !
  ! Store information on the L and U sparse matrices
  !
  call l%mv_from(ll)
  call l%set_triangle()
  call l%set_unit()
  call l%set_lower()
  call u%mv_from(uu)
  call u%set_triangle()
  call u%set_unit()
  call u%set_upper()

  !
  ! Nullify pointer / deallocate memory
  !
  if (present(blck)) then
    blck_ => null()
  else
    call blck_%free()
    deallocate(blck_,stat=info)
    if(info.ne.0) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  endif


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  !
  ! Subroutine: psb_diluk_factint
  ! Version:    real
  ! Note: internal subroutine of psb_diluk_fact
  !
  !  This routine computes either the ILU(k) or the MILU(k) factorization of the
  !  diagonal blocks of a distributed matrix. These factorizations are used to build
  !  the 'base preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a certain level of a multilevel preconditioner.
  !
  !  The local matrix is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(k)/MILU(k) factorization).
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k)/MILU(k).
  !    ialg    -  integer, input.
  !               The type of incomplete factorization to be performed.
  !               The MILU(k) factorization is computed if ialg = 2 (= psb_milu_n_);
  !               the ILU(k) factorization otherwise.
  !    m       -  integer, output.
  !               The total number of rows of the local matrix to be factorized,
  !               i.e. ma+mb.
  !    a       -  type(psb_dspmat_type), input.
  !               The sparse matrix structure containing the local matrix.
  !               Note that, if the 'base' Additive Schwarz preconditioner
  !               has overlap greater than 0 and the matrix has not been reordered
  !               (see psb_fact_bld), then a contains only the 'original' local part
  !               of the distributed matrix, i.e. the rows of the matrix held
  !               by the calling process according to the initial data distribution.
  !    b       -  type(psb_dspmat_type), input.
  !               The sparse matrix structure containing the remote rows of the
  !               distributed matrix, that have been retrieved by psb_as_bld
  !               to build an Additive Schwarz base preconditioner with overlap
  !               greater than 0. If the overlap is 0 or the matrix has been reordered
  !               (see psb_fact_bld), then b does not contain   any row.
  !    d       -  real(psb_dpk_), dimension(:), output.
  !               The inverse of the diagonal entries of the U factor in the incomplete
  !               factorization.
  !    laspk   -  real(psb_dpk_), dimension(:), input/output.
  !               The L factor in the incomplete factorization.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               according to the CSR storage format.
  !    lia2    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor in laspk, according to the CSR storage format.
  !    uval   -   real(psb_dpk_), dimension(:), input/output.
  !               The U factor in the incomplete factorization.
  !               The entries of U are stored according to the CSR format.
  !    uja    -   integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor,
  !               according to the CSR storage format.
  !    uirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor in uval, according to the CSR storage format.
  !    l1      -  integer, output
  !               The number of nonzero entries in laspk.
  !    l2      -  integer, output
  !               The number of nonzero entries in uval.
  !    info    -  integer, output.
  !               Error code.
  !
  subroutine psb_diluk_factint(fill_in,ialg,a,b,&
       & d,lval,lja,lirp,uval,uja,uirp,l1,l2,info,shft)

    use psb_base_mod

    implicit none

  ! Arguments
    integer(psb_ipk_), intent(in)                 :: fill_in, ialg
    type(psb_dspmat_type),intent(in)              :: a,b
    integer(psb_ipk_),intent(inout)               :: l1,l2,info
    integer(psb_ipk_), allocatable, intent(inout) :: lja(:),lirp(:),uja(:),uirp(:)
    real(psb_dpk_), allocatable, intent(inout) :: lval(:),uval(:)
    real(psb_dpk_), intent(inout)              :: d(:)
    real(psb_dpk_), intent(in)       :: shft

  ! Local variables
    integer(psb_ipk_) :: ma,mb,i, ktrw,err_act,nidx, m
    integer(psb_ipk_), allocatable   :: uplevs(:), rowlevs(:),idxs(:)
    real(psb_dpk_), allocatable   :: row(:)
    type(psb_i_heap)             :: heap
    type(psb_d_coo_sparse_mat) :: trw
    character(len=20), parameter :: name='psb_diluk_factint'
    character(len=20)            :: ch_err

    info=psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if


    select case(ialg)
    case(psb_ilu_n_,psb_milu_n_)
      ! Ok
    case default
      info=psb_err_input_asize_invalid_i_
      call psb_errpush(info,name,&
           & i_err=(/itwo,ialg,izero,izero,izero/))
      goto 9999
    end select
    if (fill_in < 0) then
      info=psb_err_input_asize_invalid_i_
      call psb_errpush(info,name, &
           & i_err=(/ione,fill_in,izero,izero,izero/))
      goto 9999
    end if

    ma = a%get_nrows()
    mb = b%get_nrows()
    m  = ma+mb

    !
    ! Allocate a temporary buffer for the iluk_copyin function
    !

    call trw%allocate(izero,izero,ione)
    if (info == psb_success_) call psb_ensure_size(m+1,lirp,info)
    if (info == psb_success_) call psb_ensure_size(m+1,uirp,info)

    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_sp_all')
      goto 9999
    end if

    l1=0
    l2=0
    lirp(1) = 1
    uirp(1) = 1

    !
    ! Allocate memory to hold the entries of a row and the corresponding
    ! fill levels
    !
    allocate(uplevs(size(uval)),rowlevs(m),row(m),stat=info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

    uplevs(:)  = m+1
    row(:)     = dzero
    rowlevs(:) = -(m+1)

    !
    ! Cycle over the matrix rows
    !
    do i = 1, m

      !
      ! At each iteration of the loop we keep in a heap the column indices
      ! affected by the factorization. The heap is initialized and filled
      ! in the iluk_copyin routine, and updated during the elimination, in
      ! the iluk_fact routine. The heap is ideal because at each step we need
      ! the lowest index, but we also need to insert new items, and the heap
      ! allows to do both in log time.
      !
      d(i) = dzero
      if (i<=ma) then
        !
        ! Copy into trw the i-th local row of the matrix, stored in a
        !
        call iluk_copyin(i,ma,a,ione,m,row,rowlevs,heap,ktrw,trw,info,shft)
      else
        !
        ! Copy into trw the i-th local row of the matrix, stored in b
        ! (as (i-ma)-th row)
        !
        call iluk_copyin(i-ma,mb,b,ione,m,row,rowlevs,heap,ktrw,trw,info,shft)
      endif

      ! Do an elimination step on the current row. It turns out we only
      ! need to keep track of fill levels for the upper triangle, hence we
      ! do not have a lowlevs variable.
      !
      if (info == psb_success_) call iluk_fact(fill_in,i,row,rowlevs,heap,&
           & d,uja,uirp,uval,uplevs,nidx,idxs,info)
      !
      ! Copy the row into lval/d(i)/uval
      !
      if (info == psb_success_) call iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
           & l1,l2,lja,lirp,lval,d,uja,uirp,uval,uplevs,info)
      if (info /= psb_success_) then
        info=psb_err_internal_error_
        call psb_errpush(info,name,a_err='Copy/factor loop')
        goto 9999
      end if
    end do

    !
    ! And we're done, so deallocate the memory
    !
    deallocate(uplevs,rowlevs,row,stat=info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Deallocate')
      goto 9999
    end if
    if (info == psb_success_) call trw%free()
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_diluk_factint

  !
  ! Subroutine: iluk_copyin
  ! Version:    real
  ! Note: internal subroutine of psb_diluk_fact
  !
  !  This routine copies a row of a sparse matrix A, stored in the sparse matrix
  !  structure a, into the array row and stores into a heap the column indices of
  !  the nonzero entries of the copied row. The output array row is such that it
  !  contains a full row of A, i.e. it contains also the zero entries of the row.
  !  This is useful for the elimination step performed by iluk_fact after the call
  !  to iluk_copyin (see psb_iluk_factint).
  !  The routine also sets to zero the entries of the array rowlevs corresponding
  !  to the nonzero entries of the copied row (see the description of the arguments
  !  below).
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied, row by row, into the array row, through successive calls to
  !  ilu_copyin.
  !
  !  This routine is used by psb_diluk_factint in the computation of the
  !  ILU(k)/MILU(k) factorization of a local sparse matrix.
  !
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from the
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(psb_dspmat_type), input.
  !               The sparse matrix structure containing the row to be copied.
  !    jmin    -  integer, input.
  !               The minimum valid column index.
  !    jmax    -  integer, input.
  !               The maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    row     -  real(psb_dpk_), dimension(:), input/output.
  !               In input it is the null vector (see psb_iluk_factint and
  !               iluk_copyout). In output it contains the row extracted
  !               from the matrix A. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = -(m+1) for k=1,...,m. In output
  !               rowlevs(k) = 0 for 1 <= k <= jmax and A(i,k) /= 0, for
  !               future use in iluk_fact.
  !    heap    -  type(psb_i_heap), input/output.
  !               The heap containing the column indices of the nonzero
  !               entries in the array row.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by psb_init_heap inside this
  !               routine.
  !    ktrw    -  integer, input/output.
  !               The index identifying the last entry taken from the
  !               staging buffer trw. See below.
  !    trw     -  type(psb_dspmat_type), input/output.
  !               A staging buffer. If the matrix A is not in CSR format, we use
  !               the psb_sp_getblk routine and store its output in trw; when we
  !               need to call psb_sp_getblk we do it for a block of rows, and then
  !               we consume them from trw in successive calls to this routine,
  !               until we empty the buffer. Thus we will make a call to psb_sp_getblk
  !               every nrb calls to copyin. If A is in CSR format it is unused.
  !
  subroutine iluk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw,info,shft)

    use psb_base_mod

    implicit none

  ! Arguments
    type(psb_dspmat_type), intent(in)         :: a
    type(psb_d_coo_sparse_mat), intent(inout) :: trw
    integer(psb_ipk_), intent(in)           :: i,m,jmin,jmax
    integer(psb_ipk_), intent(inout)        :: ktrw,info
    integer(psb_ipk_), intent(inout)        :: rowlevs(:)
    real(psb_dpk_), intent(inout)          :: row(:)
    type(psb_i_heap), intent(inout)         :: heap
    real(psb_dpk_), intent(in)       :: shft


  ! Local variables
    integer(psb_ipk_)             :: k,j,irb,err_act,nz
    integer(psb_ipk_), parameter  :: nrb=40
    character(len=20), parameter  :: name='iluk_copyin'
    character(len=20)             :: ch_err

    info=psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if
    call heap%init(info)

    select type (aa=> a%a)
    type is (psb_d_csr_sparse_mat)
      !
      ! Take a fast shortcut if the matrix is stored in CSR format
      !

      do j = aa%irp(i), aa%irp(i+1) - 1
        k          = aa%ja(j)
        if ((jmin<=k).and.(k<=jmax)) then
          row(k)     = aa%val(j)
          if (k==i) row(k) = row(k) + shft
          rowlevs(k) = 0
          call heap%insert(k,info)
        end if
      end do

    class default

      !
      ! Otherwise use psb_sp_getblk, slower but able (in principle) of
      ! handling any format. In this case, a block of rows is extracted
      ! instead of a single row, for performance reasons, and these
      ! rows are copied one by one into the array row, through successive
      ! calls to iluk_copyin.
      !

      if ((mod(i,nrb) == 1).or.(nrb == 1)) then
        irb = min(m-i+1,nrb)
        call aa%csget(i,i+irb-1,trw,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          ch_err='psb_sp_getblk'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        ktrw=1
      end if
      nz = trw%get_nzeros()
      do
        if (ktrw > nz) exit
        if (trw%ia(ktrw) > i) exit
        k          = trw%ja(ktrw)
        if ((jmin<=k).and.(k<=jmax)) then
          row(k)     = trw%val(ktrw)
          if (k==i) row(k) = row(k) + shft
          rowlevs(k) = 0
          call heap%insert(k,info)
        end if
        ktrw       = ktrw + 1
      enddo
    end select
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine iluk_copyin

  !
  ! Subroutine: iluk_fact
  ! Version:    real
  ! Note: internal subroutine of psb_diluk_fact
  !
  !  This routine does an elimination step of the ILU(k) factorization on a
  !  single matrix row (see the calling routine psb_iluk_factint).
  !
  !  This step is also the base for a MILU(k) elimination step on the row (see
  !  iluk_copyout). This routine is used by psb_diluk_factint in the computation
  !  of the ILU(k)/MILU(k) factorization of a local sparse matrix.
  !
  !  NOTE: it turns out we only need to keep track of the fill levels for
  !  the upper triangle.
  !
  !
  ! Arguments
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k).
  !    i       -  integer, input.
  !               The local index of the row to which the factorization is
  !               applied.
  !    row     -  real(psb_dpk_), dimension(:), input/output.
  !               In input it contains the row to which the elimination step
  !               has to be applied. In output it contains the row after the
  !               elimination step. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = 0 if the k-th entry of the row is
  !               nonzero, and rowlevs(k) = -(m+1) otherwise. In output
  !               rowlevs(k) contains the fill kevel of the k-th entry of
  !               the row after the current elimination step; rowlevs(k) = -(m+1)
  !               means that the k-th row entry is zero throughout the elimination
  !               step.
  !    heap    -  type(psb_i_heap), input/output.
  !               The heap containing the column indices of the nonzero entries
  !               in the processed row. In input it contains the indices concerning
  !               the row before the elimination step, while in output it contains
  !               the indices concerning the transformed row.
  !    d       -  real(psb_dpk_), input.
  !               The inverse of the diagonal entries of the part of the U factor
  !               above the current row (see iluk_copyout).
  !    uja    -  integer, dimension(:), input.
  !               The column indices of the nonzero entries of the part of the U
  !               factor above the current row, stored in uval row by row (see
  !               iluk_copyout, called by psb_diluk_factint), according to the CSR
  !               storage format.
  !    uirp    -  integer, dimension(:), input.
  !               The indices identifying the first nonzero entry of each row of
  !               the U factor above the current row, stored in uval row by row
  !               (see iluk_copyout, called by psb_diluk_factint), according to
  !               the CSR storage format.
  !    uval   -  real(psb_dpk_), dimension(:), input.
  !               The entries of the U factor above the current row (except the
  !               diagonal ones), stored according to the CSR format.
  !    uplevs  -  integer, dimension(:), input.
  !               The fill levels of the nonzero entries in the part of the
  !               U factor above the current row.
  !    nidx    -  integer, output.
  !               The number of entries of the array row that have been
  !               examined during the elimination step. This will be used
  !               by the routine iluk_copyout.
  !    idxs    -  integer, dimension(:), allocatable, input/output.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step.This will be used by
  !               by the routine iluk_copyout.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, done by this routine.
  !
  subroutine iluk_fact(fill_in,i,row,rowlevs,heap,d,&
       & uja,uirp,uval,uplevs,nidx,idxs,info)

    use psb_base_mod

    implicit none

  ! Arguments
    type(psb_i_heap), intent(inout)               :: heap
    integer(psb_ipk_), intent(in)                 :: i, fill_in
    integer(psb_ipk_), intent(inout)              :: nidx,info
    integer(psb_ipk_), intent(inout)              :: rowlevs(:)
    integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
    integer(psb_ipk_), intent(inout)              :: uja(:),uirp(:),uplevs(:)
    real(psb_dpk_), intent(inout)              :: row(:), uval(:),d(:)

    ! Local variables
    integer(psb_ipk_)   :: k,j,lrwk,jj,lastk, iret
    real(psb_dpk_)   :: rwk

    info = psb_success_
    if (.not.allocated(idxs)) then
      allocate(idxs(200),stat=info)
      if (info /= psb_success_) return
    endif
    nidx = 0
    lastk = -1

    !
    ! Do while there are indices to be processed
    !
    do
      ! Beware: (iret < 0) means that the heap is empty, not an error.
      call heap%get_first(k,iret)
      if (iret < 0) return

      !
      ! Just in case an index has been put on the heap more than once.
      !
      if (k == lastk) cycle

      lastk = k
      nidx = nidx + 1
      if (nidx>size(idxs)) then
        call psb_realloc(nidx+psb_heap_resize,idxs,info)
        if (info /= psb_success_) return
      end if
      idxs(nidx) = k

      if ((row(k) /= dzero).and.(rowlevs(k) <= fill_in).and.(k<i)) then
        !
        ! Note: since U is scaled while copying it out (see iluk_copyout),
        ! we can use rwk in the update below
        !
        rwk    = row(k)
        row(k) = row(k) * d(k)    ! d(k) == 1/a(k,k)
        lrwk   = rowlevs(k)

        do jj=uirp(k),uirp(k+1)-1
          j = uja(jj)
          if (j<=k) then
            info = -i
            return
          endif
          !
          ! Insert the index into the heap for further processing.
          ! The fill levels are initialized to a negative value. If we find
          ! one, it means that it is an as yet untouched index, so we need
          ! to insert it; otherwise it is already on the heap, there is no
          ! need to insert it more than once.
          !
          if (rowlevs(j)<0) then
            call heap%insert(j,info)
            if (info /= psb_success_) return
            rowlevs(j) = abs(rowlevs(j))
          end if
          !
          ! Update row(j) and the corresponding fill level
          !
          row(j)     = row(j) - rwk * uval(jj)
          rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
        end do

      end if
    end do

  end subroutine iluk_fact

  !
  ! Subroutine: iluk_copyout
  ! Version:    real
  ! Note: internal subroutine of psb_diluk_fact
  !
  !  This routine copies a matrix row, computed by iluk_fact by applying an
  !  elimination step of the ILU(k) factorization, into the arrays lval, uval,
  !  d, corresponding to the L factor, the U factor and the diagonal of U,
  !  respectively.
  !
  !  Note that
  !  - the part of the row stored into uval is scaled by the corresponding diagonal
  !    entry, according to the LDU form of the incomplete factorization;
  !  - the inverse of the diagonal entries of U is actually stored into d; this is
  !    then managed in the solve stage associated to the ILU(k)/MILU(k) factorization;
  !  - if the MILU(k) factorization has been required (ialg == psb_milu_n_), the
  !    row entries discarded because their fill levels are too high are added to
  !    the diagonal entry of the row;
  !  - the row entries are stored in lval and uval according to the CSR format;
  !  - the arrays row and rowlevs are re-initialized for future use in psb_iluk_fact
  !    (see also iluk_copyin and iluk_fact).
  !
  !  This routine is used by psb_diluk_factint in the computation of the
  !  ILU(k)/MILU(k) factorization of a local sparse matrix.
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in level k in ILU(k)/MILU(k).
  !    ialg    -  integer, input.
  !               The type of incomplete factorization considered. The MILU(k)
  !               factorization is computed if ialg = 2 (= psb_milu_n_); the
  !               ILU(k) factorization otherwise.
  !    i       -  integer, input.
  !               The local index of the row to be copied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix under factorization.
  !    row     -  real(psb_dpk_), dimension(:), input/output.
  !               It contains, input, the row to be copied, and, in output,
  !               the null vector (the latter is used in the next call to
  !               iluk_copyin in psb_iluk_fact).
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) contains the fill kevel of the k-th entry
  !               of the row to be copied. rowlevs(k) = -(m+1) indicates that
  !               this entry is zero; however, any rowlevs(k) = -(m+1) is not
  !               used by the routine. In output rowlevs(k) = -(m+1) for all k's
  !               (this is an inizialization for the next call to iluk_copyin
  !               in psb_iluk_factint).
  !    nidx    -  integer, input.
  !               The number of entries of the array row that have been examined
  !               during the elimination step carried out by the routine iluk_fact.
  !    idxs    -  integer, dimension(:), allocatable, input.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step carried out by the routine
  !               iluk_fact.
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of lval.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uval.
  !    lja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               copied in lval row by row (see psb_diluk_factint), according
  !               to the CSR storage format.
  !    lirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor, copied in lval row by row (see
  !               psb_diluk_factint), according to the CSR storage format.
  !    lval   -  real(psb_dpk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               L factor are copied.
  !    d       -  real(psb_dpk_), dimension(:), input/output.
  !               The array where the inverse of the diagonal entry of the
  !               row is copied (only d(i) is used by the routine).
  !    uja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor
  !               copied in uval row by row (see psb_diluk_factint), according
  !               to the CSR storage format.
  !    uirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor copied in uval row by row (see
  !               psb_zilu_fctint), according to the CSR storage format.
  !    uval   -  real(psb_dpk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               U factor are copied.
  !    uplevs  -  integer, dimension(:), input.
  !               The fill levels of the nonzero entries in the part of the
  !               U factor above the current row.
  !
  subroutine iluk_copyout(fill_in,ialg,i,m,row,rowlevs,nidx,idxs,&
       &  l1,l2,lja,lirp,lval,d,uja,uirp,uval,uplevs,info)

    use psb_base_mod

    implicit none

    ! Arguments
    integer(psb_ipk_), intent(in)                  :: fill_in, ialg, i, m, nidx
    integer(psb_ipk_), intent(inout)               :: l1, l2, info
    integer(psb_ipk_), intent(inout)               :: rowlevs(:), idxs(:)
    integer(psb_ipk_), allocatable, intent(inout)  :: uja(:), uirp(:), lja(:), lirp(:),uplevs(:)
    real(psb_dpk_), allocatable, intent(inout)  :: uval(:), lval(:)
    real(psb_dpk_), intent(inout)     :: row(:), d(:)

    ! Local variables
    integer(psb_ipk_)              :: j,isz,err_act,int_err(5),idxp
    character(len=20), parameter  :: name='psb_diluk_factint'
    character(len=20)             :: ch_err

    info = psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    d(i) = dzero

    do idxp=1,nidx

      j = idxs(idxp)

      if (j<i) then
        !
        ! Copy the lower part of the row
        !
        if (rowlevs(j) <= fill_in) then
          l1     = l1 + 1
          if (size(lval) < l1) then
            !
            ! Figure out a good reallocation size!
            !
            isz  = (max((l1/i)*m,int(1.2*l1),l1+100))
            call psb_realloc(isz,lval,info)
            if (info == psb_success_) call psb_realloc(isz,lja,info)
            if (info /= psb_success_) then
              info=psb_err_from_subroutine_
              call psb_errpush(info,name,a_err='Allocate')
              goto 9999
            end if
          end if
          lja(l1)   = j
          lval(l1)  = row(j)
        else if (ialg == psb_milu_n_) then
          !
          ! MILU(k): add discarded entries to the diagonal one
          !
          d(i) = d(i) + row(j)
        end if
        !
        ! Re-initialize row(j) and rowlevs(j)
        !
        row(j)     = dzero
        rowlevs(j) = -(m+1)

      else if (j == i) then
        !
        ! Copy the diagonal entry of the row and re-initialize
        ! row(j) and rowlevs(j)
        !
        d(i)       = d(i) + row(i)
        row(i)     = dzero
        rowlevs(i) = -(m+1)

      else if (j>i) then
        !
        ! Copy the upper part of the row
        !
        if (rowlevs(j) <= fill_in) then
          l2     = l2 + 1
          if (size(uval) < l2) then
            !
            ! Figure out a good reallocation size!
            !
            isz  = max((l2/i)*m,int(1.2*l2),l2+100)
            call psb_realloc(isz,uval,info)
            if (info == psb_success_) call psb_realloc(isz,uja,info)
            if (info == psb_success_) call psb_realloc(isz,uplevs,info,pad=(m+1))
            if (info /= psb_success_) then
              info=psb_err_from_subroutine_
              call psb_errpush(info,name,a_err='Allocate')
              goto 9999
            end if
          end if
          uja(l2)   = j
          uval(l2)  = row(j)
          uplevs(l2) = rowlevs(j)
        else if (ialg == psb_milu_n_) then
          !
          ! MILU(k): add discarded entries to the diagonal one
          !
          d(i) = d(i) + row(j)
        end if
        !
        ! Re-initialize row(j) and rowlevs(j)
        !
        row(j)     = dzero
        rowlevs(j) = -(m+1)
      end if
    end do

    !
    ! Store the pointers to the first non occupied entry of in
    ! lval  and uval
    !
    lirp(i+1) = l1 + 1
    uirp(i+1) = l2 + 1

    !
    ! Check the pivot size
    !
    if (abs(d(i)) < d_epstol) then
      !
      ! Too small pivot: unstable factorization
      !
      info = psb_err_pivot_too_small_
      int_err(1) = i
      write(ch_err,'(g20.10)') d(i)
      call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
      goto 9999
    else
      !
      ! Compute 1/pivot
      !
      d(i) = done/d(i)
    end if

    !
    ! Scale the upper part
    !
    do j=uirp(i), uirp(i+1)-1
      uval(j) = d(i)*uval(j)
    end do

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine iluk_copyout


end subroutine psb_diluk_fact
