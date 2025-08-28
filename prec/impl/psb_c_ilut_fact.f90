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
! File: psb_cilut_fact.f90
!
! Subroutine: psb_cilut_fact
! Version:    complex
! Contains:   psb_cilut_factint, ilut_copyin, ilut_fact, ilut_copyout
!
!  This routine computes the ILU(k,t) factorization of the diagonal blocks
!  of a distributed matrix. This factorization is used to build the 'base
!  preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
!  preconditioner) corresponding to a certain level of a multilevel preconditioner.
!
!  Details on the above factorization can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix is stored into a and blck, as specified in the description
!  of the arguments below. The storage format for both the L and U factors is
!  CSR. The diagonal of the U factor is stored separately (actually, the
!  inverse of the diagonal entries is stored; this is then managed in the
!  solve stage associated to the ILU(k,t) factorization).
!
!
! Arguments:
!    fill_in -  integer, input.
!               The fill-in parameter k in ILU(k,t).
!    thres   -  real, input.
!               The threshold t, i.e. the drop tolerance, in ILU(k,t).
!    a       -  type(psb_cspmat_type), input.
!               The sparse matrix structure containing the local matrix.
!               Note that if the 'base' Additive Schwarz preconditioner
!               has overlap greater than 0 and the matrix has not been reordered
!               (see psb_fact_bld), then a contains only the 'original' local part
!               of the distributed matrix, i.e. the rows of the matrix held
!               by the calling process according to the initial data distribution.
!    l       -  type(psb_cspmat_type), input/output.
!               The L factor in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    u       -  type(psb_cspmat_type), input/output.
!               The U factor (except its diagonal) in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    d       -  complex(psb_spk_), dimension(:), input/output.
!               The inverse of the diagonal entries of the U factor in the incomplete
!               factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    info    -  integer, output.
!               Error code.
!    blck    -  type(psb_cspmat_type), input, optional, target.
!               The sparse matrix structure containing the remote rows of the
!               distributed matrix, that have been retrieved by psb_as_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 or the matrix has been reordered
!               (see psb_fact_bld), then blck does not contain any row.
!
subroutine psb_cilut_fact(fill_in,thres,a,l,u,d,info,blck,iscale,shft)

  use psb_base_mod
  use psb_c_ilu_fact_mod, psb_protect_name => psb_cilut_fact

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)         :: fill_in
  real(psb_spk_), intent(in)             :: thres
  integer(psb_ipk_), intent(out)        :: info
  type(psb_cspmat_type),intent(in)    :: a
  type(psb_cspmat_type),intent(inout) :: l,u
  complex(psb_spk_), intent(inout)        ::  d(:)
  type(psb_cspmat_type),intent(in), optional, target :: blck
  integer(psb_ipk_), intent(in), optional       :: iscale
  complex(psb_spk_), intent(in), optional :: shft
  !     Local Variables
  integer(psb_ipk_)   ::  l1, l2, m, err_act, iscale_

  complex(psb_spk_) :: shft_
  type(psb_cspmat_type), pointer  :: blck_
  type(psb_c_csr_sparse_mat)       :: ll, uu
  real(psb_spk_)      :: scale
  character(len=20)   :: name, ch_err

  name='psb_cilut_fact'
  info = psb_success_
  call psb_erractionsave(err_act)

  if (fill_in < 0) then
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name, &
         & i_err=(/ione,fill_in,izero,izero,izero/))
    goto 9999
  end if
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
  if (present(iscale)) then
    iscale_ = iscale
  else
    iscale_ = psb_ilu_scale_none_
  end if
  if (present(shft)) then
    shft_ = shft
  else
    shft_ = czero
  end if

  select case(iscale_)
  case(psb_ilu_scale_none_)
    scale = sone
  case(psb_ilu_scale_maxval_)
    scale = max(a%maxval(),blck_%maxval())
    scale = sone/scale
  case default
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name,i_err=(/ione*9,iscale_,izero,izero,izero/))
    goto 9999
  end select

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
  ! Compute the ILU(k,t) factorization
  !
  call psb_cilut_factint(fill_in,thres,a,blck_,&
       & d,ll%val,ll%ja,ll%irp,uu%val,uu%ja,uu%irp,l1,l2,info,scale,shft_)
  if (info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_cilut_factint'
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
  ! Subroutine: psb_cilut_factint
  ! Version:    complex
  ! Note: internal subroutine of psb_cilut_fact
  !
  !  This routine computes the ILU(k,t) factorization of the diagonal blocks of a
  !  distributed matrix. This factorization is used to build the 'base
  !  preconditioner' (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a certain level of a multilevel preconditioner.
  !
  !  The local matrix to be factorized is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(k,t) factorization).
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in parameter k in ILU(k,t).
  !    thres   -  real, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    m       -  integer, output.
  !               The total number of rows of the local matrix to be factorized,
  !               i.e. ma+mb.
  !    a       -  type(psb_cspmat_type), input.
  !               The sparse matrix structure containing the local matrix.
  !               Note that, if the 'base' Additive Schwarz preconditioner
  !               has overlap greater than 0 and the matrix has not been reordered
  !               (see psb_fact_bld), then a contains only the 'original' local part
  !               of the distributed matrix, i.e. the rows of the matrix held
  !               by the calling process according to the initial data distribution.
  !    b       -  type(psb_cspmat_type), input.
  !               The sparse matrix structure containing the remote rows of the
  !               distributed matrix, that have been retrieved by psb_as_bld
  !               to build an Additive Schwarz base preconditioner with overlap
  !               greater than 0. If the overlap is 0 or the matrix has been reordered
  !               (see psb_fact_bld), then b does not contain any row.
  !    d       -  complex(psb_spk_), dimension(:), output.
  !               The inverse of the diagonal entries of the U factor in the incomplete
  !               factorization.
  !    lval   -  complex(psb_spk_), dimension(:), input/output.
  !               The L factor in the incomplete factorization.
  !    lia1    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               according to the CSR storage format.
  !    lirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor in lval, according to the CSR storage format.
  !    uval   -  complex(psb_spk_), dimension(:), input/output.
  !               The U factor in the incomplete factorization.
  !               The entries of U are stored according to the CSR format.
  !    uja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor,
  !               according to the CSR storage format.
  !    uirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor in uval, according to the CSR storage format.
  !    l1      -  integer, output
  !               The number of nonzero entries in lval.
  !    l2      -  integer, output
  !               The number of nonzero entries in uval.
  !    info    -  integer, output.
  !               Error code.
  !
  subroutine psb_cilut_factint(fill_in,thres,a,b,&
       & d,lval,lja,lirp,uval,uja,uirp,l1,l2,info,scale,shft)

    use psb_base_mod

    implicit none

  ! Arguments
    integer(psb_ipk_), intent(in)                 :: fill_in
    real(psb_spk_), intent(in)                     :: thres
    type(psb_cspmat_type),intent(in)            :: a,b
    integer(psb_ipk_),intent(inout)               :: l1,l2,info
    integer(psb_ipk_), allocatable, intent(inout) :: lja(:),lirp(:),uja(:),uirp(:)
    complex(psb_spk_), allocatable, intent(inout)   :: lval(:),uval(:)
    complex(psb_spk_), intent(inout)                :: d(:)
    real(psb_spk_), intent(in), optional          :: scale
    complex(psb_spk_), intent(in)       :: shft

    ! Local Variables
    integer(psb_ipk_) :: i, ktrw,err_act,nidx,nlw,nup,jmaxup, ma, mb, m
    real(psb_spk_)               ::  nrmi
    real(psb_spk_)               ::  weight
    integer(psb_ipk_), allocatable         :: idxs(:)
    complex(psb_spk_), allocatable  :: row(:)
    type(psb_i_heap)             :: heap
    type(psb_c_coo_sparse_mat)   :: trw
    character(len=20), parameter :: name='psb_cilut_factint'
    character(len=20)            :: ch_err

    info = psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if


    ma = a%get_nrows()
    mb = b%get_nrows()
    m  = ma+mb

    !
    ! Allocate a temporary buffer for the ilut_copyin function
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
    ! Allocate memory to hold the entries of a row
    !
    allocate(row(m),stat=info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

    row(:) = czero
    weight = sone
    if (present(scale)) weight = abs(scale)
    !
    ! Cycle over the matrix rows
    !
    do i = 1, m

      !
      ! At each iteration of the loop we keep in a heap the column indices
      ! affected by the factorization. The heap is initialized and filled
      ! in the ilut_copyin function, and updated during the elimination, in
      ! the ilut_fact routine. The heap is ideal because at each step we need
      ! the lowest index, but we also need to insert new items, and the heap
      ! allows to do both in log time.
      !
      d(i) = czero
      if (i<=ma) then
        call ilut_copyin(i,ma,a,i,ione,m,nlw,nup,jmaxup,nrmi,weight,&
             & row,heap,ktrw,trw,info,shft)
      else
        call ilut_copyin(i-ma,mb,b,i,ione,m,nlw,nup,jmaxup,nrmi,weight,&
             & row,heap,ktrw,trw,info,shft)
      endif

      !
      ! Do an elimination step on current row
      !
      if (info == psb_success_) call ilut_fact(thres,i,nrmi,row,heap,&
           & d,uja,uirp,uval,nidx,idxs,info)
      !
      ! Copy the row into lval/d(i)/uval
      !
      if (info == psb_success_) call ilut_copyout(fill_in,thres,i,m,&
           & nlw,nup,jmaxup,nrmi,row,nidx,idxs,&
           & l1,l2,lja,lirp,lval,d,uja,uirp,uval,info)

      if (info /= psb_success_) then
        info=psb_err_internal_error_
        call psb_errpush(info,name,a_err='Copy/factor loop')
        goto 9999
      end if

    end do
    !
    ! Adjust diagonal accounting for scale factor
    !
    if (weight /= sone) then
      d(1:m) = d(1:m)*weight
    end if

    !
    ! And we're sone, so deallocate the memory
    !
    deallocate(row,idxs,stat=info)
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

  end subroutine psb_cilut_factint

  !
  ! Subroutine: ilut_copyin
  ! Version:    complex
  ! Note: internal subroutine of psb_cilut_fact
  !
  !  This routine performs the following tasks:
  !  - copying a row of a sparse matrix A, stored in the sparse matrix structure a,
  !    into the array row;
  !  - storing into a heap the column indices of the nonzero entries of the copied
  !    row;
  !  - computing the column index of the first entry with maximum absolute value
  !    in the part of the row belonging to the upper triangle;
  !  - computing the 2-norm of the row.
  !  The output array row is such that it contains a full row of A, i.e. it contains
  !  also the zero entries of the row. This is useful for the elimination step
  !  performed by ilut_fact after the call to ilut_copyin (see psb_ilut_factint).
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied, row by row, into the array row, through successive calls to
  !  ilut_copyin.
  !
  !  This routine is used by psb_cilut_factint in the computation of the ILU(k,t)
  !  factorization of a local sparse matrix.
  !
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from the
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(psb_cspmat_type), input.
  !               The sparse matrix structure containing the row to be
  !               copied.
  !    jd      -  integer, input.
  !               The column index of the diagonal entry of the row to be
  !               copied.
  !    jmin    -  integer, input.
  !               The minimum valid column index.
  !    jmax    -  integer, input.
  !               The maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    nlw     -  integer, output.
  !               The number of nonzero entries in the part of the row
  !               belonging to the lower triangle of the matrix.
  !    nup     -  integer, output.
  !               The number of nonzero entries in the part of the row
  !               belonging to the upper triangle of the matrix.
  !    jmaxup  -  integer, output.
  !               The column index of the first entry with maximum absolute
  !               value in the part of the row belonging to the upper triangle
  !    nrmi    -  real(psb_spk_), output.
  !               The 2-norm of the current row.
  !    row     -  complex(psb_spk_), dimension(:), input/output.
  !               In input it is the null vector (see psb_ilut_factint and
  !               ilut_copyout). In output it contains the row extracted
  !               from the matrix A. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    rowlevs -  integer, dimension(:), input/output.
  !               In input rowlevs(k) = -(m+1) for k=1,...,m. In output
  !               rowlevs(k) = 0 for 1 <= k <= jmax and A(i,k) /= 0, for
  !               future use in ilut_fact.
  !    heap    -  type(psb_int_heap), input/output.
  !               The heap containing the column indices of the nonzero
  !               entries in the array row.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, sone by psb_init_heap inside this
  !               routine.
  !    ktrw    -  integer, input/output.
  !               The index identifying the last entry taken from the
  !               staging buffer trw. See below.
  !    trw     -  type(psb_cspmat_type), input/output.
  !               A staging buffer. If the matrix A is not in CSR format, we use
  !               the psb_sp_getblk routine and store its output in trw; when we
  !               need to call psb_sp_getblk we do it for a block of rows, and then
  !               we consume them from trw in successive calls to this routine,
  !               until we empty the buffer. Thus we will make a call to psb_sp_getblk
  !               every nrb calls to copyin. If A is in CSR format it is unused.
  !
  subroutine ilut_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,&
       & nrmi,weight,row,heap,ktrw,trw,info,shft)
    use psb_base_mod
    implicit none
    type(psb_cspmat_type), intent(in)         :: a
    type(psb_c_coo_sparse_mat), intent(inout) :: trw
    integer(psb_ipk_), intent(in)               :: i, m,jmin,jmax,jd
    integer(psb_ipk_), intent(inout)            :: ktrw,nlw,nup,jmaxup,info
    real(psb_spk_), intent(inout)                :: nrmi
    complex(psb_spk_), intent(inout)              :: row(:)
    real(psb_spk_), intent(in)                 :: weight
    type(psb_i_heap), intent(inout)             :: heap
    complex(psb_spk_), intent(in)       :: shft

    integer(psb_ipk_)               :: k,j,irb,kin,nz
    integer(psb_ipk_), parameter    :: nrb=40
    real(psb_spk_)        :: dmaxup
    real(psb_spk_), external    :: dnrm2
    character(len=20), parameter  :: name='psb_cilut_factint'

    info = psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    call heap%init(info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_init_heap')
      goto 9999
    end if

    !
    ! nrmi is the norm of the current sparse row (for the time being,
    ! we use the 2-norm).
    ! NOTE: the 2-norm below includes also elements that are outside
    ! [jmin:jmax] strictly. Is this really important? TO BE CHECKED.
    !

    nlw    = 0
    nup    = 0
    jmaxup = 0
    dmaxup = szero
    nrmi   = szero

    select type (aa=> a%a)
    type is (psb_c_csr_sparse_mat)
      !
      ! Take a fast shortcut if the matrix is stored in CSR format
      !

      do j = aa%irp(i), aa%irp(i+1) - 1
        k          = aa%ja(j)
        if ((jmin<=k).and.(k<=jmax)) then
          row(k)     = aa%val(j)*weight
          call heap%insert(k,info)
          if (info /= psb_success_) exit
          if (k<jd) nlw = nlw + 1
          if (k == jd) row(k) = row(k) + (shft*weight)
          if (k>jd) then
            nup = nup + 1
            if (abs(row(k))>dmaxup) then
              jmaxup = k
              dmaxup = abs(row(k))
            end if
          end if
        end if
      end do
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if

      nz   = aa%irp(i+1) - aa%irp(i)
      nrmi = weight*dnrm2(nz,aa%val(aa%irp(i)),ione)


    class default

      !
      ! Otherwise use psb_sp_getblk, slower but able (in principle) of
      ! handling any format. In this case, a block of rows is extracted
      ! instead of a single row, for performance reasons, and these
      ! rows are copied one by one into the array row, through successive
      ! calls to ilut_copyin.
      !

      if ((mod(i,nrb) == 1).or.(nrb == 1)) then
        irb = min(m-i+1,nrb)
        call aa%csget(i,i+irb-1,trw,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_sp_getblk')
          goto 9999
        end if
        ktrw=1
      end if

      kin = ktrw
      nz = trw%get_nzeros()
      do
        if (ktrw > nz) exit
        if (trw%ia(ktrw) > i) exit
        k          = trw%ja(ktrw)
        if ((jmin<=k).and.(k<=jmax)) then
          row(k)     = trw%val(ktrw)*weight
          call heap%insert(k,info)
          if (info /= psb_success_) exit
          if (k<jd) nlw = nlw + 1
          if (k == jd) row(k) = row(k) + (shft*weight)
          if (k>jd) then
            nup = nup + 1
            if (abs(row(k))>dmaxup) then
              jmaxup = k
              dmaxup = abs(row(k))
            end if
          end if
        end if
        ktrw       = ktrw + 1
      enddo
      nz = ktrw - kin
      nrmi = weight*dnrm2(nz,trw%val(kin),ione)
    end select

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine ilut_copyin

  !
  ! Subroutine: ilut_fact
  ! Version:    complex
  ! Note: internal subroutine of psb_cilut_fact
  !
  !  This routine does an elimination step of the ILU(k,t) factorization on a single
  !  matrix row (see the calling routine psb_ilut_factint). Actually, only the dropping
  !  rule based on the threshold is applied here. The dropping rule based on the
  !  fill-in is applied by ilut_copyout.
  !
  !  The routine is used by psb_cilut_factint in the computation of the ILU(k,t)
  !  factorization of a local sparse matrix.
  !
  !
  ! Arguments
  !    thres   -  real, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    i       -  integer, input.
  !               The local index of the row to which the factorization is applied.
  !    nrmi    -  real(psb_spk_), input.
  !               The 2-norm of the row to which the elimination step has to be
  !               applied.
  !    row     -  complex(psb_spk_), dimension(:), input/output.
  !               In input it contains the row to which the elimination step
  !               has to be applied. In output it contains the row after the
  !               elimination step. It actually contains a full row, i.e.
  !               it contains also the zero entries of the row.
  !    heap    -  type(psb_i_heap), input/output.
  !               The heap containing the column indices of the nonzero entries
  !               in the processed row. In input it contains the indices concerning
  !               the row before the elimination step, while in output it contains
  !               the previous indices plus the ones corresponding to transformed
  !               entries in the 'upper part' that have not been dropped.
  !    d       -  complex(psb_spk_), input.
  !               The inverse of the diagonal entries of the part of the U factor
  !               above the current row (see ilut_copyout).
  !    uja    -  integer, dimension(:), input.
  !               The column indices of the nonzero entries of the part of the U
  !               factor above the current row, stored in uval row by row (see
  !               ilut_copyout, called by psb_cilut_factint), according to the CSR
  !               storage format.
  !    uirp    -  integer, dimension(:), input.
  !               The indices identifying the first nonzero entry of each row of
  !               the U factor above the current row, stored in uval row by row
  !               (see ilut_copyout, called by psb_cilut_factint), according to
  !               the CSR storage format.
  !    uval   -  complex(psb_spk_), dimension(:), input.
  !               The entries of the U factor above the current row (except the
  !               diagonal ones), stored according to the CSR format.
  !    nidx    -  integer, output.
  !               The number of entries of the array row that have been
  !               examined during the elimination step. This will be used
  !               by the routine ilut_copyout.
  !    idxs    -  integer, dimension(:), allocatable, input/output.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step.This will be used by
  !               by the routine ilut_copyout.
  !               Note: this argument is intent(inout) and not only intent(out)
  !               to retain its allocation, sone by this routine.
  !
  subroutine ilut_fact(thres,i,nrmi,row,heap,d,uja,uirp,uval,nidx,idxs,info)

    use psb_base_mod

    implicit none

  ! Arguments
    type(psb_i_heap), intent(inout)               :: heap
    integer(psb_ipk_), intent(in)                 :: i
    integer(psb_ipk_), intent(inout)              :: nidx,info
    real(psb_spk_), intent(in)                     :: thres,nrmi
    integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
    integer(psb_ipk_), intent(inout)              :: uja(:),uirp(:)
    complex(psb_spk_), intent(inout)       :: row(:), uval(:),d(:)

    ! Local Variables
    integer(psb_ipk_)               :: k,j,jj,lastk,iret
    complex(psb_spk_)      :: rwk

    info  = psb_success_
    call psb_ensure_size(200*ione,idxs,info)
    if (info /= psb_success_) return
    nidx  = 0
    lastk = -1
    !
    ! Do while there are indices to be processed
    !
    do

      call heap%get_first(k,iret)
      if (iret < 0) exit

      !
      ! An index may have been put on the heap more than once.
      !
      if (k == lastk) cycle

      lastk = k
      lowert: if (k<i)  then

        !
        ! Dropping rule based on the threshold: compare the absolute
        ! value of each updated entry of row with thres * 2-norm of row.
        !
        rwk    = row(k)
        row(k) = row(k) * d(k)
        if (abs(row(k)) < thres*nrmi) then
          !
          ! Drop the entry.
          !
          row(k) = czero
          cycle
        else
          !
          ! Note: since U is scaled while copying it out (see ilut_copyout),
          ! we can use rwk in the update below.
          !
          do jj=uirp(k),uirp(k+1)-1
            j = uja(jj)
            if (j<=k) then
              info = -i
              return
            endif
            !
            ! Update row(j) and, if it is not to be discarded, insert
            ! its index into the heap for further processing.
            !
            row(j)     = row(j) - rwk * uval(jj)
            if (abs(row(j)) < thres*nrmi) then
              !
              ! Drop the entry.
              !
              row(j) = czero
            else
              !
              ! Do the insertion.
              !
              call heap%insert(j,info)
              if (info /= psb_success_) return
            endif
          end do
        end if
      end if lowert

      !
      ! If we get here it is an index we need to keep on copyout.
      !
      nidx       = nidx + 1
      call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)
      if (info /= psb_success_) return
      idxs(nidx) = k

    end do

  end subroutine ilut_fact

  !
  ! Subroutine: ilut_copyout
  ! Version:    complex
  ! Note: internal subroutine of psb_cilut_fact
  !
  !  This routine copies a matrix row, computed by ilut_fact by applying an
  !  elimination step of the ILU(k,t) factorization, into the arrays lval,
  !  uval, d, corresponding to the L factor, the U factor and the diagonal
  !  of U, respectively.
  !
  !  Note that
  !  - the dropping rule based on the fill-in is applied here and not in ilut_fact;
  !    it consists in keeping the nlw+k entries with largest absolute value in
  !    the 'lower part' of the row, and the nup+k ones in the 'upper part';
  !  - the entry in the upper part of the row which has maximum absolute value
  !    in the original matrix is included in the above nup+k entries anyway;
  !  - the part of the row stored into uval is scaled by the corresponding
  !    diagonal entry, according to the LDU form of the incomplete factorization;
  !  - the inverse of the diagonal entries of U is actually stored into d; this
  !    is then managed in the solve stage associated to the ILU(k,t) factorization;
  !  - the row entries are stored in lval and uval according to the CSR format;
  !  - the array row is re-initialized for future use in psb_ilut_fact(see also
  !    ilut_copyin and ilut_fact).
  !
  !  This routine is used by psb_cilut_factint in the computation of the ILU(k,t)
  !  factorization of a local sparse matrix.
  !
  !
  ! Arguments:
  !    fill_in -  integer, input.
  !               The fill-in parameter k in ILU(k,t).
  !    thres   -  real, input.
  !               The threshold t, i.e. the drop tolerance, in ILU(k,t).
  !    i       -  integer, input.
  !               The local index of the row to be copied.
  !    m       -  integer, input.
  !               The number of rows of the local matrix under factorization.
  !    nlw     -  integer, input.
  !               The number of nonzero entries of the 'lower part' of the row
  !               in the initial matrix (i.e. the matrix before the factorization).
  !    nup     -  integer, input.
  !               The number of nonzero entries in the 'upper part' of the row
  !               in the initial matrix.
  !    jmaxup  -  integer, input.
  !               The column index of the first entry with maximum absolute
  !               value in the 'upper part' of the row in the initial matrix.
  !    nrmi    -  real(psb_spk_), input.
  !               The 2-norm of the current row in the initial matrix.
  !    row     -  complex(psb_spk_), dimension(:), input/output.
  !               It contains, input, the row to be copied, and, in output,
  !               the null vector (the latter is used in the next call to
  !               ilut_copyin in psb_ilut_fact).
  !    nidx    -  integer, input.
  !               The number of entries of the array row that have been examined
  !               during the elimination step carried out by the routine ilut_fact.
  !    idxs    -  integer, dimension(:), allocatable, input.
  !               The indices of the entries of the array row that have been
  !               examined during the elimination step carried out by the routine
  !               ilut_fact.
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of lval.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uval.
  !    lja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               copied in lval row by row (see psb_cilut_factint), according
  !               to the CSR storage format.
  !    lirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor, copied in lval row by row (see
  !               psb_cilut_factint), according to the CSR storage format.
  !    lval   -  complex(psb_spk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               L factor are copied.
  !    d       -  complex(psb_spk_), dimension(:), input/output.
  !               The array where the inverse of the diagonal entry of the
  !               row is copied (only d(i) is used by the routine).
  !    uja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor
  !               copied in uval row by row (see psb_cilut_factint), according
  !               to the CSR storage format.
  !    uirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor copied in uval row by row (see
  !               psb_dilu_fctint), according to the CSR storage format.
  !    uval   -  complex(psb_spk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               U factor are copied.
  !
  subroutine ilut_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,&
       & nrmi,row, nidx,idxs,l1,l2,lja,lirp,lval,&
       & d,uja,uirp,uval,info)

    use psb_base_mod

    implicit none

    ! Arguments
    integer(psb_ipk_), intent(in)                       :: fill_in,i,m,nidx,nlw,nup,jmaxup
    integer(psb_ipk_), intent(in)                       :: idxs(:)
    integer(psb_ipk_), intent(inout)                    :: l1,l2, info
    integer(psb_ipk_), allocatable, intent(inout)       :: uja(:),uirp(:), lja(:),lirp(:)
    real(psb_spk_), intent(in)                :: thres,nrmi
    complex(psb_spk_),allocatable, intent(inout) :: uval(:), lval(:)
    complex(psb_spk_), intent(inout)             :: row(:), d(:)

    ! Local variables
    complex(psb_spk_),allocatable :: xw(:)
    integer(psb_ipk_), allocatable          :: xwid(:), indx(:)
    complex(psb_spk_)             :: witem
    integer(psb_ipk_)                       :: widx
    integer(psb_ipk_)                       :: k,isz,err_act,int_err(5),idxp, nz
    type(psb_c_idx_heap)        :: heap
    character(len=20), parameter  :: name='ilut_copyout'
    character(len=20)             :: ch_err
    logical                       :: fndmaxup

    info=psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    !
    ! Here we need to apply also the dropping rule base on the fill-in.
    ! We do it by putting into a heap the elements that are not dropped
    ! by using the 2-norm rule, and then copying them out.
    !
    ! The heap goes down on the entry absolute value, so the first item
    ! is the largest absolute value.
    !

    call heap%init(info,dir=psb_asort_down_)

    if (info == psb_success_) allocate(xwid(nidx),xw(nidx),indx(nidx),stat=info)
    if (info /= psb_success_) then
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/3*nidx,izero,izero,izero,izero/),&
           & a_err='complex(psb_spk_)')
      goto 9999
    end if

    !
    ! First the lower part
    !

    nz   = 0
    idxp = 0

    do

      idxp = idxp + 1
      if (idxp > nidx) exit
      if (idxs(idxp) >= i) exit
      widx      = idxs(idxp)
      witem     = row(widx)
      !
      ! Dropping rule based on the 2-norm
      !
      if (abs(witem) < thres*nrmi) cycle

      nz       = nz + 1
      xw(nz)   = witem
      xwid(nz) = widx
      call heap%insert(witem,widx,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if
    end do

    !
    ! Now we have to take out the first nlw+fill_in entries
    !
    if (nz <= nlw+fill_in) then
      !
      ! Just copy everything from xw, and it is already ordered
      !
    else
      nz = nlw+fill_in
      do k=1,nz
        call heap%get_first(witem,widx,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='psb_heap_get_first')
          goto 9999
        end if

        xw(k)   = witem
        xwid(k) = widx
      end do
    end if

    !
    ! Now put things back into ascending column order
    !
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)

    !
    ! Copy out the lower part of the row
    !
    do k=1,nz
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
      lja(l1)   = xwid(k)
      lval(l1)  = xw(indx(k))
    end do

    !
    ! Make sure idxp points to the diagonal entry
    !
    if (idxp <= size(idxs)) then
      if (idxs(idxp) < i) then
        do
          idxp = idxp + 1
          if (idxp > nidx) exit
          if (idxs(idxp) >= i) exit
        end do
      end if
    end if
    if (idxp > size(idxs)) then
!!$      write(0,*) 'Warning: missing diagonal element in the row '
    else
      if (idxs(idxp) > i) then
!!$        write(0,*) 'Warning: missing diagonal element in the row '
      else if (idxs(idxp) /= i) then
!!$        write(0,*) 'Warning: impossible error: diagonal has vanished'
      else
        !
        ! Copy the diagonal entry
        !
        widx      = idxs(idxp)
        witem     = row(widx)
        d(i)      = witem
        if (abs(d(i)) < s_epstol) then
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
          d(i) = cone/d(i)
        end if
      end if
    end if

    !
    ! Now the upper part
    !

    call heap%init(info,dir=psb_asort_down_)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_init_heap')
      goto 9999
    end if

    nz       = 0
    do

      idxp = idxp + 1
      if (idxp > nidx) exit
      widx      = idxs(idxp)
      if (widx <= i) then
!!$        write(0,*) 'Warning: lower triangle in upper copy',widx,i,idxp,idxs(idxp)
        cycle
      end if
      if (widx > m) then
!!$        write(0,*) 'Warning: impossible value',widx,i,idxp,idxs(idxp)
        cycle
      end if
      witem     = row(widx)
      !
      ! Dropping rule based on the 2-norm. But keep the jmaxup-th entry anyway.
      !
      if ((widx /= jmaxup) .and. (abs(witem) < thres*nrmi)) then
        cycle
      end if

      nz       = nz + 1
      xw(nz)   = witem
      xwid(nz) = widx
      call heap%insert(witem,widx,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if

    end do

    !
    ! Now we have to take out the first nup-fill_in entries. But make sure
    ! we include entry jmaxup.
    !
    if (nz <= nup+fill_in) then
      !
      ! Just copy everything from xw
      !
      fndmaxup=.true.
    else
      fndmaxup = .false.
      nz = nup+fill_in
      do k=1,nz
        call heap%get_first(witem,widx,info)
        xw(k)   = witem
        xwid(k) = widx
        if (widx == jmaxup) fndmaxup=.true.
      end do
    end if
    if ((i<jmaxup).and.(jmaxup<=m)) then
      if (.not.fndmaxup) then
        !
        ! Include entry jmaxup, if it is not already there.
        ! Put it in the place of the smallest coefficient.
        !
        xw(nz)   = row(jmaxup)
        xwid(nz) = jmaxup
      endif
    end if

    !
    ! Now we put things back into ascending column order
    !
    call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)

    !
    ! Copy out the upper part of the row
    !
    do k=1,nz
      l2     = l2 + 1
      if (size(uval) < l2) then
        !
        ! Figure out a good reallocation size!
        !
        isz  = max((l2/i)*m,int(1.2*l2),l2+100)
        call psb_realloc(isz,uval,info)
        if (info == psb_success_) call psb_realloc(isz,uja,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='Allocate')
          goto 9999
        end if
      end if
      uja(l2)   = xwid(k)
      uval(l2)  = d(i)*xw(indx(k))
    end do
    !
    ! Set row to zero
    !
    do idxp=1,nidx
      row(idxs(idxp)) = czero
    end do

    !
    ! Store the pointers to the first non occupied entry of in
    ! lval and uval
    !
    lirp(i+1) = l1 + 1
    uirp(i+1) = l2 + 1

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine ilut_copyout


end subroutine psb_cilut_fact
