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
! File: psb_zilu0_fact.f90
!
! Subroutine: psb_zilu0_fact
! Version:    complex
! Contains:   psb_zilu0_factint, ilu_copyin
!
!  This routine computes either the ILU(0) or the MILU(0) factorization of
!  the diagonal blocks of a distributed matrix. These factorizations are used
!  to build the 'base preconditioner' (block-Jacobi preconditioner/solver,
!  Additive Schwarz preconditioner) corresponding to a given level of a
!  multilevel preconditioner.
!
!  Details on the above factorizations can be found in
!    Y. Saad, Iterative Methods for Sparse Linear Systems, Second Edition,
!    SIAM, 2003, Chapter 10.
!
!  The local matrix is stored into a and blck, as specified in the description
!  of the arguments below. The storage format for both the L and U factors is CSR.
!  The diagonal of the U factor is stored separately (actually, the inverse of the
!  diagonal entries is stored; this is then managed in the solve stage associated
!  to the ILU(0)/MILU(0) factorization).
!
!  The routine copies and factors "on the fly" from a and blck into l (L factor),
!  u (U factor, except its diagonal) and d (diagonal of U).
!
!  This implementation of ILU(0)/MILU(0) is faster than the implementation in
!  psb_ziluk_fct (the latter routine performs the more general ILU(k)/MILU(k)).
!
!
! Arguments:
!    ialg    -  integer, input.
!               The type of incomplete factorization to be performed.
!               The MILU(0) factorization is computed if ialg = 2 (= psb_milu_n_);
!               the ILU(0) factorization otherwise.
!    a       -  type(psb_zspmat_type), input.
!               The sparse matrix structure containing the local matrix.
!               Note that if the 'base' Additive Schwarz preconditioner
!               has overlap greater than 0 and the matrix has not been reordered
!               (see psb_as_bld), then a contains only the 'original' local part
!               of the distributed matrix, i.e. the rows of the matrix held
!               by the calling process according to the initial data distribution.
!    l       -  type(psb_zspmat_type), input/output.
!               The L factor in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    u       -  type(psb_zspmat_type), input/output.
!               The U factor (except its diagonal) in the incomplete factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    d       -  complex(psb_dpk_), dimension(:), input/output.
!               The inverse of the diagonal entries of the U factor in the incomplete
!               factorization.
!               Note: its allocation is managed by the calling routine psb_ilu_bld,
!               hence it cannot be only intent(out).
!    info    -  integer, output.
!               Error code.
!    blck    -  type(psb_zspmat_type), input, optional, target.
!               The sparse matrix structure containing the remote rows of the
!               distributed matrix, that have been retrieved by psb_as_bld
!               to build an Additive Schwarz base preconditioner with overlap
!               greater than 0. If the overlap is 0 or the matrix has been reordered
!               (see psb_fact_bld), then blck is empty.
!
subroutine psb_zilu0_fact(ialg,a,l,u,d,info,blck, upd,shft)

  use psb_base_mod
  use psb_z_ilu_fact_mod, psb_protect_name => psb_zilu0_fact

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)       :: ialg
  type(psb_zspmat_type),intent(inout) :: a
  type(psb_zspmat_type),intent(inout) :: l,u
  complex(psb_dpk_), intent(inout)    :: d(:)
  integer(psb_ipk_), intent(out)      :: info
  type(psb_zspmat_type),intent(in), optional, target :: blck
  character, intent(in), optional     :: upd
  complex(psb_dpk_), intent(in), optional :: shft

  ! Local variables
  integer(psb_ipk_)   :: l1, l2, m, err_act
  type(psb_zspmat_type), pointer  :: blck_
  type(psb_z_csr_sparse_mat)      :: ll, uu, aa, bb
  complex(psb_dpk_) :: shft_
  character                       :: upd_
  character(len=20)    :: name, ch_err

  name='psb_zilu0_fact'
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
  if (present(upd)) then
    upd_ = psb_toupper(upd)
  else
    upd_ = 'F'
  end if
  if (present(shft)) then
    shft_ = shft
  else
    shft_ = zzero
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
  call a%cp_to(aa)
  call blck_%cp_to(bb)
  !
  ! Compute the ILU(0) or the MILU(0) factorization, depending on ialg
  !
  call psb_zilu0_factint(ialg,aa,bb,&
       & d,ll%val,ll%ja,ll%irp,uu%val,uu%ja,uu%irp,l1,l2,upd_,shft_,info)
  if(info.ne.0) then
    info=psb_err_from_subroutine_
    ch_err='psb_zilu0_factint'
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
    if(info.ne.0) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    deallocate(blck_)
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  !
  ! Subroutine: psb_zilu0_factint
  ! Version:    complex
  ! Note: internal subroutine of psb_zilu0_fact.
  !
  !  This routine computes either the ILU(0) or the MILU(0) factorization of the
  !  diagonal blocks of a distributed matrix.
  !  These factorizations are used to build the 'base preconditioner'
  !  (block-Jacobi preconditioner/solver, Additive Schwarz
  !  preconditioner) corresponding to a given level of a multilevel preconditioner.
  !
  !  The local matrix is stored into a and b, as specified in the
  !  description of the arguments below. The storage format for both the L and U
  !  factors is CSR. The diagonal of the U factor is stored separately (actually,
  !  the inverse of the diagonal entries is stored; this is then managed in the
  !  solve stage associated to the ILU(0)/MILU(0) factorization).
  !
  !  The routine copies and factors "on the fly" from the sparse matrix structures a
  !  and b into the arrays lval, uval, d (L, U without its diagonal, diagonal of U).
  !
  !
  ! Arguments:
  !    ialg    -  integer, input.
  !               The type of incomplete factorization to be performed.
  !               The ILU(0) factorization is computed if ialg = 1 (= psb_ilu_n_),
  !               the MILU(0) one if ialg = 2 (= psb_milu_n_); other values
  !               are not allowed.
  !    m       -  integer, output.
  !               The total number of rows of the local matrix to be factorized,
  !               i.e. ma+mb.
  !    ma      -  integer, input
  !               The number of rows of the local submatrix stored into a.
  !    a       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the local matrix.
  !               Note that, if the 'base' Additive Schwarz preconditioner
  !               has overlap greater than 0 and the matrix has not been reordered
  !               (see psb_fact_bld), then a contains only the 'original' local part
  !               of the distributed matrix, i.e. the rows of the matrix held
  !               by the calling process according to the initial data distribution.
  !    mb      -  integer, input.
  !               The number of rows of the local submatrix stored into b.
  !    b       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the remote rows of the
  !               distributed matrix, that have been retrieved by psb_as_bld
  !               to build an Additive Schwarz base preconditioner with overlap
  !               greater than 0. If the overlap is 0 or the matrix has been
  !               reordered (see psb_fact_bld), then b does not contain any row.
  !    d       -  complex(psb_dpk_), dimension(:), output.
  !               The inverse of the diagonal entries of the U factor in the
  !               incomplete factorization.
  !    lval   -   complex(psb_dpk_), dimension(:), input/output.
  !               The entries of U are stored according to the CSR format.
  !               The L factor in the incomplete factorization.
  !    lja    -   integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the L factor,
  !               according to the CSR storage format.
  !    lirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the L factor in lval, according to the CSR storage format.
  !    uval   -   complex(psb_dpk_), dimension(:), input/output.
  !               The U factor in the incomplete factorization.
  !               The entries of U are stored according to the CSR format.
  !    uja    -   integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the U factor,
  !               according to the CSR storage format.
  !    uirp    -  integer, dimension(:), input/output.
  !               The indices identifying the first nonzero entry of each row
  !               of the U factor in uval, according to the CSR storage format.
  !    l1      -  integer, output.
  !               The number of nonzero entries in lval.
  !    l2      -  integer, output.
  !               The number of nonzero entries in uval.
  !    info    -  integer, output.
  !               Error code.
  !
  subroutine psb_zilu0_factint(ialg,a,b,&
       & d,lval,lja,lirp,uval,uja,uirp,l1,l2,upd,shft,info)

    implicit none

    ! Arguments
    integer(psb_ipk_), intent(in)     :: ialg
    class(psb_z_base_sparse_mat),intent(inout) :: a,b
    integer(psb_ipk_),intent(inout)   :: l1,l2,info
    integer(psb_ipk_), intent(inout)  :: lja(:),lirp(:),uja(:),uirp(:)
    complex(psb_dpk_), intent(inout)  :: lval(:),uval(:),d(:)
    character, intent(in)             :: upd
    complex(psb_dpk_), intent(in)       :: shft

    ! Local variables
    integer(psb_ipk_) :: i,j,k,l,low1,low2,kk,jj,ll, ktrw,err_act, m
    integer(psb_ipk_) :: ma,mb
    complex(psb_dpk_) :: dia,temp
    integer(psb_ipk_), parameter :: nrb=16
    type(psb_z_coo_sparse_mat) :: trw
    integer(psb_ipk_)   :: int_err(5)
    character(len=20)   :: name, ch_err

    name='psb_zilu0_factint'
    info=psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if
    ma = a%get_nrows()
    mb = b%get_nrows()

    select case(ialg)
    case(psb_ilu_n_,psb_milu_n_)
      ! Ok
    case default
      info=psb_err_input_asize_invalid_i_
      call psb_errpush(info,name,&
           & i_err=(/ione,ialg,izero,izero,izero/))
      goto 9999
    end select

    call trw%allocate(izero,izero,ione)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    m = ma+mb

    if (psb_toupper(upd) == 'F' ) then
      lirp(1) = 1
      uirp(1) = 1
      l1      = 0
      l2      = 0

      !
      ! Cycle over the matrix rows
      !
      do i = 1, m

      d(i) = zzero

        if (i <= ma) then
          !
          ! Copy the i-th local row of the matrix, stored in a,
          ! into lval/d(i)/uval
          !
          call ilu_copyin(i,ma,a,i,ione,m,l1,lja,lval,&
               & d(i),l2,uja,uval,ktrw,trw,upd,shft_)
        else
          !
          ! Copy the i-th local row of the matrix, stored in b
          ! (as (i-ma)-th row), into lval/d(i)/uval
          !
          call ilu_copyin(i-ma,mb,b,i,ione,m,l1,lja,lval,&
               & d(i),l2,uja,uval,ktrw,trw,upd,shft_)
        endif

        lirp(i+1) = l1 + 1
        uirp(i+1) = l2 + 1

        dia = d(i)
        do kk = lirp(i), lirp(i+1) - 1
          !
          ! Compute entry l(i,k) (lower factor L) of the incomplete
          ! factorization
          !
          temp      = lval(kk)
          k         = lja(kk)
          lval(kk) = temp*d(k)
          !
          ! Update the rest of row i (lower and upper factors L and U)
          ! using l(i,k)
          !
          low1 = kk + 1
          low2 = uirp(i)
          !
          updateloop: do  jj = uirp(k), uirp(k+1) - 1
            !
            j = uja(jj)
            !
            if (j < i) then
              !
              ! search l(i,*) (i-th row of L) for a matching index j
              !
              do  ll = low1, lirp(i+1) - 1
                l = lja(ll)
                if (l > j) then
                  low1 = ll
                  exit
                else if (l == j) then
                  lval(ll) = lval(ll) - temp*uval(jj)
                  low1 = ll + 1
                  cycle updateloop
                end if
              enddo

            else if (j == i) then
              !
              ! j=i: update the diagonal
              !
              dia = dia - temp*uval(jj)
              cycle updateloop
              !
            else if (j > i) then
              !
              ! search u(i,*) (i-th row of U) for a matching index j
              !
              do ll = low2, uirp(i+1) - 1
                l = uja(ll)
                if (l > j) then
                  low2 = ll
                  exit
                else if (l == j) then
                  uval(ll) = uval(ll) - temp*uval(jj)
                  low2 = ll + 1
                  cycle updateloop
                end if
              enddo
            end if
            !
            ! If we get here we missed the cycle updateloop, which means
            ! that this entry does not match; thus we accumulate on the
            ! diagonal for MILU(0).
            !
            if (ialg == psb_milu_n_) then
              dia = dia - temp*uval(jj)
            end if
          enddo updateloop
        enddo
        !
        ! Check the pivot size
        !
        if (abs(dia) < d_epstol) then
          !
          ! Too small pivot: unstable factorization
          !
          info = psb_err_pivot_too_small_
          int_err(1) = i
          write(ch_err,'(g20.10)') abs(dia)
          call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
          goto 9999
        else
          !
          ! Compute 1/pivot
          !
          dia = zone/dia
        end if
        d(i) = dia
        !
        ! Scale row i of upper triangle
        !
        do  kk = uirp(i), uirp(i+1) - 1
          uval(kk) = uval(kk)*dia
        enddo
      enddo
    else
      write(0,*) 'Update not implemented '
      info = 31
      call psb_errpush(info,name,&
           & i_err=(/ione*13,izero,izero,izero,izero/),a_err=upd)
      goto 9999

    end if

    call trw%free()

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_zilu0_factint

  !
  ! Subroutine: ilu_copyin
  ! Version:    complex
  ! Note: internal subroutine of psb_zilu0_fact
  !
  !  This routine copies a row of a sparse matrix A, stored in the psb_zspmat_type
  !  data structure a, into the arrays lval and uval and into the scalar variable
  !  dia, corresponding to the lower and upper triangles of A and to the diagonal
  !  entry of the row, respectively. The entries in lval and uval are stored
  !  according to the CSR format; the corresponding column indices are stored in
  !  the arrays lja and uja.
  !
  !  If the sparse matrix is in CSR format, a 'straight' copy is performed;
  !  otherwise psb_sp_getblk is used to extract a block of rows, which is then
  !  copied into lval, dia, uval row by row, through successive calls to
  !  ilu_copyin.
  !
  !  The routine is used by psb_zilu0_factint in the computation of the ILU(0)/MILU(0)
  !  factorization of a local sparse matrix.
  !
  !  TODO: modify the routine to allow copying into output L and U that are
  !  already filled with indices; this would allow computing an ILU(k) pattern,
  !  then use the ILU(0) internal for subsequent calls with the same pattern.
  !
  ! Arguments:
  !    i       -  integer, input.
  !               The local index of the row to be extracted from  the
  !               sparse matrix structure a.
  !    m       -  integer, input.
  !               The number of rows of the local matrix stored into a.
  !    a       -  type(psb_zspmat_type), input.
  !               The sparse matrix structure containing the row to be copied.
  !    jd      -  integer, input.
  !               The column index of the diagonal entry of the row to be
  !               copied.
  !    jmin    -  integer, input.
  !               Minimum valid column index.
  !    jmax    -  integer, input.
  !               Maximum valid column index.
  !               The output matrix will contain a clipped copy taken from
  !               a(1:m,jmin:jmax).
  !    l1      -  integer, input/output.
  !               Pointer to the last occupied entry of lval.
  !    lja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the lower triangle
  !               copied in lval row by row (see psb_zilu0_factint), according
  !               to the CSR storage format.
  !    lval   -  complex(psb_dpk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               lower triangle are copied.
  !    dia     -  complex(psb_dpk_), output.
  !               The diagonal entry of the copied row.
  !    l2      -  integer, input/output.
  !               Pointer to the last occupied entry of uval.
  !    uja    -  integer, dimension(:), input/output.
  !               The column indices of the nonzero entries of the upper triangle
  !               copied in uval row by row (see psb_zilu0_factint), according
  !               to the CSR storage format.
  !    uval   -  complex(psb_dpk_), dimension(:), input/output.
  !               The array where the entries of the row corresponding to the
  !               upper triangle are copied.
  !    ktrw    -  integer, input/output.
  !               The index identifying the last entry taken from the
  !               staging buffer trw. See below.
  !    trw     -  type(psb_zspmat_type), input/output.
  !               A staging buffer. If the matrix A is not in CSR format, we use
  !               the psb_sp_getblk routine and store its output in trw; when we
  !               need to call psb_sp_getblk we do it for a block of rows, and then
  !               we consume them from trw in successive calls to this routine,
  !               until we empty the buffer. Thus we will make a call to psb_sp_getblk
  !               every nrb calls to copyin. If A is in CSR format it is unused.
  !
  subroutine ilu_copyin(i,m,aa,jd,jmin,jmax,l1,lja,lval,&
       & dia,l2,uja,uval,ktrw,trw,upd,shft)

    use psb_base_mod

    implicit none

    ! Arguments
    class(psb_z_base_sparse_mat),intent(inout)  :: aa
    type(psb_z_coo_sparse_mat), intent(inout) :: trw
    integer(psb_ipk_), intent(in)        :: i,m,jd,jmin,jmax
    integer(psb_ipk_), intent(inout)     :: ktrw,l1,l2
    integer(psb_ipk_), intent(inout)     :: lja(:), uja(:)
    complex(psb_dpk_), intent(inout)       :: lval(:), uval(:), dia
    character, intent(in)                :: upd
    complex(psb_dpk_), intent(in)       :: shft
    ! Local variables
    integer(psb_ipk_)             :: k,j,info,irb, nz
    integer(psb_ipk_), parameter  :: nrb=40
    character(len=20), parameter  :: name='ilu_copyin'
    character(len=20)             :: ch_err

    info=psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if
    if (psb_toupper(upd) == 'F') then

      select type(aa)
      type is (psb_z_csr_sparse_mat)

        !
        ! Take a fast shortcut if the matrix is stored in CSR format
        !

        do j = aa%irp(i), aa%irp(i+1) - 1
          k = aa%ja(j)
          !           write(0,*)'KKKKK',k
          if ((k < jd).and.(k >= jmin)) then
            l1 = l1 + 1
            lval(l1) = aa%val(j)
            lja(l1) = k
          else if (k == jd) then
            dia = aa%val(j) + shft
          else if ((k > jd).and.(k <= jmax)) then
            l2 = l2 + 1
            uval(l2) = aa%val(j)
            uja(l2) = k
          end if
        enddo

      class default

        !
        ! Otherwise use psb_sp_getblk, slower but able (in principle) of
        ! handling any format. In this case, a block of rows is extracted
        ! instead of a single row, for performance reasons, and these
        ! rows are copied one by one into lval, dia, uval, through
        ! successive calls to ilu_copyin.
        !

        if ((mod(i,nrb) == 1).or.(nrb == 1)) then
          irb = min(m-i+1,nrb)
          call aa%csget(i,i+irb-1,trw,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='csget'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          ktrw=1
        end if

        nz = trw%get_nzeros()
        do
          if (ktrw > nz) exit
          if (trw%ia(ktrw) > i) exit
          k = trw%ja(ktrw)
          if ((k < jd).and.(k >= jmin)) then
            l1 = l1 + 1
            lval(l1) = trw%val(ktrw)
            lja(l1) = k
          else if (k == jd) then
            dia = trw%val(ktrw) + shft
          else if ((k > jd).and.(k <= jmax)) then
            l2 = l2 + 1
            uval(l2) = trw%val(ktrw)
            uja(l2) = k
          end if
          ktrw = ktrw + 1
        enddo

      end select

    else

      write(0,*) 'Update not implemented '
      info = 31
      call psb_errpush(info,name,&
           & i_err=(/ione*13,izero,izero,izero,izero/),a_err=upd)
      goto 9999

    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine ilu_copyin

end subroutine psb_zilu0_fact
