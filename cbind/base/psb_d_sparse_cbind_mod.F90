! Distributed sparse operations for C callers. The numerical operations are
! delegated to PSBLAS; only handles, descriptors and storage are adapted here.
module psb_d_sparse_cbind_mod
  use iso_c_binding
  use psb_base_mod
  use psb_objhandle_mod
  use psb_cpenv_mod
  use psb_sparse_cbind_mod
  implicit none
  private
  public :: psb_c_dspall_2desc, psb_c_dspins_2desc, psb_c_dspasb_2desc
  public :: psb_c_dglob_transpose, psb_c_dpar_spspmm
  public :: psb_c_dspget_coo, psb_c_dsprowscl
  public :: psb_c_dspmm_2desc

contains

  function psb_c_dspall_2desc(ah, drh, dch) bind(c) result(res)
    type(psb_c_dspmat) :: ah
    type(psb_c_descriptor) :: drh, dch
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc
    type(psb_dspmat_type), pointer :: ap
    integer(psb_ipk_) :: info, stat

    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    res = info
    if (info /= psb_success_) return
    res = psb_c_sparse_check_pair(dr, dc, .false.)
    if (res /= psb_success_) return
    res = psb_err_invalid_cd_state_
    if (c_associated(drh%item, dch%item)) return
    if (.not.dc%is_bld()) return
    res = psb_err_invalid_mat_state_
    if (c_associated(ah%item)) return

    allocate(ap, stat=stat)
    res = psb_err_alloc_dealloc_
    if (stat /= 0) return
    call psb_spall(ap, dr, info)
    if (info == psb_success_) then
      call ap%set_ncols(dc%get_local_cols())
      ah%item = c_loc(ap)
    else
      call ap%free()
      deallocate(ap)
    end if
    res = info
  end function psb_c_dspall_2desc

  function psb_c_dspins_2desc(nz, irw, icl, val, ah, drh, dch) bind(c) result(res)
    integer(psb_c_ipk_), value :: nz
    integer(psb_c_lpk_), intent(in) :: irw(*), icl(*)
    real(c_double), intent(in) :: val(*)
    type(psb_c_dspmat) :: ah
    type(psb_c_descriptor) :: drh, dch
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc
    type(psb_dspmat_type), pointer :: ap
    integer(psb_ipk_) :: info, ixb, stat
    integer(psb_ipk_), allocatable :: local_rows(:)
    integer(psb_lpk_), allocatable :: rows(:), cols(:)

    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    res = info
    if (info /= psb_success_) return
    res = psb_c_sparse_check_pair(dr, dc, .false.)
    if (res /= psb_success_) return
    res = psb_err_invalid_cd_state_
    if (c_associated(drh%item, dch%item)) return
    if (.not.dc%is_bld()) return
    res = psb_err_invalid_mat_state_
    if (.not.c_associated(ah%item)) return
    call c_f_pointer(ah%item, ap)
    if (.not.ap%is_bld()) return
    res = psb_err_invalid_matrix_sizes_
    if (ap%get_nrows() /= dr%get_local_rows()) return
    res = psb_err_iarg_invalid_value_
    if (nz < 0) return
    if (nz == 0) then
      res = psb_success_
      return
    end if

    allocate(rows(nz), cols(nz), local_rows(nz), stat=stat)
    res = psb_err_alloc_dealloc_
    if (stat /= 0) return
    ixb = psb_c_get_index_base()
    rows = irw(1:nz) + int(1-ixb, psb_lpk_)
    cols = icl(1:nz) + int(1-ixb, psb_lpk_)
    res = psb_err_iarg_invalid_value_
    if (any(rows < 1).or.any(rows > dr%get_global_rows())) return
    if (any(cols < 1).or.any(cols > dc%get_global_rows())) return
    call dr%g2l(rows, local_rows, info, owned=.true.)
    res = info
    if (info /= psb_success_) return
    ! The native two-descriptor insertion silently discards remote rows.
    ! Reject them here: callers must first route entries to their row owners.
    res = psb_err_iarg_invalid_value_
    if (any(local_rows < 1)) return
    call psb_spins(nz, rows, cols, val(1:nz), ap, dr, dc, info)
    res = info
  end function psb_c_dspins_2desc

  function psb_c_dspasb_2desc(ah, drh, dch) bind(c) result(res)
    type(psb_c_dspmat) :: ah
    type(psb_c_descriptor) :: drh, dch
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc
    type(psb_dspmat_type), pointer :: ap
    integer(psb_ipk_) :: info

    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    res = info
    if (info /= psb_success_) return
    res = psb_c_sparse_check_pair(dr, dc, .true.)
    if (res /= psb_success_) return
    res = psb_err_invalid_mat_state_
    if (.not.c_associated(ah%item)) return
    call c_f_pointer(ah%item, ap)
    if (.not.ap%is_bld()) return
    res = psb_err_invalid_matrix_sizes_
    if (ap%get_nrows() /= dr%get_local_rows()) return
    call ap%set_ncols(dc%get_local_cols())
    call ap%cscnv(info, type='CSR', dupl=psb_dupl_add_)
    res = info
  end function psb_c_dspasb_2desc

  function psb_c_dglob_transpose(ah, drh, dch, bh, dbch, trans) bind(c) result(res)
    type(psb_c_dspmat) :: ah, bh
    type(psb_c_descriptor) :: drh, dch, dbch
    character(c_char), intent(in) :: trans(*)
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc, dbc
    type(psb_dspmat_type), pointer :: ap, bp
    type(psb_d_coo_sparse_mat) :: acoo, bcoo
    integer(psb_ipk_) :: info, stat, cleanup_info
    character :: trans_

    nullify(bp, dbc)
    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    if (info == psb_success_) info = psb_c_sparse_check_pair(dr, dc, .true.)
    if (info == psb_success_) then
      if (.not.c_associated(ah%item)) then
        info = psb_err_invalid_mat_state_
      else
        call c_f_pointer(ah%item, ap)
        if (.not.ap%is_asb()) info = psb_err_invalid_mat_state_
        if ((ap%get_nrows() /= dr%get_local_rows()).or. &
            (ap%get_ncols() /= dc%get_local_cols())) info = psb_err_invalid_matrix_sizes_
      end if
    end if
    if (c_associated(bh%item).or.c_associated(dbch%item)) info = psb_err_invalid_mat_state_
    trans_ = psb_toupper(trans(1))
    if ((trans_ /= 'T').and.(trans_ /= 'C')) info = psb_err_iarg_invalid_value_
    call psb_c_sparse_agree(dr, info)
    if (info /= psb_success_) goto 900

    allocate(bp, dbc, stat=stat)
    if (stat /= 0) info = psb_err_alloc_dealloc_
    if (info == psb_success_) call ap%cscnv(acoo, info)
    call psb_c_sparse_agree(dr, info)
    if (info /= psb_success_) goto 900

    ! The output-column descriptor contains the new transpose halo; neither
    ! input descriptor nor the source matrix is changed.
    if (c_associated(drh%item, dch%item)) then
      call psb_glob_transpose(acoo, dr, info, atrans=bcoo, desc_rx=dbc)
    else
      call psb_glob_transpose(acoo, dr, info, atrans=bcoo, desc_c=dc, desc_rx=dbc)
    end if
    call psb_c_sparse_agree(dr, info)
    if (info /= psb_success_) goto 900
    ! Transpose and conjugate transpose coincide for real entries.
    call bp%mv_from(bcoo)
    call bp%cscnv(info, type='CSR')
    call psb_c_sparse_agree(dr, info)
    if (info /= psb_success_) goto 900
    bh%item = c_loc(bp)
    dbch%item = c_loc(dbc)
    nullify(bp, dbc)

900 continue
    if (associated(bp)) then
      call bp%free()
      deallocate(bp)
    end if
    if (associated(dbc)) then
      call dbc%free(cleanup_info)
      deallocate(dbc)
    end if
    call acoo%free()
    call bcoo%free()
    res = info
  end function psb_c_dglob_transpose

  function psb_c_dpar_spspmm(ah, darh, dach, bh, dbrh, dbch, ch, dcch) bind(c) result(res)
    type(psb_c_dspmat) :: ah, bh, ch
    type(psb_c_descriptor) :: darh, dach, dbrh, dbch, dcch
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dar, dac, dbr, dbc, dcc
    type(psb_dspmat_type), pointer :: ap, bp, cp
    type(psb_d_csr_sparse_mat) :: acsr, bcsr, ccsr
    integer(psb_ipk_) :: info, stat, cleanup_info

    nullify(cp, dcc)
    call psb_c_sparse_get_desc(darh, dar, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dach, dac, info)
    if (info == psb_success_) call psb_c_sparse_get_desc(dbrh, dbr, info)
    if (info == psb_success_) call psb_c_sparse_get_desc(dbch, dbc, info)
    if (info == psb_success_) info = psb_c_sparse_check_pair(dar, dac, .true.)
    if (info == psb_success_) info = psb_c_sparse_check_pair(dac, dbr, .true.)
    if (info == psb_success_) info = psb_c_sparse_check_pair(dbr, dbc, .true.)
    if (info == psb_success_) then
      if (.not.psb_c_sparse_same_owned(dac, dbr, info)) then
        if (info == psb_success_) info = psb_err_invalid_matrix_sizes_
      end if
      if (.not.(c_associated(ah%item).and.c_associated(bh%item))) then
        info = psb_err_invalid_mat_state_
      else
        call c_f_pointer(ah%item, ap)
        call c_f_pointer(bh%item, bp)
        if (.not.(ap%is_asb().and.bp%is_asb())) info = psb_err_invalid_mat_state_
        if ((ap%get_nrows() /= dar%get_local_rows()).or. &
            (ap%get_ncols() /= dac%get_local_cols()).or. &
            (bp%get_nrows() /= dbr%get_local_rows()).or. &
            (bp%get_ncols() /= dbc%get_local_cols())) info = psb_err_invalid_matrix_sizes_
      end if
    end if
    if (c_associated(ch%item).or.c_associated(dcch%item)) info = psb_err_invalid_mat_state_
    call psb_c_sparse_agree(dar, info)
    if (info /= psb_success_) goto 900

    allocate(cp, dcc, stat=stat)
    if (stat /= 0) info = psb_err_alloc_dealloc_
    if (info == psb_success_) call ap%cscnv(acsr, info)
    if (info == psb_success_) call bp%cscnv(bcsr, info)
    if (info == psb_success_) call dbc%clone(dcc, info)
    call psb_c_sparse_agree(dar, info)
    if (info /= psb_success_) goto 900
    call psb_cd_reinit(dcc, info)
    call psb_c_sparse_agree(dar, info)
    if (info /= psb_success_) goto 900

    ! Native SpGEMM extends B with halo rows and expands its column map. Both
    ! are private copies, so callers retain their original operators/layouts.
    call psb_par_spspmm(acsr, dac, bcsr, ccsr, dcc, info)
    call psb_c_sparse_agree(dar, info)
    if (info /= psb_success_) goto 900
    call psb_cdasb(dcc, info)
    call psb_c_sparse_agree(dar, info)
    if (info /= psb_success_) goto 900
    call ccsr%set_nrows(dar%get_local_rows())
    call ccsr%set_ncols(dcc%get_local_cols())
    call cp%mv_from(ccsr)
    call cp%set_asb()
    ch%item = c_loc(cp)
    dcch%item = c_loc(dcc)
    nullify(cp, dcc)

900 continue
    if (associated(cp)) then
      call cp%free()
      deallocate(cp)
    end if
    if (associated(dcc)) then
      call dcc%free(cleanup_info)
      deallocate(dcc)
    end if
    call acsr%free()
    call bcsr%free()
    call ccsr%free()
    res = info
  end function psb_c_dpar_spspmm

  function psb_c_dspmm_2desc(alpha, ah, xh, beta, yh, drh, dch, trans) bind(c) result(res)
    real(c_double), value :: alpha, beta
    type(psb_c_dspmat) :: ah
    type(psb_c_dvector) :: xh, yh
    type(psb_c_descriptor) :: drh, dch
    character(c_char), intent(in) :: trans(*)
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc
    type(psb_dspmat_type), pointer :: ap
    type(psb_d_vect_type), pointer :: xp, yp
    integer(psb_ipk_) :: info, nr, nc, nowned
    character :: trans_

    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    if (info == psb_success_) info = psb_c_sparse_check_pair(dr, dc, .true.)
    trans_ = psb_toupper(trans(1))
    if ((trans_ /= 'N').and.(trans_ /= 'T').and.(trans_ /= 'C')) info = psb_err_iarg_invalid_value_
    if (info == psb_success_) then
      if (.not.(c_associated(ah%item).and.c_associated(xh%item).and.c_associated(yh%item))) then
        info = psb_err_invalid_mat_state_
      else if (c_associated(xh%item, yh%item)) then
        info = psb_err_iarg_invalid_value_
      else
        call c_f_pointer(ah%item, ap)
        call c_f_pointer(xh%item, xp)
        call c_f_pointer(yh%item, yp)
        if (.not.ap%is_asb()) info = psb_err_invalid_mat_state_
        if (.not.(xp%is_asb().and.yp%is_asb())) info = psb_err_invalid_vect_state_
        nr = dr%get_local_rows()
        nc = dc%get_local_cols()
        nowned = dc%get_local_rows()
        if ((ap%get_nrows() /= nr).or.(ap%get_ncols() /= nc)) info = psb_err_invalid_matrix_sizes_
        if (trans_ == 'N') then
          if ((xp%get_nrows() < nowned).or.(yp%get_nrows() < nr)) info = psb_err_invalid_matrix_sizes_
          if ((info == psb_success_).and.(xp%get_nrows() < nc)) call xp%reall(nc, info)
        else
          if ((xp%get_nrows() < nr).or.(yp%get_nrows() < nowned)) info = psb_err_invalid_matrix_sizes_
          if ((info == psb_success_).and.(yp%get_nrows() < nc)) call yp%reall(nc, info)
        end if
      end if
    end if
    call psb_c_sparse_agree(dr, info)
    if (info /= psb_success_) then
      res = info
      return
    end if

    if (trans_ == 'N') then
      call psb_halo(xp, dc, info)
      if (info == psb_success_) call psb_csmm(alpha, ap, xp, beta, yp, info)
    else
      ! Only owned y entries carry beta*y. Old halo values must not be
      ! multiplied or added back to their owners on a repeated application.
      if (nc > nowned) call yp%set(0.0_psb_dpk_, nowned+1, nc)
      call psb_csmm(alpha, ap, xp, beta, yp, info, trans=trans_)
      call psb_c_sparse_agree(dr, info)
      if (info == psb_success_) call psb_halo(yp, dc, info, tran='T')
    end if
    res = info
  end function psb_c_dspmm_2desc

  function psb_c_dspget_coo(ah, drh, dch, capacity, nz, irw, icl, val) bind(c) result(res)
    type(psb_c_dspmat) :: ah
    type(psb_c_descriptor) :: drh, dch
    integer(psb_c_ipk_), value :: capacity
    integer(psb_c_ipk_), intent(out) :: nz
    integer(psb_c_lpk_), intent(out) :: irw(*), icl(*)
    real(c_double), intent(out) :: val(*)
    integer(psb_c_ipk_) :: res
    type(psb_desc_type), pointer :: dr, dc
    type(psb_dspmat_type), pointer :: ap
    type(psb_d_coo_sparse_mat) :: acoo
    integer(psb_lpk_), allocatable :: rows(:), cols(:)
    integer(psb_ipk_) :: info, nnz, stat, ixb

    nz = 0
    call psb_c_sparse_get_desc(drh, dr, info)
    res = info
    if (info /= psb_success_) return
    call psb_c_sparse_get_desc(dch, dc, info)
    res = info
    if (info /= psb_success_) return
    res = psb_c_sparse_check_pair(dr, dc, .true.)
    if (res /= psb_success_) return
    res = psb_err_iarg_invalid_value_
    if (capacity < 0) return
    res = psb_err_invalid_mat_state_
    if (.not.c_associated(ah%item)) return
    call c_f_pointer(ah%item, ap)
    if (.not.ap%is_asb()) return
    res = psb_err_invalid_matrix_sizes_
    if ((ap%get_nrows() /= dr%get_local_rows()).or. &
        (ap%get_ncols() /= dc%get_local_cols())) return
    nnz = ap%get_nzeros()
    nz = nnz
    if ((capacity == 0).or.(nnz == 0)) then
      res = psb_success_
      return
    end if
    res = psb_err_iarg_invalid_value_
    if (capacity < nnz) return

    call ap%cscnv(acoo, info)
    res = info
    if (info /= psb_success_) return
    nnz = acoo%get_nzeros()
    nz = nnz
    allocate(rows(nnz), cols(nnz), stat=stat)
    res = psb_err_alloc_dealloc_
    if (stat /= 0) return
    rows = int(acoo%ia(1:nnz), psb_lpk_)
    cols = int(acoo%ja(1:nnz), psb_lpk_)
    call dr%l2gip(rows, info, owned=.true.)
    if (info == psb_success_) call dc%l2gip(cols, info, owned=.false.)
    res = info
    if (info /= psb_success_) return
    ixb = psb_c_get_index_base()
    irw(1:nnz) = rows - int(1-ixb, psb_lpk_)
    icl(1:nnz) = cols - int(1-ixb, psb_lpk_)
    val(1:nnz) = acoo%val(1:nnz)
    call acoo%free()
  end function psb_c_dspget_coo

  function psb_c_dsprowscl(ah, nr, val) bind(c) result(res)
    type(psb_c_dspmat) :: ah
    integer(psb_c_ipk_), value :: nr
    real(c_double), intent(in) :: val(*)
    integer(psb_c_ipk_) :: res
    type(psb_dspmat_type), pointer :: ap
    integer(psb_ipk_) :: info

    res = psb_err_invalid_mat_state_
    if (.not.c_associated(ah%item)) return
    call c_f_pointer(ah%item, ap)
    if (.not.ap%is_asb()) return
    res = psb_err_invalid_matrix_sizes_
    if ((nr < 0).or.(nr /= ap%get_nrows())) return
    if (nr == 0) then
      res = psb_success_
      return
    end if
    call ap%scal(val(1:nr), info, side='L')
    res = info
  end function psb_c_dsprowscl

end module psb_d_sparse_cbind_mod
