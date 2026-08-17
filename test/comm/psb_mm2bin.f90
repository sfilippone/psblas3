!
! Convert a Matrix Market file into a row-sorted binary CSR-like layout that can
! be read in parallel.
!
! Why this exists: Matrix Market is ASCII with variable-length lines, so the byte
! offset of nonzero k cannot be computed without scanning the file. Reading it in
! parallel is therefore impossible without a redistribution step. Converting once
! to a fixed-record binary layout lets every rank seek straight to its own rows,
! which also removes the need for psb_matdist -- at 448 ranks that scatter cost
! 28 s of MPI time on its own.
!
! Layout, all 8-byte native-endian:
!
!   header   magic, version, ipk_bytes, rpk_bytes, nrows, ncols, nnz, reserved
!   row_ptr  nrows+1 entries, 1-based, row i occupies row_ptr(i)..row_ptr(i+1)-1
!   col_idx  nnz entries, 1-based global column indices
!   values   nnz doubles
!
! Native endianness on purpose: this is an intermediate for local runs, and the
! header records the kinds so a mismatched reader fails loudly instead of
! producing garbage.
!
program psb_mm2bin
  use psb_base_mod
  use psb_util_mod
  implicit none

  integer(psb_lpk_), parameter :: bin_magic = 6510534141251461121_psb_lpk_
  integer(psb_lpk_), parameter :: bin_version = 1_psb_lpk_

  type(psb_ldspmat_type)    :: aux
  type(psb_ld_coo_sparse_mat) :: coo
  character(len=1024)       :: infile, outfile
  integer(psb_ipk_)         :: info, iout
  integer(psb_lpk_)         :: nr, nc, nz, k, irow
  integer(psb_lpk_), allocatable :: row_ptr(:)
  real(psb_dpk_)            :: t0, t1

  if (command_argument_count() < 2) then
    write(psb_err_unit,*) 'usage: psb_mm2bin <in.mtx> <out.psb>'
    stop 1
  end if
  call get_command_argument(1, infile)
  call get_command_argument(2, outfile)

  call cpu_time(t0)
  call mm_mat_read(aux, info, filename=trim(infile))
  if (info /= psb_success_) then
    write(psb_err_unit,*) 'cannot read ', trim(infile), ' info=', info
    stop 2
  end if
  call cpu_time(t1)
  write(psb_out_unit,'("read      : ",f8.2," s")') t1-t0

  nr = aux%get_nrows()
  nc = aux%get_ncols()

  ! Move the storage into a COO object and sort it row-major: fix() is what
  ! guarantees the entries of a row are contiguous, which is the whole premise
  ! of the layout below.
  call cpu_time(t0)
  call aux%mv_to(coo)
  call coo%fix(info)
  if (info /= psb_success_) then
    write(psb_err_unit,*) 'fix failed, info=', info
    stop 3
  end if
  nz = coo%get_nzeros()
  call cpu_time(t1)
  write(psb_out_unit,'("sort      : ",f8.2," s")') t1-t0
  write(psb_out_unit,'("rows ",i0,"  cols ",i0,"  nonzeros ",i0)') nr, nc, nz

  ! row_ptr from the sorted row indices
  allocate(row_ptr(nr+1), stat=info)
  if (info /= 0) stop 4
  row_ptr = 0
  do k = 1, nz
    irow = coo%ia(k)
    row_ptr(irow+1) = row_ptr(irow+1) + 1
  end do
  row_ptr(1) = 1
  do k = 1, nr
    row_ptr(k+1) = row_ptr(k+1) + row_ptr(k)
  end do
  if (row_ptr(nr+1) /= nz+1) then
    write(psb_err_unit,*) 'row_ptr inconsistent: ', row_ptr(nr+1), nz+1
    stop 5
  end if

  call cpu_time(t0)
  iout = 91
  open(unit=iout, file=trim(outfile), access='stream', form='unformatted', &
       & status='replace', action='write', iostat=info)
  if (info /= 0) stop 6
  write(iout) bin_magic, bin_version, &
       & int(storage_size(nz)/8, psb_lpk_), int(storage_size(coo%val(1))/8, psb_lpk_), &
       & nr, nc, nz, 0_psb_lpk_
  write(iout) row_ptr
  write(iout) coo%ja(1:nz)
  write(iout) coo%val(1:nz)
  close(iout)
  call cpu_time(t1)
  write(psb_out_unit,'("write     : ",f8.2," s")') t1-t0
  write(psb_out_unit,'("output    : ",a)') trim(outfile)

end program psb_mm2bin
