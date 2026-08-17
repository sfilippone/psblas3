!
! Parallel read of the binary layout written by psb_mm2bin.
!
! Every rank opens the file, reads the header and only its own slice of row_ptr,
! then seeks straight to the nonzeros of the rows it owns. Nothing is read twice
! and nothing is redistributed: the block-row partition is implicit in the file
! order, so psb_matdist -- 28 s of MPI time at 448 ranks -- is not needed at all.
!
! Positioned stream reads are used rather than MPI-IO. The ranks touch disjoint
! byte ranges, so there is nothing to coordinate; if collective buffering ever
! turns out to matter at high rank counts, swapping these for
! MPI_File_read_at_all is a change local to this file.
!
module psb_binmat_mod
  use psb_base_mod
  implicit none

  integer(psb_lpk_), parameter :: psb_bin_magic   = 6510534141251461121_psb_lpk_
  integer(psb_lpk_), parameter :: psb_bin_version = 1_psb_lpk_

contains

  subroutine psb_bin_mat_read(ctxt, fname, afmt, a, desc_a, info)
    type(psb_ctxt_type), intent(in)      :: ctxt
    character(len=*), intent(in)         :: fname, afmt
    type(psb_dspmat_type), intent(out)   :: a
    type(psb_desc_type), intent(out)     :: desc_a
    integer(psb_ipk_), intent(out)       :: info

    integer(psb_ipk_) :: iam, np, iu, nloc, ist
    integer(psb_lpk_) :: magic, vers, ipkb, rpkb, nr, nc, nz, res
    integer(psb_lpk_) :: r1, r2, k1, k2, lnz, i, k, base, cbase, vbase, nb, rem
    integer(psb_lpk_), allocatable :: rp(:), ia(:), ja(:)
    real(psb_dpk_), allocatable    :: val(:)

    info = psb_success_
    call psb_info(ctxt, iam, np)

    iu = 71 + iam
    open(unit=iu, file=trim(fname), access='stream', form='unformatted', &
         & status='old', action='read', iostat=ist)
    if (ist /= 0) then
      if (iam == psb_root_) write(psb_err_unit,*) 'cannot open ', trim(fname)
      info = psb_err_internal_error_
      return
    end if

    read(iu, pos=1, iostat=ist) magic, vers, ipkb, rpkb, nr, nc, nz, res
    if ((ist /= 0).or.(magic /= psb_bin_magic)) then
      if (iam == psb_root_) write(psb_err_unit,*) 'not a psb binary matrix: ', trim(fname)
      info = psb_err_internal_error_ ; close(iu) ; return
    end if
    ! The header records the kinds the file was written with: a mismatch has to
    ! fail here rather than silently reinterpret the bytes.
    if ((ipkb /= 8).or.(rpkb /= 8)) then
      if (iam == psb_root_) write(psb_err_unit,*) 'kind mismatch in ', trim(fname), ipkb, rpkb
      info = psb_err_internal_error_ ; close(iu) ; return
    end if
    if (nr /= nc) then
      if (iam == psb_root_) write(psb_err_unit,*) 'matrix is not square'
      info = psb_err_internal_error_ ; close(iu) ; return
    end if

    ! Block-row split, in rank order: the same assignment psb_cdall(nl=) makes,
    ! which is what lets the file order stand in for the partition.
    nb  = nr/np
    rem = mod(nr,np)
    if (iam < rem) then
      nloc = int(nb+1, psb_ipk_)
      r1   = iam*(nb+1) + 1
    else
      nloc = int(nb, psb_ipk_)
      r1   = rem*(nb+1) + (iam-rem)*nb + 1
    end if
    r2 = r1 + nloc - 1

    base  = 65                          ! row_ptr starts here (1-based stream pos)
    cbase = base + (nr+1)*8             ! col_idx
    vbase = cbase + nz*8                ! values

    allocate(rp(nloc+1), stat=ist)
    if (ist /= 0) then
      info = psb_err_alloc_dealloc_ ; close(iu) ; return
    end if
    read(iu, pos=base + (r1-1)*8, iostat=ist) rp
    if (ist /= 0) then
      info = psb_err_internal_error_ ; close(iu) ; return
    end if

    k1  = rp(1)
    k2  = rp(nloc+1) - 1
    lnz = k2 - k1 + 1
    if (lnz < 0) lnz = 0

    allocate(ia(max(lnz,1_psb_lpk_)), ja(max(lnz,1_psb_lpk_)), &
         &   val(max(lnz,1_psb_lpk_)), stat=ist)
    if (ist /= 0) then
      info = psb_err_alloc_dealloc_ ; close(iu) ; return
    end if

    if (lnz > 0) then
      read(iu, pos=cbase + (k1-1)*8, iostat=ist) ja(1:lnz)
      if (ist == 0) read(iu, pos=vbase + (k1-1)*8, iostat=ist) val(1:lnz)
      if (ist /= 0) then
        info = psb_err_internal_error_ ; close(iu) ; return
      end if
      ! Global row index for every nonzero, from the row_ptr slice.
      do i = 1, nloc
        do k = rp(i), rp(i+1)-1
          ia(k-k1+1) = r1 + i - 1
        end do
      end do
    end if
    close(iu)

    call psb_cdall(ctxt, desc_a, info, nl=nloc)
    if (info == psb_success_) call psb_spall(a, desc_a, info, nnz=int(lnz,psb_ipk_))
    if ((info == psb_success_).and.(lnz > 0)) &
         & call psb_spins(int(lnz,psb_ipk_), ia, ja, val, a, desc_a, info)
    if (info == psb_success_) call psb_cdasb(desc_a, info)
    if (info == psb_success_) call psb_spasb(a, desc_a, info, afmt=afmt)

  end subroutine psb_bin_mat_read

end module psb_binmat_mod
