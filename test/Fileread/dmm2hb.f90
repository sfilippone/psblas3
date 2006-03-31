program dmm2hb
  use psb_sparse_mod
  use mmio
  use hbio
  type(psb_dspmat_type) :: a
  
  integer n, nnz,info,i,j,k


  nrhs        = 0
  nrhsix      = 0

  call mm_mat_read(a,info)

  call hb_write(a,info)

  stop
end program dmm2hb
  
