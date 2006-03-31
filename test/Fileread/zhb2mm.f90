program zhb2mm
  use psb_sparse_mod
  use mmio
  use hbio
  type(psb_zspmat_type) :: a
  
  integer n, nnz,info,i,j,k
  INTEGER :: iwflag,IOUT,NCOL,NELTVL,NNZERO,NRHS,NRHSIX,NROW,&
       &  iter
  CHARACTER  :: RHSDATATYPE,DATATYPE*3,KEY*8,OUTFILE*20,MTITLE*72

  
  call hb_read(a,info,mtitle=mtitle)

  call mm_mat_write(a,mtitle,info)

  stop


end program zhb2mm
  
