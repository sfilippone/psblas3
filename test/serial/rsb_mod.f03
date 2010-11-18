module rsb_mod
  use iso_c_binding

! module constants:

interface
integer(c_int) function &
  &rsb_init&
  &()&
  &bind(c,name='rsb_init')
use iso_c_binding
 end function rsb_init
end interface

interface
integer(c_int) function &
  &rsb_was_initialized&
  &()&
  &bind(c,name='rsb_was_initialized')
use iso_c_binding
 end function rsb_was_initialized
end interface

interface
integer(c_int) function &
  &rsb_exit&
  &()&
  &bind(c,name='rsb_exit')
use iso_c_binding
 end function rsb_exit
end interface

interface
type(c_ptr) function &
  &rsb_allocate_rsb_sparse_matrix_const&
  &(VAc,IAc,JAc,nnz,typecode,m,k,Mb,Kb,flags,errvalp)&
  &bind(c,name='rsb_allocate_rsb_sparse_matrix_const')
use iso_c_binding
 real(c_double) :: VAc(*)
 integer(c_int) :: IAc(*)
 integer(c_int) :: JAc(*)
 integer(c_int), value  :: nnz
 integer(c_int), value  :: typecode
 integer(c_int), value  :: m
 integer(c_int), value  :: k
 integer(c_int), value  :: Mb
 integer(c_int), value  :: Kb
 integer(c_int), value  :: flags
 integer(c_int) :: errvalp
 end function rsb_allocate_rsb_sparse_matrix_const
end interface

interface
type(c_ptr) function &
  &rsb_free_sparse_matrix&
  &(matrix)&
  &bind(c,name='rsb_free_sparse_matrix')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_free_sparse_matrix
end interface

interface
type(c_ptr) function &
  &rsb_clone&
  &(matrix)&
  &bind(c,name='rsb_clone')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_clone
end interface

interface
integer(c_int) function &
  &rsb_mark_matrix_with_type_flags&
  &(matrix)&
  &bind(c,name='rsb_mark_matrix_with_type_flags')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_mark_matrix_with_type_flags
end interface

interface
integer(c_int) function &
  &rsb_meminfo&
  &()&
  &bind(c,name='rsb_meminfo')
use iso_c_binding
 end function rsb_meminfo
end interface

interface
integer(c_int) function &
  &rsb_check_leak&
  &()&
  &bind(c,name='rsb_check_leak')
use iso_c_binding
 end function rsb_check_leak
end interface

interface
integer(c_int) function &
  &rsb_spmv&
  &(matrix,x,y,alphap,betap,incx,incy,transa)&
  &bind(c,name='rsb_spmv')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: incx
 integer(c_int), value  :: incy
 integer(c_int), value  :: transa
 end function rsb_spmv
end interface

interface
integer(c_int) function &
  &rsb_spmv_aa&
  &(matrix,x,y,transa)&
  &bind(c,name='rsb_spmv_aa')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 integer(c_int), value  :: transa
 end function rsb_spmv_aa
end interface

interface
integer(c_int) function &
  &rsb_spmv_sa&
  &(matrix,x,y,alphap,transa)&
  &bind(c,name='rsb_spmv_sa')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 real(c_double) :: alphap
 integer(c_int), value  :: transa
 end function rsb_spmv_sa
end interface

interface
integer(c_int) function &
  &rsb_spmv_unua&
  &(matrix,x,y,transa)&
  &bind(c,name='rsb_spmv_unua')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 integer(c_int), value  :: transa
 end function rsb_spmv_unua
end interface

interface
integer(c_int) function &
  &rsb_spmv_az&
  &(matrix,x,y,transa)&
  &bind(c,name='rsb_spmv_az')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 integer(c_int), value  :: transa
 end function rsb_spmv_az
end interface

interface
integer(c_int) function &
  &rsb_spmv_uxux&
  &(matrix,x,y,alphap,betap,transa)&
  &bind(c,name='rsb_spmv_uxux')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: transa
 end function rsb_spmv_uxux
end interface

interface
integer(c_int) function &
  &rsb_spmv_sxsx&
  &(matrix,x,y,alphap,betap,transa,incx,incy)&
  &bind(c,name='rsb_spmv_sxsx')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: transa
 integer(c_int), value  :: incx
 integer(c_int), value  :: incy
 end function rsb_spmv_sxsx
end interface

interface
integer(c_int) function &
  &rsb_infinity_norm&
  &(matrix,infinity_norm,transa)&
  &bind(c,name='rsb_infinity_norm')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: infinity_norm
 integer(c_int), value  :: transa
 end function rsb_infinity_norm
end interface

interface
integer(c_int) function &
  &rsb_one_norm&
  &(matrix,infinity_norm,transa)&
  &bind(c,name='rsb_one_norm')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: infinity_norm
 integer(c_int), value  :: transa
 end function rsb_one_norm
end interface

interface
integer(c_int) function &
  &rsb_rows_sums&
  &(matrix,d)&
  &bind(c,name='rsb_rows_sums')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 end function rsb_rows_sums
end interface

interface
integer(c_int) function &
  &rsb_columns_sums&
  &(matrix,d)&
  &bind(c,name='rsb_columns_sums')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 end function rsb_columns_sums
end interface

interface
integer(c_int) function &
  &rsb_absolute_rows_sums&
  &(matrix,d)&
  &bind(c,name='rsb_absolute_rows_sums')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 end function rsb_absolute_rows_sums
end interface

interface
integer(c_int) function &
  &rsb_absolute_columns_sums&
  &(matrix,d)&
  &bind(c,name='rsb_absolute_columns_sums')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 end function rsb_absolute_columns_sums
end interface

interface
integer(c_int) function &
  &rsb_spsv&
  &(matrix,x,y,alphap,incx,incy,transl)&
  &bind(c,name='rsb_spsv')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: x(*)
 real(c_double) :: y(*)
 real(c_double) :: alphap
 integer(c_int), value  :: incx
 integer(c_int), value  :: incy
 integer(c_int), value  :: transl
 end function rsb_spsv
end interface

interface
integer(c_int) function &
  &rsb_spmm_az&
  &(matrix,mrhs,mout,bstride,cstride,nrhs,transa)&
  &bind(c,name='rsb_spmm_az')
use iso_c_binding
 type(c_ptr), value  :: matrix
 type(c_ptr), value  :: mrhs
 type(c_ptr), value  :: mout
 integer(c_int), value  :: bstride
 integer(c_int), value  :: cstride
 integer(c_int), value  :: nrhs
 integer(c_int), value  :: transa
 end function rsb_spmm_az
end interface

interface
integer(c_int) function &
  &rsb_spmm_sxsx&
  &(matrix,b,c,ldb,ldc,nrhs,transa,alphap,betap,order)&
  &bind(c,name='rsb_spmm_sxsx')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: b(*)
 real(c_double) :: c(*)
 integer(c_int), value  :: ldb
 integer(c_int), value  :: ldc
 integer(c_int), value  :: nrhs
 integer(c_int), value  :: transa
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: order
 end function rsb_spmm_sxsx
end interface

interface
integer(c_int) function &
  &rsb_spmm&
  &(matrix,b,c,ldb,ldc,nrhs,transa,alphap,betap,order)&
  &bind(c,name='rsb_spmm')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: b(*)
 real(c_double) :: c(*)
 integer(c_int), value  :: ldb
 integer(c_int), value  :: ldc
 integer(c_int), value  :: nrhs
 integer(c_int), value  :: transa
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: order
 end function rsb_spmm
end interface

interface
integer(c_int) function &
  &rsb_spsm&
  &(matrix,b,ldb,nrhs,transt,alphap,betap,order)&
  &bind(c,name='rsb_spsm')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: b(*)
 integer(c_int), value  :: ldb
 integer(c_int), value  :: nrhs
 integer(c_int), value  :: transt
 real(c_double) :: alphap
 real(c_double) :: betap
 integer(c_int), value  :: order
 end function rsb_spsm
end interface

interface
type(c_ptr) function &
  &rsb_matrix_sum&
  &(matrixa,alphap,transa,matrixb,betap,transb,errvalp)&
  &bind(c,name='rsb_matrix_sum')
use iso_c_binding
 type(c_ptr), value  :: matrixa
 real(c_double) :: alphap
 integer(c_int), value  :: transa
 type(c_ptr), value  :: matrixb
 real(c_double) :: betap
 integer(c_int), value  :: transb
 integer(c_int) :: errvalp
 end function rsb_matrix_sum
end interface

interface
type(c_ptr) function &
  &rsb_matrix_mul&
  &(matrixa,alphap,transa,matrixb,betap,transb,errvalp)&
  &bind(c,name='rsb_matrix_mul')
use iso_c_binding
 type(c_ptr), value  :: matrixa
 real(c_double) :: alphap
 integer(c_int), value  :: transa
 type(c_ptr), value  :: matrixb
 real(c_double) :: betap
 integer(c_int), value  :: transb
 integer(c_int) :: errvalp
 end function rsb_matrix_mul
end interface

interface
integer(c_int) function &
  &rsb_matrix_add_to_dense&
  &(matrixa,alphap,transa,matrixb,ldb,nr,nc,rowmajor)&
  &bind(c,name='rsb_matrix_add_to_dense')
use iso_c_binding
 type(c_ptr), value  :: matrixa
 real(c_double) :: alphap
 integer(c_int), value  :: transa
 type(c_ptr), value  :: matrixb
 integer(c_int), value  :: ldb
 integer(c_int), value  :: nr
 integer(c_int), value  :: nc
 integer(c_int), value  :: rowmajor
 end function rsb_matrix_add_to_dense
end interface

interface
integer(c_int) function &
  &rsb_negation&
  &(matrix)&
  &bind(c,name='rsb_negation')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_negation
end interface

interface
integer(c_int) function &
  &rsb_scal&
  &(matrix,d,transa)&
  &bind(c,name='rsb_scal')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 integer(c_int), value  :: transa
 end function rsb_scal
end interface

interface
integer(c_int) function &
  &rsb_scale_rows&
  &(matrix,d)&
  &bind(c,name='rsb_scale_rows')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: d(*)
 end function rsb_scale_rows
end interface

interface
integer(c_int) function &
  &rsb_cest&
  &(IA,JA,nnz,typecode,m,k,p_r,p_c,M_b,K_b,flags)&
  &bind(c,name='rsb_cest')
use iso_c_binding
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int), value  :: nnz
 integer(c_int), value  :: typecode
 integer(c_int), value  :: m
 integer(c_int), value  :: k
 type(c_ptr), value  :: p_r
 type(c_ptr), value  :: p_c
 integer(c_int), value  :: M_b
 integer(c_int), value  :: K_b
 integer(c_int), value  :: flags
 end function rsb_cest
end interface

interface
integer(c_int) function &
  &rsb_perror&
  &(errval)&
  &bind(c,name='rsb_perror')
use iso_c_binding
 integer(c_int), value  :: errval
 end function rsb_perror
end interface

interface
integer(c_int) function &
  &rsb_sizeof&
  &(matrix)&
  &bind(c,name='rsb_sizeof')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_sizeof
end interface

interface
integer(c_int) function &
  &rsb_get_coo&
  &(matrix,VA,IA,JA,flags)&
  &bind(c,name='rsb_get_coo')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: VA(*)
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int), value  :: flags
 end function rsb_get_coo
end interface

interface
integer(c_int) function &
  &rsb_reinit&
  &(matrix)&
  &bind(c,name='rsb_reinit')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_reinit
end interface

interface
integer(c_int) function &
  &rsb_getdiag&
  &(matrix,diagonal)&
  &bind(c,name='rsb_getdiag')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: diagonal(*)
 end function rsb_getdiag
end interface

interface
integer(c_int) function &
  &rsb_get_sub_diag&
  &(matrix,diagonal,loffset)&
  &bind(c,name='rsb_get_sub_diag')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: diagonal(*)
 integer(c_int), value  :: loffset
 end function rsb_get_sub_diag
end interface

interface
integer(c_int) function &
  &rsb_get_supra_diag&
  &(matrix,diagonal,uoffset)&
  &bind(c,name='rsb_get_supra_diag')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: diagonal(*)
 integer(c_int), value  :: uoffset
 end function rsb_get_supra_diag
end interface

interface
integer(c_int) function &
  &rsb_get_row_dense&
  &(matrix,row,i)&
  &bind(c,name='rsb_get_row_dense')
use iso_c_binding
 type(c_ptr), value  :: matrix
 type(c_ptr), value  :: row
 integer(c_int), value  :: i
 end function rsb_get_row_dense
end interface

interface
integer(c_int) function &
  &rsb_get_rows_nnz&
  &(matrix,fr,lr,flags,errvalp)&
  &bind(c,name='rsb_get_rows_nnz')
use iso_c_binding
 type(c_ptr), value  :: matrix
 integer(c_int), value  :: fr
 integer(c_int), value  :: lr
 integer(c_int), value  :: flags
 integer(c_int) :: errvalp
 end function rsb_get_rows_nnz
end interface

interface
integer(c_int) function &
  &rsb_get_block_nnz&
  &(matrix,fr,lr,fc,lc,flags,errvalp)&
  &bind(c,name='rsb_get_block_nnz')
use iso_c_binding
 type(c_ptr), value  :: matrix
 integer(c_int), value  :: fr
 integer(c_int), value  :: lr
 integer(c_int), value  :: fc
 integer(c_int), value  :: lc
 integer(c_int), value  :: flags
 integer(c_int) :: errvalp
 end function rsb_get_block_nnz
end interface

interface
integer(c_int) function &
  &rsb_get_rows_sparse&
  &(matrix,VA,fr,lr,IA,JA,rnz,flags)&
  &bind(c,name='rsb_get_rows_sparse')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: VA(*)
 integer(c_int), value  :: fr
 integer(c_int), value  :: lr
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int) :: rnz
 integer(c_int), value  :: flags
 end function rsb_get_rows_sparse
end interface

interface
integer(c_int) function &
  &rsb_get_block_sparse_pattern&
  &(matrix,fr,lr,fc,lc,IA,JA,IREN,JREN,rnz,flags)&
  &bind(c,name='rsb_get_block_sparse_pattern')
use iso_c_binding
 type(c_ptr), value  :: matrix
 integer(c_int), value  :: fr
 integer(c_int), value  :: lr
 integer(c_int), value  :: fc
 integer(c_int), value  :: lc
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 type(c_ptr), value  :: IREN
 type(c_ptr), value  :: JREN
 integer(c_int) :: rnz
 integer(c_int), value  :: flags
 end function rsb_get_block_sparse_pattern
end interface

interface
integer(c_int) function &
  &rsb_get_block_sparse&
  &(matrix,VA,fr,lr,fc,lc,IA,JA,IREN,JREN,rnz,flags)&
  &bind(c,name='rsb_get_block_sparse')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: VA(*)
 integer(c_int), value  :: fr
 integer(c_int), value  :: lr
 integer(c_int), value  :: fc
 integer(c_int), value  :: lc
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 type(c_ptr), value  :: IREN
 type(c_ptr), value  :: JREN
 integer(c_int) :: rnz
 integer(c_int), value  :: flags
 end function rsb_get_block_sparse
end interface

interface
integer(c_int) function &
  &rsb_get_columns_sparse&
  &(matrix,columns,fc,lc,IA,JA,rnz,flags)&
  &bind(c,name='rsb_get_columns_sparse')
use iso_c_binding
 type(c_ptr), value  :: matrix
 type(c_ptr), value  :: columns
 integer(c_int), value  :: fc
 integer(c_int), value  :: lc
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int) :: rnz
 integer(c_int), value  :: flags
 end function rsb_get_columns_sparse
end interface

interface
integer(c_int) function &
  &rsb_assign&
  &(new_matrix,matrix)&
  &bind(c,name='rsb_assign')
use iso_c_binding
 type(c_ptr), value  :: new_matrix
 type(c_ptr), value  :: matrix
 end function rsb_assign
end interface

interface
integer(c_int) function &
  &rsb_sym_transpose&
  &(matrix)&
  &bind(c,name='rsb_sym_transpose')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_sym_transpose
end interface

interface
integer(c_int) function &
  &rsb_transpose&
  &(matrix)&
  &bind(c,name='rsb_transpose')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_transpose
end interface

interface
integer(c_int) function &
  &rsb_htranspose&
  &(matrix)&
  &bind(c,name='rsb_htranspose')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_htranspose
end interface

interface
integer(c_int) function &
  &rsb_get_matrix_nnz&
  &(matrix)&
  &bind(c,name='rsb_get_matrix_nnz')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_get_matrix_nnz
end interface

interface
integer(c_int) function &
  &rsb_elemental_scale&
  &(matrix,alphap)&
  &bind(c,name='rsb_elemental_scale')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: alphap
 end function rsb_elemental_scale
end interface

interface
integer(c_int) function &
  &rsb_elemental_scale_inv&
  &(matrix,alphap)&
  &bind(c,name='rsb_elemental_scale_inv')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: alphap
 end function rsb_elemental_scale_inv
end interface

interface
integer(c_int) function &
  &rsb_elemental_pow&
  &(matrix,alphap)&
  &bind(c,name='rsb_elemental_pow')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: alphap
 end function rsb_elemental_pow
end interface

interface
integer(c_int) function &
  &rsb_set_elements&
  &(matrix,VA,IA,JA,nnz)&
  &bind(c,name='rsb_set_elements')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: VA(*)
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int), value  :: nnz
 end function rsb_set_elements
end interface

interface
integer(c_int) function &
  &rsb_update_elements&
  &(matrix,VA,IA,JA,nnz,flags)&
  &bind(c,name='rsb_update_elements')
use iso_c_binding
 type(c_ptr), value  :: matrix
 real(c_double) :: VA(*)
 integer(c_int) :: IA(*)
 integer(c_int) :: JA(*)
 integer(c_int), value  :: nnz
 integer(c_int), value  :: flags
 end function rsb_update_elements
end interface

interface
integer(c_int) function &
  &rsb_psblas_trans_to_rsb_trans&
  &(trans)&
  &bind(c,name='rsb_psblas_trans_to_rsb_trans')
use iso_c_binding
 character(c_char), value  :: trans
 end function rsb_psblas_trans_to_rsb_trans
end interface

interface
integer(c_int) function &
  &rsb_print_matrix_t&
  &(matrix)&
  &bind(c,name='rsb_print_matrix_t')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_print_matrix_t
end interface

interface
integer(c_int) function &
  &rsb_print_matrix_unsorted_coo&
  &(matrix)&
  &bind(c,name='rsb_print_matrix_unsorted_coo')
use iso_c_binding
 type(c_ptr), value  :: matrix
 end function rsb_print_matrix_unsorted_coo
end interface

interface
type(c_ptr) function &
  &rsb_load_matrix_file_as_binary&
  &(filename,errvalp)&
  &bind(c,name='rsb_load_matrix_file_as_binary')
use iso_c_binding
 type(c_ptr), value  :: filename
 integer(c_int) :: errvalp
 end function rsb_load_matrix_file_as_binary
end interface

interface
integer(c_int) function &
  &rsb_save_matrix_file_as_binary&
  &(matrix,filename)&
  &bind(c,name='rsb_save_matrix_file_as_binary')
use iso_c_binding
 type(c_ptr), value  :: matrix
 type(c_ptr), value  :: filename
 end function rsb_save_matrix_file_as_binary
end interface

interface
integer(c_int) function &
  &rsb_save_matrix_file_as_matrix_market&
  &(matrix,filename)&
  &bind(c,name='rsb_save_matrix_file_as_matrix_market')
use iso_c_binding
 type(c_ptr), value  :: matrix
 type(c_ptr), value  :: filename
 end function rsb_save_matrix_file_as_matrix_market
end interface

interface
type(c_ptr) function &
  &rsb_load_matrix_file_as_matrix_market&
  &(filename,flags,typecode,errvalp)&
  &bind(c,name='rsb_load_matrix_file_as_matrix_market')
use iso_c_binding
 type(c_ptr), value  :: filename
 integer(c_int), value  :: flags
 integer(c_int), value  :: typecode
 integer(c_int) :: errvalp
 end function rsb_load_matrix_file_as_matrix_market
end interface
end module rsb_mod
