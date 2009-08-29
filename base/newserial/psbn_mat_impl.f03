subroutine psbn_d_spcnv(a,b,info,type,mold,upd,dupl)
  use psbn_d_mat_mod, psb_protect_name => psbn_d_spcnv
  use psb_realloc_mod
  use psb_sort_mod
  type(psbn_d_sparse_mat), intent(in)    :: a
  type(psbn_d_sparse_mat), intent(out)   :: b
  integer, intent(out)                   :: info
  integer,optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold
  
end subroutine psbn_d_spcnv

subroutine psbn_d_spcnv_ip(a,info,type,mold,dupl)
  use psbn_d_mat_mod, psb_protect_name => psbn_d_spcnv_ip
  use psb_realloc_mod
  use psb_sort_mod
  
  type(psbn_d_sparse_mat), intent(inout)  :: a
  integer, intent(out)                    :: info
  integer,optional, intent(in)            :: dupl
  character(len=*), optional, intent(in)  :: type
  class(psbn_d_base_sparse_mat), intent(in), optional :: mold

end subroutine psbn_d_spcnv_ip
