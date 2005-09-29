! Module containing interfaces for subroutine in SRC/F90/INTERNALS

module psi_mod

  use psb_descriptor_type

  interface
     subroutine psi_compute_size(desc_data,&
          & index_in, dl_lda, info)
       integer  :: info, dl_lda
       integer  :: desc_data(:), index_in(:)
     end subroutine psi_compute_size
  end interface

  interface
     subroutine psi_crea_bnd_elem(desc_a,info)
       use psb_descriptor_type
       type(psb_desc_type)  :: desc_a
       integer, intent(out) :: info
     end subroutine psi_crea_bnd_elem
  end interface

  interface
     subroutine psi_crea_index(desc_a,index_in,index_out,glob_idx,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)  :: desc_a
       integer, intent(out)             :: info
       integer, intent(in)              :: index_in(:)
       integer, pointer                 :: index_out(:)
       logical                          :: glob_idx
     end subroutine psi_crea_index
  end interface

  interface
     subroutine psi_crea_ovr_elem(desc_overlap,ovr_elem)
       integer :: desc_overlap(:)
       integer, pointer :: ovr_elem(:)
     end subroutine psi_crea_ovr_elem
  end interface
  
  interface
     subroutine psi_desc_index(desc_data,index_in,dep_list,&
          & length_dl,loc_to_glob,glob_to_loc,desc_index,&
          & isglob_in,info)
       integer :: desc_data(:),index_in(:),dep_list(:)
       integer :: loc_to_glob(:),glob_to_loc(:)
       integer,pointer :: desc_index(:)
       integer :: length_dl, info
       logical :: isglob_in
     end subroutine psi_desc_index
  end interface
  
  interface
     subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
       integer :: np,dep_list(:,:), l_dep_list(:), info
     end subroutine psi_sort_dl
  end interface

  interface psi_swapdata
     subroutine psi_dswapdatam(flag,n,beta,y,desc_a,work,info,data)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:,:), beta, work(:)
       type(psb_desc_type)  :: desc_a
       integer, optional    :: data
     end subroutine psi_dswapdatam
     subroutine psi_dswapdatav(flag,beta,y,desc_a,work,info,data)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:), beta, work(:)
       type(psb_desc_type)  :: desc_a
       integer, optional    :: data
     end subroutine psi_dswapdatav
     subroutine psi_iswapdatam(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       integer              :: y(:,:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_iswapdatam
     subroutine psi_iswapdatav(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       integer              :: y(:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_iswapdatav
  end interface


  interface psi_swaptran
     subroutine psi_dswaptranm(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:,:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_dswaptranm
     subroutine psi_dswaptranv(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_dswaptranv
     subroutine psi_iswaptranm(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       integer              :: y(:,:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_iswaptranm
     subroutine psi_iswaptranv(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       integer              :: y(:), beta, work(:)
       type(psb_desc_type)  :: desc_a
     end subroutine psi_iswaptranv
  end interface


  interface psi_gth
     subroutine psi_dgthm(n,k,idx,x,y)
       integer :: n, k, idx(:)
       real(kind(1.d0)) :: x(:,:), y(:)
     end subroutine psi_dgthm
     subroutine psi_dgthv(n,idx,x,y)
       integer :: n, idx(:)
       real(kind(1.d0)) :: x(:), y(:)
     end subroutine psi_dgthv
     subroutine psi_igthm(n,k,idx,x,y)
       integer :: n, k, idx(:)
       integer :: x(:,:), y(:)
     end subroutine psi_igthm
     subroutine psi_igthv(n,idx,x,y)
       integer :: n, idx(:)
       integer :: x(:), y(:)
     end subroutine psi_igthv
  end interface

  interface psi_sct
     subroutine psi_dsctm(n,k,idx,x,beta,y)
       integer :: n, k, idx(:)
       real(kind(1.d0)) :: beta, x(:), y(:,:)
     end subroutine psi_dsctm
     subroutine psi_dsctv(n,idx,x,beta,y)
       integer :: n, idx(:)
       real(kind(1.d0)) :: beta, x(:), y(:)
     end subroutine psi_dsctv
     subroutine psi_isctm(n,k,idx,x,beta,y)
       integer :: n, k, idx(:)
       integer :: beta, x(:), y(:,:)
     end subroutine psi_isctm
     subroutine psi_isctv(n,idx,x,beta,y)
       integer :: n, idx(:)
       integer :: beta, x(:), y(:)
     end subroutine psi_isctv
  end interface

end module psi_mod
