module psb_renum_mod
  use psb_base_mod

  integer, parameter :: psb_mat_renum_identity_ = 0
  integer, parameter :: psb_mat_renum_gps_ = 456
  integer, parameter :: psb_mat_renum_amd_ = psb_mat_renum_gps_ + 1 
  
  
  interface psb_mat_renum
    subroutine psb_d_mat_renums(alg,mat,info,perm)
      import  psb_dspmat_type
      character(len=*), intent(in) :: alg
      type(psb_dspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_d_mat_renums
    subroutine psb_d_mat_renum(alg,mat,info,perm)
      import  psb_dspmat_type
      integer, intent(in) :: alg
      type(psb_dspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_d_mat_renum
    subroutine psb_s_mat_renums(alg,mat,info,perm)
      import  psb_sspmat_type
      character(len=*), intent(in) :: alg
      type(psb_sspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_s_mat_renums
    subroutine psb_s_mat_renum(alg,mat,info,perm)
      import  psb_sspmat_type
      integer, intent(in) :: alg
      type(psb_sspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_s_mat_renum
    subroutine psb_z_mat_renums(alg,mat,info,perm)
      import  psb_zspmat_type
      character(len=*), intent(in) :: alg
      type(psb_zspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_z_mat_renums
    subroutine psb_z_mat_renum(alg,mat,info,perm)
      import  psb_zspmat_type
      integer, intent(in) :: alg
      type(psb_zspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_z_mat_renum
    subroutine psb_c_mat_renums(alg,mat,info,perm)
      import  psb_cspmat_type
      character(len=*), intent(in) :: alg
      type(psb_cspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_c_mat_renums
    subroutine psb_c_mat_renum(alg,mat,info,perm)
      import  psb_cspmat_type
      integer, intent(in) :: alg
      type(psb_cspmat_type), intent(inout) :: mat
      integer, intent(out) :: info
      integer, allocatable, optional, intent(out) :: perm(:)
    end subroutine psb_c_mat_renum
  end interface psb_mat_renum

  
  interface psb_cmp_bwpf
    subroutine psb_s_cmp_bwpf(mat,bwl,bwu,prf,info)
      import  psb_sspmat_type
      type(psb_sspmat_type), intent(in) :: mat
      integer, intent(out) :: bwl, bwu
      integer, intent(out) :: prf
      integer, intent(out) :: info
    end subroutine psb_s_cmp_bwpf
    subroutine psb_d_cmp_bwpf(mat,bwl,bwu,prf,info)
      import  psb_dspmat_type
      type(psb_dspmat_type), intent(in) :: mat
      integer, intent(out) :: bwl, bwu
      integer, intent(out) :: prf
      integer, intent(out) :: info
    end subroutine psb_d_cmp_bwpf
    subroutine psb_c_cmp_bwpf(mat,bwl,bwu,prf,info)
      import  psb_cspmat_type
      type(psb_cspmat_type), intent(in) :: mat
      integer, intent(out) :: bwl, bwu
      integer, intent(out) :: prf
      integer, intent(out) :: info
    end subroutine psb_c_cmp_bwpf
    subroutine psb_z_cmp_bwpf(mat,bwl,bwu,prf,info)
      import  psb_zspmat_type
      type(psb_zspmat_type), intent(in) :: mat
      integer, intent(out) :: bwl, bwu
      integer, intent(out) :: prf
      integer, intent(out) :: info
    end subroutine psb_z_cmp_bwpf
  end interface psb_cmp_bwpf


end module psb_renum_mod
