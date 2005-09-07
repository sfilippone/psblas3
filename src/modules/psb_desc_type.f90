!
!	Module to   define desc_a,
!      structure for coomunications.
!
! Typedef: psb_desc_type
!    Defines a communication descriptor


module psb_descriptor_type
  use psb_const_mod

  ! desc_type contains data for communications.
  type psb_desc_type
     ! contain decomposition informations
     integer, pointer :: matrix_data(:)=>null()
     ! contain index of halo elements to send/receive
     integer, pointer :: halo_index(:)=>null()
     ! contain indices of boundary  elements 
     integer, pointer :: bnd_elem(:)=>null()
     ! contain index of overlap elements to send/receive
     integer, pointer :: ovrlap_elem(:)=>null()
     ! contain for each local overlap element, the number of times
     ! that is duplicated
     integer, pointer :: ovrlap_index(:)=>null()
     ! contain for each local element the corresponding global index
     integer, pointer :: loc_to_glob(:)=>null()
     ! contain for each global element the corresponding local index,
     ! if exist.
     integer, pointer :: glob_to_loc (:)=>null()
     ! local renumbering induced by sparse matrix storage. 
     integer, pointer :: lprm(:)=>null()
     ! index space in case it is not just the contiguous range 1:n
     integer, pointer :: idx_space(:)=>null()
  end type psb_desc_type
  
contains 

  subroutine psb_nullify_desc(desc)
    type(psb_desc_type), intent(inout) :: desc
    
    nullify(desc%matrix_data,desc%loc_to_glob,desc%glob_to_loc,&
         &desc%halo_index,desc%bnd_elem,desc%ovrlap_elem,&
         &desc%ovrlap_index, desc%lprm, desc%idx_space)

  end subroutine psb_nullify_desc

  logical function psb_is_ok_dec(dectype)
    integer :: dectype
    
    psb_is_ok_dec = ((dectype == desc_asb).or.(dectype == desc_bld).or.&
         &(dectype == desc_upd).or.(dectype== desc_upd_asb))
    
  end function psb_is_ok_dec

  logical function psb_is_bld_dec(dectype)
    integer :: dectype

    psb_is_bld_dec = (dectype == desc_bld)
  end function psb_is_bld_dec

  logical function psb_is_upd_dec(dectype)          
    integer :: dectype

    psb_is_upd_dec = (dectype == desc_upd)

  end function psb_is_upd_dec

  logical function psb_is_asb_upd_dec(dectype)
    integer :: dectype

    psb_is_asb_upd_dec = (dectype == desc_upd_asb)

  end function psb_is_asb_upd_dec

  logical function psb_is_asb_dec(dectype)
    integer :: dectype

    psb_is_asb_dec = (dectype == desc_asb)

  end function psb_is_asb_dec

end module psb_descriptor_type
