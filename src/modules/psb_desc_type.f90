!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!
!	Module to   define desc_a,
!      structure for coomunications.
!
! Typedef: psb_desc_type
!    Defines a communication descriptor


module psb_descriptor_type
  use psb_const_mod
  implicit none

  ! desc_type contains data for communications.
  type psb_desc_type
    ! contain decomposition informations
    integer, pointer :: matrix_data(:)=>null()
    ! contain index of halo elements to send/receive
    integer, pointer :: halo_index(:)=>null()
    ! contain indices of boundary  elements 
    integer, pointer :: bnd_elem(:)=>null()
    ! contain index of overlap elements to send/receive
    integer, pointer :: ovrlap_index(:)=>null()
    ! contain for each local overlap element, the number of times
    ! that is duplicated
    integer, pointer :: ovrlap_elem(:)=>null()
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
         & desc%halo_index,desc%bnd_elem,desc%ovrlap_elem,&
         & desc%ovrlap_index, desc%lprm, desc%idx_space)!,&
    !         & desc%halo_pt,desc%ovrlap_pt)

  end subroutine psb_nullify_desc

  logical function psb_is_ok_dec(dectype)
    integer :: dectype

    psb_is_ok_dec = ((dectype == psb_desc_asb_).or.(dectype == psb_desc_bld_).or.&
         & (dectype == psb_desc_upd_).or.(dectype== psb_desc_upd_asb_).or.&
         & (dectype== psb_desc_repl_))

  end function psb_is_ok_dec

  logical function psb_is_bld_dec(dectype)
    integer :: dectype

    psb_is_bld_dec = (dectype == psb_desc_bld_)
  end function psb_is_bld_dec

  logical function psb_is_upd_dec(dectype)          
    integer :: dectype

    psb_is_upd_dec = (dectype == psb_desc_upd_)

  end function psb_is_upd_dec

  logical function psb_is_asb_upd_dec(dectype)
    integer :: dectype

    psb_is_asb_upd_dec = (dectype == psb_desc_upd_asb_)

  end function psb_is_asb_upd_dec

  logical function psb_is_asb_dec(dectype)
    integer :: dectype

    psb_is_asb_dec = (dectype == psb_desc_asb_).or.&
         & (dectype== psb_desc_repl_)

  end function psb_is_asb_dec

end module psb_descriptor_type
