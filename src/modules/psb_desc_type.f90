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
     integer, allocatable :: matrix_data(:)
     ! contain index of halo elements to send/receive
     integer, allocatable :: halo_index(:)
     ! contain indices of boundary  elements 
     integer, allocatable :: bnd_elem(:)
     ! contain index of overlap elements to send/receive
     integer, allocatable :: ovrlap_index(:)
     ! contain for each local overlap element, the number of times
     ! that is duplicated
     integer, allocatable :: ovrlap_elem(:)
     ! contain for each local element the corresponding global index
     integer, allocatable :: loc_to_glob(:)
     ! contain for each global element the corresponding local index,
     ! if exist.
     integer, allocatable :: glob_to_loc (:)
     ! local renumbering induced by sparse matrix storage. 
     integer, allocatable :: lprm(:)
     ! index space in case it is not just the contiguous range 1:n
     integer, allocatable :: idx_space(:)
  end type psb_desc_type

contains 

  subroutine psb_nullify_desc(desc)
    type(psb_desc_type), intent(inout) :: desc
    
!!$    nullify(desc%matrix_data,desc%loc_to_glob,desc%glob_to_loc,&
!!$         &desc%halo_index,desc%bnd_elem,desc%ovrlap_elem,&
!!$         &desc%ovrlap_index, desc%lprm, desc%idx_space)

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

  integer function psb_get_local_rows(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_get_local_rows = desc%matrix_data(psb_n_row_)
  end function psb_get_local_rows

  integer function psb_get_local_cols(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_get_local_cols = desc%matrix_data(psb_n_col_)
  end function psb_get_local_cols

  integer function psb_get_global_rows(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_get_global_rows = desc%matrix_data(psb_m_)
  end function psb_get_global_rows

  integer function psb_get_global_cols(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_get_global_cols = desc%matrix_data(psb_n_)
  end function psb_get_global_cols

  integer function psb_get_context(desc)
    type(psb_desc_type), intent(in) :: desc
    
    psb_get_context = desc%matrix_data(psb_ctxt_)
  end function psb_get_context
    
end module psb_descriptor_type
