!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
module psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_sort_mod

  use psi_serial_mod, &
       & psb_gth => psi_gth,&
       & psb_sct => psi_sct


  interface psb_symbmm
    subroutine psb_ssymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none 
      type(psb_sspmat_type), intent(in)  :: a,b
      type(psb_sspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_ssymbmm
    subroutine psb_sbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_sbase_symbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none 
      type(psb_dspmat_type), intent(in)  :: a,b
      type(psb_dspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_dsymbmm
    subroutine psb_dbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_dbase_symbmm
    subroutine psb_csymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none 
      type(psb_cspmat_type), intent(in)  :: a,b
      type(psb_cspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_csymbmm
    subroutine psb_cbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_cbase_symbmm
    subroutine psb_zsymbmm(a,b,c,info)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none 
      type(psb_zspmat_type), intent(in)  :: a,b
      type(psb_zspmat_type), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_zsymbmm
    subroutine psb_zbase_symbmm(a,b,c,info)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_zbase_symbmm
  end interface

  interface psb_numbmm
    subroutine psb_snumbmm(a,b,c)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none 
      type(psb_sspmat_type), intent(in) :: a,b
      type(psb_sspmat_type), intent(inout)  :: c
    end subroutine psb_snumbmm
    subroutine psb_sbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_s_base_sparse_mat, psb_s_csr_sparse_mat
      implicit none 
      class(psb_s_base_sparse_mat), intent(in) :: a,b
      type(psb_s_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_sbase_numbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none 
      type(psb_dspmat_type), intent(in) :: a,b
      type(psb_dspmat_type), intent(inout)  :: c
    end subroutine psb_dnumbmm
    subroutine psb_dbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_d_base_sparse_mat, psb_d_csr_sparse_mat
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_dbase_numbmm
    subroutine psb_cnumbmm(a,b,c)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none 
      type(psb_cspmat_type), intent(in) :: a,b
      type(psb_cspmat_type), intent(inout)  :: c
    end subroutine psb_cnumbmm
    subroutine psb_cbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_c_base_sparse_mat, psb_c_csr_sparse_mat
      implicit none 
      class(psb_c_base_sparse_mat), intent(in) :: a,b
      type(psb_c_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_cbase_numbmm
    subroutine psb_znumbmm(a,b,c)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none 
      type(psb_zspmat_type), intent(in) :: a,b
      type(psb_zspmat_type), intent(inout)  :: c
    end subroutine psb_znumbmm
    subroutine psb_zbase_numbmm(a,b,c)
      use psb_mat_mod, only : psb_z_base_sparse_mat, psb_z_csr_sparse_mat
      implicit none 
      class(psb_z_base_sparse_mat), intent(in) :: a,b
      type(psb_z_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_zbase_numbmm
  end interface

  interface psb_rwextd
    subroutine psb_srwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_sspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_sspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_sspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_srwextd
    subroutine psb_sbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_s_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_s_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_sbase_rwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_dspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_dspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_dspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_drwextd
    subroutine psb_dbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_d_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_d_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_dbase_rwextd
    subroutine psb_crwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_cspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_cspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_cspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_crwextd
    subroutine psb_cbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_c_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_c_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_c_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_cbase_rwextd
    subroutine psb_zrwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_zspmat_type
      implicit none
      integer, intent(in)                          :: nr
      type(psb_zspmat_type), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_zspmat_type), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_zrwextd
    subroutine psb_zbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod, only : psb_z_base_sparse_mat
      implicit none
      integer, intent(in)                                :: nr
      class(psb_z_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_zbase_rwextd
  end interface
  
  
end module psb_serial_mod

