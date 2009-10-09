!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
module psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_sort_mod

  use psi_serial_mod, &
       & psb_gth => psi_gth,&
       & psb_sct => psi_sct

  use psb_mat_mod

  interface psb_symbmm
    subroutine psb_dsymbmm(a,b,c,info)
      use psb_mat_mod
      implicit none 
      type(psb_d_sparse_mat), intent(in)  :: a,b
      type(psb_d_sparse_mat), intent(out) :: c
      integer, intent(out)                :: info
    end subroutine psb_dsymbmm
    subroutine psb_dbase_symbmm(a,b,c,info)
      use psb_mat_mod
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(out)  :: c
      integer, intent(out)                     :: info
    end subroutine psb_dbase_symbmm
  end interface
  interface psb_numbmm
    subroutine psb_dnumbmm(a,b,c)
      use psb_mat_mod
      implicit none 
      type(psb_d_sparse_mat), intent(in) :: a,b
      type(psb_d_sparse_mat), intent(inout)  :: c
    end subroutine psb_dnumbmm
    subroutine psb_dbase_numbmm(a,b,c)
      use psb_mat_mod
      implicit none 
      class(psb_d_base_sparse_mat), intent(in) :: a,b
      type(psb_d_csr_sparse_mat), intent(inout)  :: c
    end subroutine psb_dbase_numbmm
  end interface

  interface psb_rwextd
    subroutine psb_drwextd(nr,a,info,b,rowscale)
      use psb_mat_mod
      implicit none
      integer, intent(in)                          :: nr
      type(psb_d_sparse_mat), intent(inout)        :: a
      integer,intent(out)                          :: info
      type(psb_d_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                 :: rowscale
    end subroutine psb_drwextd
    subroutine psb_dbase_rwextd(nr,a,info,b,rowscale)
      use psb_mat_mod
      implicit none
      integer, intent(in)                                :: nr
      class(psb_d_base_sparse_mat), intent(inout)        :: a
      integer,intent(out)                                :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: b
      logical,intent(in), optional                       :: rowscale
    end subroutine psb_dbase_rwextd
  end interface
  
  
end module psb_serial_mod

