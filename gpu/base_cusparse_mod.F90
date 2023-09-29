!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
  

module base_cusparse_mod
  use iso_c_binding 
  ! Interface to CUSPARSE. 

  enum, bind(c)
    enumerator cusparse_status_success
    enumerator cusparse_status_not_initialized
    enumerator cusparse_status_alloc_failed
    enumerator cusparse_status_invalid_value
    enumerator cusparse_status_arch_mismatch
    enumerator cusparse_status_mapping_error
    enumerator cusparse_status_execution_failed
    enumerator cusparse_status_internal_error
    enumerator cusparse_status_matrix_type_not_supported
  end enum

  enum, bind(c)
    enumerator cusparse_matrix_type_general 
    enumerator cusparse_matrix_type_symmetric     
    enumerator cusparse_matrix_type_hermitian 
    enumerator cusparse_matrix_type_triangular 
  end enum
  
  enum, bind(c)
    enumerator cusparse_fill_mode_lower 
    enumerator cusparse_fill_mode_upper
  end enum
  
  enum, bind(c)
    enumerator cusparse_diag_type_non_unit 
    enumerator cusparse_diag_type_unit
  end enum
  
  enum, bind(c)
    enumerator cusparse_index_base_zero 
    enumerator cusparse_index_base_one
  end enum
  
  enum, bind(c)
    enumerator cusparse_operation_non_transpose  
    enumerator cusparse_operation_transpose
    enumerator cusparse_operation_conjugate_transpose
  end enum
  
  enum, bind(c)
    enumerator cusparse_direction_row
    enumerator cusparse_direction_column
  end enum


#if defined(HAVE_CUDA) && defined(HAVE_SPGPU)
  
  interface 
    function FcusparseCreate() &
         & bind(c,name="FcusparseCreate") result(res)
      use iso_c_binding
      integer(c_int) :: res
    end function FcusparseCreate
  end interface

  interface 
    function FcusparseDestroy() &
         & bind(c,name="FcusparseDestroy") result(res)
      use iso_c_binding
      integer(c_int) :: res
    end function FcusparseDestroy
  end interface

contains
  
  function initFcusparse() result(res)
    implicit none 
    integer(c_int) :: res
    
    res = FcusparseCreate()
  end function initFcusparse

  function closeFcusparse() result(res)
    implicit none 
    integer(c_int) :: res
    res = FcusparseDestroy()
  end function closeFcusparse

#endif
end module base_cusparse_mod
