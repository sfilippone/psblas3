module psb_objhandle_mod
  use iso_c_binding
  use psb_cbind_const_mod
  
  type, bind(c) :: psb_c_object_type
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_object_type

  type, bind(c) :: psb_c_descriptor
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_descriptor

  type, bind(c) :: psb_c_svector
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_svector

  type, bind(c) :: psb_c_sspmat
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_sspmat

  type, bind(c) :: psb_c_dvector
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_dvector

  type, bind(c) :: psb_c_dspmat
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_dspmat

  type, bind(c) :: psb_c_cvector
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_cvector

  type, bind(c) :: psb_c_cspmat
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_cspmat

  type, bind(c) :: psb_c_zvector
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_zvector

  type, bind(c) :: psb_c_zspmat
    type(c_ptr) :: item = c_null_ptr
  end type psb_c_zspmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sp3mm c code structs
! TODO : rename to conventions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type, bind(C, name=spmat) :: spmat_t
    ! number of non zeros and dimensions
    integer(c_size_t)               :: nz, m, n
    ! value array
    real(c_float), allocatable      :: as(:)
    ! columns array
    integer(c_size_t), allocatable  :: ja(:)
    ! row index pointers array
    integer(c_size_t), allocatable  :: irp(:)
    ! lengths of the rows
    integer(c_size_t), allocatable  :: rl(:)
    ! max value of rl
    integer(c_size_t)               :: max_row_nz
  end type spmat_t

  type, bind(C, name=CONFIG) :: config_t
    ! dimensions of the grid
    integer(c_short)  :: grid_rows, grid_cols
    ! how to compute symb mul (if required)
    integer(c_int)    :: symb_mm_row_impl_id
    ! thread num to use in OMP parallel region
    integer(c_int)    :: thread_num
    ! CHUNKS_DISTR_INTERF func pntr
    type(c_ptr)       :: chunk_distrb_func
  end type config_t
end module psb_objhandle_mod
