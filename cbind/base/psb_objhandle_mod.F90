module psb_objhandle_mod
  use iso_c_binding

#if defined(LONG_INTEGERS)
  integer, parameter :: psb_c_int = c_int64_t
#else
  integer, parameter :: psb_c_int = c_int32_t
#endif

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

end module psb_objhandle_mod
