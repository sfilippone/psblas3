module psb_objhandle_mod
  use iso_c_binding

  integer, parameter :: psb_c_mpk = c_int32_t
#if defined(IPK4) &&  defined(LPK4)
  integer, parameter :: psb_c_ipk = c_int32_t
  integer, parameter :: psb_c_lpk = c_int32_t
#elif defined(IPK4) &&  defined(LPK8)
  integer, parameter :: psb_c_ipk = c_int32_t
  integer, parameter :: psb_c_lpk = c_int64_t
#elif defined(IPK8) &&  defined(LPK8)
  integer, parameter :: psb_c_ipk = c_int64_t
  integer, parameter :: psb_c_lpk = c_int64_t
#else
  integer, parameter :: psb_c_ipk = -1
  integer, parameter :: psb_c_lpk = -1
#endif
  integer, parameter :: psb_c_epk = c_int64_t

  
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
