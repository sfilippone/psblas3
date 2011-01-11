module elldev_mod
  use iso_c_binding 

  integer(c_int), parameter :: elldev_float       = 0
  integer(c_int), parameter :: elldev_double      = 1 
  integer(c_int), parameter :: elldev_success     = 0
  integer(c_int), parameter :: elldev_nomem       = -1
  integer(c_int), parameter :: elldev_unsupported = -2
  
  type, bind(c) :: elldev_parms
    integer(c_int) :: element_type
    integer(c_int) :: pitch
    integer(c_int) :: rows
    integer(c_int) :: maxRowSize
    integer(c_int) :: firstIndex
  end type elldev_parms

#ifdef HAVE_ELL_GPU  
  interface 
    function getEllDeviceParams(rows,maxRowSize,elementType,firstIndex) &
         & result(val) bind(c,name='getEllDeviceParams')
      use iso_c_binding
      import :: elldev_parms
      type(elldev_parms)    :: val
      integer(c_int), value :: rows,maxRowSize,elementType,firstIndex
    end function getEllDeviceParams
  end interface

  interface 
    function allocEllDevice(deviceMat,parms) &
         & result(val) bind(c,name='allocEllDevice')
      use iso_c_binding
      import :: elldev_parms
      integer(c_int)            :: val
      type(elldev_parms), value :: parms
      type(c_ptr)               :: deviceMat
    end function allocEllDevice
  end interface

  interface 
    subroutine  freeEllDevice(deviceMat) &
         & bind(c,name='freeEllDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeEllDevice
  end interface

  interface 
    function allocVecDevice(deviceVec, size, datatype) &
         & result(val) bind(c,name='allocVecDevice')
      use iso_c_binding
      import :: elldev_parms
      integer(c_int)            :: val
      integer(c_int), value     :: size, datatype
      type(c_ptr)               :: deviceVec
    end function allocVecDevice
  end interface

  interface 
    subroutine  freeVecDevice(deviceVec) &
         & bind(c,name='freeVecDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
    end subroutine freeVecDevice
  end interface

  interface writeVecDevice
    function writeVecDeviceFloat(deviceVec,hostVec) &
         & result(val) bind(c,name='writeVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: val
      type(c_ptr), value  :: deviceVec
      real(c_float)       :: hostVec(*)
    end function writeVecDeviceFloat
    function writeVecDeviceDouble(deviceVec,hostVec) &
         & result(val) bind(c,name='writeVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: val
      type(c_ptr), value  :: deviceVec
      real(c_double)      :: hostVec(*)
    end function writeVecDeviceDouble
  end interface writeVecDevice

  interface readVecDevice
    function readVecDeviceFloat(deviceVec,hostVec) &
         & result(val) bind(c,name='readVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: val
      type(c_ptr), value  :: deviceVec
      real(c_float)      :: hostVec(*)
    end function readVecDeviceFloat
    function readVecDeviceDouble(deviceVec,hostVec) &
         & result(val) bind(c,name='readVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: val
      type(c_ptr), value  :: deviceVec
      real(c_double)      :: hostVec(*)
    end function readVecDeviceDouble
  end interface readVecDevice

  interface spmvEllDevice
    function spmvEllDeviceFloat(deviceMat,alpha,x,beta,y) &
         & result(val) bind(c,name='spmvEllDeviceFloat')
      use iso_c_binding
      integer(c_int)     :: val
      type(c_ptr), value :: deviceMat, x, y
      real(c_float)     :: alpha, beta
    end function spmvEllDeviceFloat
    function spmvEllDeviceDouble(deviceMat,alpha,x,beta,y) &
         & result(val) bind(c,name='spmvEllDeviceDouble')
      use iso_c_binding
      integer(c_int)     :: val
      type(c_ptr), value :: deviceMat, x, y 
      real(c_double)     :: alpha,  beta
    end function spmvEllDeviceDouble
  end interface spmvEllDevice
    
  interface dotVecDevice
    function dotVecDeviceFloat(res,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: val
      real(c_float)       :: res
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function dotVecDeviceFloat
    function dotVecDeviceDouble(res,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: val
      real(c_double)      :: res
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function dotVecDeviceDouble
  end interface dotVecDevice
  
  interface axpbyVecDevice
    function axpbyVecDeviceFloat(alpha,deviceVecA,beta,deviceVecB) &
         & result(val) bind(c,name='axpbyVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: val
      real(c_float)       :: alpha, beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyVecDeviceFloat
    function axpbyVecDeviceDouble(alpha,deviceVecA,beta,deviceVecB) &
         & result(val) bind(c,name='axpbyVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: val
      real(c_double)      :: alpha,beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyVecDeviceDouble
  end interface axpbyVecDevice
  

#endif  


end module elldev_mod
