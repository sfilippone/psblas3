module zacmp_mod
  interface operator(<)
    module procedure zalt
  end interface
  interface operator(<=)
    module procedure zale
  end interface
  interface operator(>)
    module procedure zagt
  end interface
  interface operator(>=)
    module procedure zage
  end interface

contains

  function zalt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zalt
    
    zalt = (abs(a) < abs(b))
  end function zalt
  function zale(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zale
    
    zale = (abs(a) <= abs(b))
  end function zale

  function zagt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zagt
    
    zagt = (abs(a) > abs(b))
  end function zagt
  function zage(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zage
    
    zage = (abs(a) >= abs(b))
  end function zage

end module zacmp_mod
  
