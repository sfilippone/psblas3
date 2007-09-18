module zalcmp_mod
  interface operator(<)
    module procedure zallt
  end interface
  interface operator(<=)
    module procedure zalle
  end interface
  interface operator(>)
    module procedure zalgt
  end interface
  interface operator(>=)
    module procedure zalge
  end interface

contains

  function zallt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zallt
    
    zallt = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a))==abs(real(b))).and.(abs(aimag(a))<abs(aimag(b))))
  end function zallt
  function zalle(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zalle
    
    zalle = (abs(real(a))<abs(real(b))).or. &
         & ((abs(real(a))==abs(real(b))).and.(abs(aimag(a))<=abs(aimag(b))))
  end function zalle

  function zalgt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zalgt
    
    zalgt = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a))==abs(real(b))).and.(abs(aimag(a))>abs(aimag(b))))
  end function zalgt
  function zalge(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zalge
    
    zalge = (abs(real(a))>abs(real(b))).or. &
         & ((abs(real(a))==abs(real(b))).and.(abs(aimag(a))>=abs(aimag(b))))
  end function zalge

end module zalcmp_mod
  
