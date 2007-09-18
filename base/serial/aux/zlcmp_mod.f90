module zlcmp_mod
  interface operator(<)
    module procedure zllt
  end interface
  interface operator(<=)
    module procedure zlle
  end interface
  interface operator(>)
    module procedure zlgt
  end interface
  interface operator(>=)
    module procedure zlge
  end interface

contains

  function zllt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zllt
    
    zllt = (real(a)<real(b)).or.((real(a)==real(b)).and.(aimag(a)<aimag(b)))
  end function zllt
  function zlle(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zlle
    
    zlle = (real(a)<real(b)).or.((real(a)==real(b)).and.(aimag(a)<=aimag(b)))
  end function zlle

  function zlgt(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zlgt
    
    zlgt = (real(a)>real(b)).or.((real(a)==real(b)).and.(aimag(a)>aimag(b)))
  end function zlgt
  function zlge(a,b)
    complex(kind(1.d0)), intent(in) :: a,b
    logical :: zlge
    
    zlge = (real(a)>real(b)).or.((real(a)==real(b)).and.(aimag(a)>=aimag(b)))
  end function zlge

end module zlcmp_mod
  
