module psb_base_krylov_cbind_mod

  use iso_c_binding
  type, bind(c) :: solveroptions
    integer(c_int) :: iter, itmax, itrace, irst, istop
    real(c_double) ::  eps, err
  end type solveroptions

contains 

  function psb_c_DefaultSolverOptions(options)&
       &   bind(c,name='psb_c_DefaultSolverOptions') result(res)
    implicit none 
    type(solveroptions)   :: options
    integer(c_int) :: res

    options%itmax  = 1000
    options%itrace = 0
    options%istop  = 2
    options%irst   = 10
    options%eps    = 1.d-6

    res            = 0
  end function psb_c_DefaultSolverOptions
    

end module psb_base_krylov_cbind_mod
