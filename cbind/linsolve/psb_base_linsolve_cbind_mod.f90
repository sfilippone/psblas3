module psb_base_linsolve_cbind_mod

  use iso_c_binding
  use psb_objhandle_mod
  
  type, bind(c) :: solveroptions
    integer(psb_c_ipk_) :: iter, itmax, itrace, irst, istop
    real(c_double) ::  eps, err
    type(psb_c_object_type)    ::  s1, s2
  end type solveroptions

contains 

  function psb_c_DefaultSolverOptions(options)&
       &   bind(c,name='psb_c_DefaultSolverOptions') result(res)
    implicit none 
    type(solveroptions)   :: options
    integer(psb_c_ipk_) :: res
    options%itmax  = 1000
    options%itrace = 0
    options%istop  = 2
    options%irst   = 10
    options%eps    = 1.d-6
    options%s1     = psb_c_get_new_object()
    options%s2     = psb_c_get_new_object()
    res            = 0
  end function psb_c_DefaultSolverOptions

  function psb_c_PrintSolverOptions(options)&
       &   bind(c,name='psb_c_PrintSolverOptions') result(res)
    implicit none 
    type(solveroptions), value   :: options
    integer(psb_c_ipk_) :: res

    write(*,*) 'PSBLAS C Interface Solver Options '
    write(*,*) '   Maximum number of iterations  :', options%itmax
    write(*,*) '   Tracing                       :', options%itrace
    write(*,*) '   Stopping Criterion            :', options%istop
    write(*,*) '   Restart                       :', options%irst
    write(*,*) '   EPS (tolerance)               :', options%eps  
    write(*,*) '   S1 scaling                    :', c_associated(options%s1%item)
    write(*,*) '   S2 scaling                    :', c_associated(options%s2%item)
    res            = 0
  end function psb_c_PrintSolverOptions
    

end module psb_base_linsolve_cbind_mod
