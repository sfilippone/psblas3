!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
module psb_error_mod

  integer, parameter, public :: psb_act_ret_=0, psb_act_abort_=1, psb_no_err_=0
  integer, parameter, public :: psb_debug_ext_=1, psb_debug_outer_=2
  integer, parameter, public :: psb_debug_comp_=3, psb_debug_inner_=4
  integer, parameter, public :: psb_debug_serial_=8, psb_debug_serial_comp_=9
  !
  !     Error handling 
  !
  public psb_errpush, psb_error, psb_get_errstatus,&
       & psb_get_errverbosity, psb_set_errverbosity,psb_errcomm, &
       & psb_erractionsave, psb_erractionrestore, &
       & psb_get_erraction, psb_set_erraction, &
       & psb_get_debug_level, psb_set_debug_level,&
       & psb_get_debug_unit, psb_set_debug_unit,&
       & psb_get_serial_debug_level, psb_set_serial_debug_level

  interface psb_error
    module procedure psb_serror
    module procedure psb_perror
  end interface


  private

  type psb_errstack_node

    !  the error code
    integer                  ::   err_code=0         
    !  the name of the routine generating the error
    character(len=20)        ::   routine=''       
    !  array of integer data to complete the error msg   
    integer,dimension(5)     ::   i_err_data=0     
    !     real(psb_dpk_)(dim=10) ::   r_err_data=0.d0    !  array of real data to complete the error msg
    !     complex(dim=10)          ::   c_err_data=0.c0    !  array of complex data to complete the error msg
    !  array of character data to complete the error msg
    character(len=40)        ::   a_err_data=''      
    !  pointer to the next element in the stack 
    type(psb_errstack_node), pointer :: next         

  end type psb_errstack_node


  type psb_errstack

    !  pointer to the top element of the stack
    type(psb_errstack_node), pointer :: top => null()    
    !  number of entries in the stack
    integer                          :: n_elems=0        

  end type psb_errstack


  type(psb_errstack), save :: error_stack         !  the PSBLAS-2.0 error stack
  integer, save            :: error_status=0      !  the error status (maybe not here)
  integer, save            :: verbosity_level=1   !  the verbosity level (maybe not here)
  integer, save            :: err_action=psb_act_abort_
  integer, save            :: debug_level=0, debug_unit=0, serial_debug_level=0

contains


  ! saves action to support error traceback
  ! also changes error action to "return"
  subroutine psb_erractionsave(err_act)
    integer, intent(out) :: err_act
    err_act    = err_action
    err_action = psb_act_ret_
  end subroutine psb_erractionsave


  ! return the action to take upon error occurrence
  subroutine psb_get_erraction(err_act)
    integer, intent(out) :: err_act
    err_act=err_action
  end subroutine psb_get_erraction

  ! sets the action to take upon error occurrence
  subroutine psb_set_erraction(err_act)
    integer, intent(in) :: err_act
    err_action=err_act
  end subroutine psb_set_erraction


  ! restores error action previously saved with psb_erractionsave
  subroutine psb_erractionrestore(err_act)
    integer, intent(in) :: err_act
    err_action=err_act
  end subroutine psb_erractionrestore


  function  psb_get_debug_level()
    integer :: psb_get_debug_level
    psb_get_debug_level = debug_level
  end function psb_get_debug_level

  subroutine psb_set_debug_level(level)
    integer, intent(in) :: level
    if (level >= 0) then
      debug_level = level
    else
      debug_level = 0
    end if
  end subroutine psb_set_debug_level

  function  psb_get_debug_unit()
    integer :: psb_get_debug_unit
    psb_get_debug_unit = debug_unit
  end function psb_get_debug_unit

  subroutine psb_set_debug_unit(unit)
    integer, intent(in) :: unit
    if (unit >= 0) then
      debug_unit = unit
    else
      debug_unit = 0
    end if
  end subroutine psb_set_debug_unit

  function  psb_get_serial_debug_level()
    integer :: psb_get_serial_debug_level
    psb_get_serial_debug_level = serial_debug_level
  end function psb_get_serial_debug_level

  subroutine psb_set_serial_debug_level(level)
    integer, intent(in) :: level
    if (level >= 0) then
      serial_debug_level = level
    else
      serial_debug_level = 0
    end if
  end subroutine psb_set_serial_debug_level


  ! checks wether an error has occurred on one of the porecesses in the execution pool
  subroutine psb_errcomm(ictxt, err)
    integer, intent(in)   :: ictxt
    integer, intent(inout):: err
    integer :: temp(2)
    integer, parameter :: ione=1
    ! Cannot use psb_amx or otherwise we have a recursion in module usage
#if !defined(SERIAL_MPI)
    call igamx2d(ictxt, 'A', ' ', ione, ione, err, ione,&
         &temp ,temp,-ione ,-ione,-ione)
#endif    
  end subroutine psb_errcomm



  ! sets verbosity of the error message
  subroutine psb_set_errverbosity(v)
    integer, intent(in) :: v
    verbosity_level=v
  end subroutine psb_set_errverbosity



  ! returns verbosity of the error message
  function psb_get_errverbosity()
    integer :: psb_get_errverbosity

    psb_get_errverbosity=verbosity_level
  end function psb_get_errverbosity



  ! checks the status of the error condition
  function psb_get_errstatus()
    integer :: psb_get_errstatus
    psb_get_errstatus=error_status
  end function psb_get_errstatus



  ! pushes an error on the error stack
  subroutine psb_errpush(err_c, r_name, i_err, a_err)

    integer, intent(in)              ::  err_c
    character(len=*), intent(in)     ::  r_name
    character(len=*), optional       ::  a_err
    integer, optional                ::  i_err(5)

    type(psb_errstack_node), pointer     ::  new_node


    allocate(new_node)

    new_node%err_code   = err_c
    new_node%routine    = r_name
    if(present(i_err)) then
      new_node%i_err_data = i_err
    end if
    if(present(a_err)) then 
      new_node%a_err_data = a_err
    end if
    new_node%next       => error_stack%top
    error_stack%top     => new_node
    error_stack%n_elems = error_stack%n_elems+1
    if(error_status == 0) error_status=1
    nullify(new_node)

  end subroutine psb_errpush


  ! pops an error from the error stack
  subroutine psb_errpop(err_c, r_name, i_e_d, a_e_d)

    integer, intent(out)             ::  err_c
    character(len=20), intent(out)        ::  r_name
    character(len=40), intent(out)        ::  a_e_d
    integer, intent(out)             ::  i_e_d(5)

    type(psb_errstack_node), pointer     ::  old_node

    err_c      =  error_stack%top%err_code
    r_name     =  error_stack%top%routine
    i_e_d      =  error_stack%top%i_err_data
    a_e_d      =  error_stack%top%a_err_data

    old_node   => error_stack%top
    error_stack%top  => old_node%next
    error_stack%n_elems = error_stack%n_elems - 1
    if(error_stack%n_elems == 0) error_status=0

    deallocate(old_node)

  end subroutine psb_errpop



  ! handles the occurence of an error in a parallel routine
  subroutine psb_perror(ictxt)

    integer, intent(in)     ::  ictxt
    integer                 ::  err_c
    character(len=20)       ::  r_name
    character(len=40)       ::  a_e_d
    integer                 ::  i_e_d(5)
    integer                 ::  nprow, npcol, me, mypcol
    integer, parameter      ::  ione=1, izero=0

#if defined(SERIAL_MPI)
    me = -1
#else        
    call blacs_gridinfo(ictxt,nprow,npcol,me,mypcol)
#endif


    if(error_status > 0) then
      if(verbosity_level > 1) then

        do while (error_stack%n_elems > izero)
          write(0,'(50("="))')
          call psb_errpop(err_c, r_name, i_e_d, a_e_d)
          call psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)
          !            write(0,'(50("="))')
        end do
#if defined(SERIAL_MPI)
        stop 
#else        
        call blacs_abort(ictxt,-1)
#endif
      else

        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)
        do while (error_stack%n_elems > 0)
          call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        end do
#if defined(SERIAL_MPI)
        stop 
#else        
        call blacs_abort(ictxt,-1)
#endif
      end if
    end if

    if(error_status > izero) then
#if defined(SERIAL_MPI)
        stop 
#else        
      call blacs_abort(ictxt,err_c)
#endif
    end if


  end subroutine psb_perror


  ! handles the occurence of an error in a serial routine
  subroutine psb_serror()

    integer                 ::  err_c
    character(len=20)       ::  r_name
    character(len=40)       ::  a_e_d
    integer                 ::  i_e_d(5)
    integer, parameter      ::  ione=1, izero=0

    if(error_status > 0) then
      if(verbosity_level > 1) then

        do while (error_stack%n_elems > izero)
          write(0,'(50("="))')
          call psb_errpop(err_c, r_name, i_e_d, a_e_d)
          call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
          !            write(0,'(50("="))')
        end do

      else

        call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        call psb_errmsg(err_c, r_name, i_e_d, a_e_d)
        do while (error_stack%n_elems > 0)
          call psb_errpop(err_c, r_name, i_e_d, a_e_d)
        end do
      end if
    end if

  end subroutine psb_serror


  ! prints the error msg associated to a specific error code
  subroutine psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)

    integer, intent(in)              ::  err_c
    character(len=20), intent(in)    ::  r_name
    character(len=40), intent(in)    ::  a_e_d
    integer, intent(in)              ::  i_e_d(5)
    integer, optional                ::  me

    if(present(me)) then
      write(0,'("Process: ",i0,".  PSBLAS Error (",i0,") in subroutine: ",a20)')me,err_c,r_name
    else
      write(0,'("PSBLAS Error (",i0,") in subroutine: ",a20)')err_c,r_name
    end if


    select case (err_c)
    case(:0)
      write (0,'("error on calling sperror. err_c must be greater than 0")')
    case(2)
      write (0,'("pivot too small: ",i0,1x,a)')i_e_d(1),a_e_d
    case(3)
      write (0,'("Invalid number of ovr:",i0)')i_e_d(1)
    case(5)
      write (0,'("Invalid input")')

    case(10)
      write (0,'("input argument n. ",i0," cannot be less than 0")')i_e_d(1)
      write (0,'("current value is ",i0)')i_e_d(2)

    case(20)
      write (0,'("input argument n. ",i0," cannot be greater than 0")')i_e_d(1)
      write (0,'("current value is ",i0)')i_e_d(2)
    case(30)
      write (0,'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write (0,'("current value is ",i0)')i_e_d(2)
    case(31)
      write (0,'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write (0,'("current value is ",a)')a_e_d
    case(35)
      write (0,'("Size of input array argument n. ",i0," is invalid.")')i_e_d(1)
      write (0,'("Current value is ",i0)')i_e_d(2)
    case(36)
      write (0,'("Size of input array argument n. ",i0," must be ")')i_e_d(1)
      write (0,'("at least  ",i0)')i_e_d(2)
    case(40)
      write (0,'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write (0,'("current value is ",a)')a_e_d(2:2)
    case(50)
      write (0,'("input argument n. ",i0," must be equal or greater than input argument n. ",i0)') i_e_d(1), i_e_d(3)
      write (0,'("current values are ",i0," < ",i0)') i_e_d(2),i_e_d(5)
    case(60)
      write (0,'("input argument n. ",i0," must be greater than or equal to ",i0)')i_e_d(1),i_e_d(2)
      write (0,'("current value is ",i0," < ",i0)')i_e_d(3), i_e_d(2)
    case(70)
      write (0,'("input argument n. ",i0," in entry # ",i0," has an invalid value")')i_e_d(1:2)
      write (0,'("current value is ",a)')a_e_d
    case(71)
      write (0,'("Impossible error in ASB: nrow>ncol,")')
      write (0,'("Actual values are ",i0," > ",i0)')i_e_d(1:2)
      !        ... csr format error ...
    case(80)
      write (0,'("input argument ia2(1) is less than 0")')
      write (0,'("current value is ",i0)')i_e_d(1)
      !        ... csr format error ...
    case(90)
      write (0,'("indices in ia2 array are not in  increasing order")')
    case(91)
      write (0,'("indices in ia1 array are not in increasing order")')
      !        ... csr format error ...
    case(100)
      write (0,'("indices in ia1 array are not within problem dimension")')
      write (0,'("problem dimension is ",i0)')i_e_d(1)
    case(110)
      write (0,'("invalid combination of input arguments")')
    case(115)
      write (0,'("Invalid process identifier in input array argument n. ",i0,".")')i_e_d(1)
      write (0,'("Current value is ",i0)')i_e_d(2)
    case(120)
      write (0,'("input argument n. ",i0," must be greater than input argument n. ",i0)')i_e_d(1:2)
      write (0,'("current values are ",i0," < ",i0)') i_e_d(3:4)
      !        ... coo format error ...
    case(130)
      write (0,'("there are duplicated elements in coo format")')
      write (0,'("and you have chosen psb_dupl_err_ ")')
    case(134)
      write (0,'("Invalid input format ",a3)')a_e_d(1:3)
    case(135)
      write (0,'("Format ",a3," not yet supported here")')a_e_d(1:3)
    case(136)
      write (0,'("Format ",a3," is unknown")')a_e_d(1:3)
    case(140)
      write (0,'("indices in input array are not within problem dimension ",2(i0,2x))')i_e_d(1:2)
    case(150)
      write (0,'("indices in input array are not belonging to the calling process ",i0)')i_e_d(1)
    case(290)
      write (0,'("To call this routine you must first call psb_geall on the same matrix")')
    case(295)
      write (0,'("To call this routine you must first call psb_spall on the same matrix")')
    case(300)
      write (0,'("Input argument n. ",i0," must be equal to entry n. ",i0," in array input argument n.",i0)') &
           & i_e_d(1),i_e_d(4),i_e_d(3)
      write (0,'("Current values are ",i0," != ",i0)')i_e_d(2), i_e_d(5)
    case(400)
      write (0,'("MPI error:",i0)')i_e_d(1)
    case(550)
      write (0,'("Parameter n. ",i0," must be equal on all processes. ",i0)')i_e_d(1)
    case(551)
      write (0,'("Entry n. ",i0," out of ",i0," should be between 1 and ",i0," but is ",i0)')i_e_d(1),i_e_d(3),i_e_d(4),i_e_d(2)
    case(552)
      write (0,'("Index lists are inconsistent: some indices are orphans")')
    case(570)
      write (0,'("partition function passed as input argument n. ",i0," returns number of processes")')i_e_d(1)
      write (0,'("greater than No of grid s processes on global point ",i0,". Actual number of grid s ")')i_e_d(4)
      write (0,'("processes is ",i0,", number returned is ",i0)')i_e_d(2),i_e_d(3)
    case(575)
      write (0,'("partition function passed as input argument n. ",i0," returns number of processes")')i_e_d(1)
      write (0,'("less or equal to 0 on global point ",i0,". Number returned is ",i0)')i_e_d(3),i_e_d(2)
    case(580)
      write (0,'("partition function passed as input argument n. ",i0," returns wrong processes identifier")')i_e_d(1)
      write (0,'("on global point ",i0,". Current value returned is : ",i0)')i_e_d(3),i_e_d(2)
    case(581)
      write (0,'("Exactly one of the optional arguments  ",a," must be present")')a_e_d
    case(582)
      write (0,'("Argument M is required when argument PARTS is specified")')
    case(583)
      write (0,'("No more than one of the optional arguments  ",a," must be present")')a_e_d
    case(600)
      write (0,'("Sparse Matrix and descriptors are in an invalid state for this subroutine call: ",i0)')i_e_d(1)
    case(700)
      write (0,'("The base version of subroutine ''",a,"'' has been called.",/,&
           &"The class implementation for ''",a,"'' may be incomplete!")') &
           & trim(r_name), trim(a_e_d)

    case (1121)
      write (0,'("Invalid state for sparse matrix A")')
    case (1122)
      write (0,'("Invalid state for communication descriptor")')
    case (1123)
      write (0,'("Invalid combined state for A and DESC_A")')
    case (1124)
      write (0,'("Invalid state for object:",a)') trim(a_e_d)
    case(1125:1999)
      write (0,'("computational error. code: ",i0)')err_c
    case(2010)
      write (0,'("BLACS error. Number of processes=-1")')
    case(2011)
      write (0,'("Initialization error: not enough processes available in the parallel environment")')
    case(2030)
      write (0,'("BLACS ERROR: Number of grid columns must be equal to 1\nCurrent value is ",i4," != 1.")')i_e_d(1)
    case(2231)
      write (0,'("Invalid input state for matrix.")')
    case(2232)
      write (0,'("Input state for matrix is not adequate for regeneration.")')
    case (2233:2999)
      write(0,'("resource error. code: ",i0)')err_c
    case(3000:3009)
      write (0,'("sparse matrix representation ",a3," not yet implemented")')a_e_d(1:3)
    case(3010)
      write (0,'("Case lld not equal matrix_data[N_COL_] is not yet implemented.")')
    case(3015)
      write (0,'("transpose option for sparse matrix representation ",a3," not implemented")')a_e_d(1:3)
    case(3020)
      write (0,'("Case trans = C is not yet implemented.")') 
    case(3021)
      write (0,'("Case trans /= N is not yet implemented.")') 
    case(3022)
      write (0,'("Only unit diagonal so far for triangular matrices. ")') 
    case(3023)
      write (0,'("Cases DESCRA(1:1)=S  DESCRA(1:1)=T not yet implemented. ")') 
    case(3024)
      write (0,'("Cases DESCRA(1:1)=G not yet implemented. ")') 
    case(3030)
      write (0,'("Case ja /= ix or ia/=iy is not yet implemented.")')
    case(3040)
      write (0,'("Case ix /= 1 or iy /= 1 is not yet implemented.")')
    case(3050)
      write (0,'("Case ix /= iy is not yet implemented.")')
    case(3060)
      write (0,'("Case ix /= 1 is not yet implemented.")')
    case(3070)
      write (0,'("This operation is only implemented with no overlap.")')
    case(3080)
      write (0,'("Decompostion type ",i0," not yet supported.")')i_e_d(1)
    case(3090)
      write (0,'("Insert matrix mode not yet implemented.")')
    case(3100)
      write (0,'("Error on index. Element has not been inserted")')
      write (0,'("local index is: ",i0," and global index is:",i0)')i_e_d(1:2)
    case(3110)
      write (0,'("Before you call this routine, you must assembly sparse matrix")')
    case(3111)
      write (0,'("Before you call this routine, you must initialize the preconditioner")')
    case(3112)
      write (0,'("Before you call this routine, you must build the preconditioner")')
    case(3113:3999)
      write(0,'("miscellaneus error. code: ",i0)')err_c
    case(4000)
      write(0,'("Allocation/deallocation error")')
    case(4001)
      write(0,'("Internal error: ",a)')a_e_d
    case(4010)
      write (0,'("Error from call to subroutine ",a)')a_e_d
    case(4011)
      write (0,'("Error from call to a subroutine ")')
    case(4012)
      write (0,'("Error ",i0," from call to a subroutine ")')i_e_d(1)
    case(4013)
      write (0,'("Error from call to subroutine ",a," ",i0)')a_e_d,i_e_d(1)
    case(4025)
      write (0,'("Error on allocation request for ",i0," items of type ",a)')i_e_d(1),a_e_d
    case(4110)
      write (0,'("Error ",i0," from call to an external package in subroutine ",a)')i_e_d(1),a_e_d
    case (5001)
      write (0,'("Invalid ISTOP: ",i0)')i_e_d(1)
    case (5002)
      write (0,'("Invalid PREC: ",i0)')i_e_d(1)
    case (5003)
      write (0,'("Invalid PREC: ",a3)')a_e_d(1:3)
    case default
      write(0,'("unknown error (",i0,") in subroutine ",a)')err_c,r_name
      write(0,'(5(i0,2x))') i_e_d
      write(0,'(a)') a_e_d

    end select

  end subroutine psb_errmsg



end module psb_error_mod
