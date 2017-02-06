!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
  use psb_const_mod

  integer(psb_ipk_), parameter, public :: psb_act_ret_=0
  integer(psb_ipk_), parameter, public :: psb_act_print_=1
  integer(psb_ipk_), parameter, public :: psb_act_abort_=2
  integer(psb_ipk_), parameter, public :: psb_debug_ext_=1, psb_debug_outer_=2
  integer(psb_ipk_), parameter, public :: psb_debug_comp_=3, psb_debug_inner_=4
  integer(psb_ipk_), parameter, public :: psb_debug_serial_=8, psb_debug_serial_comp_=9
  
  integer(psb_ipk_), parameter, public ::  psb_no_err_      = 0
  integer(psb_ipk_), parameter, public ::  psb_err_warning_ = 1
  integer(psb_ipk_), parameter, public ::  psb_err_fatal_   = 2

  integer(psb_ipk_), parameter, public ::  psb_max_errmsg_len_   = 132
  
  !
  !     Error handling 
  !
  public psb_errpush, psb_error, psb_get_errstatus,&
       & psb_errstatus_fatal, psb_errstatus_warning,&
       & psb_errstatus_ok, psb_warning_push,&
       & psb_errpop, psb_errcomm, psb_get_numerr, &
       & psb_get_errverbosity, psb_set_errverbosity, &
       & psb_erractionsave, psb_erractionrestore, &
       & psb_get_erraction, psb_set_erraction, &
       & psb_set_erract_return, psb_set_erract_print, psb_set_erract_abort,&
       & psb_is_erract_return, psb_is_erract_print, psb_is_erract_abort,&
       & psb_get_debug_level, psb_set_debug_level,&
       & psb_get_debug_unit, psb_set_debug_unit,&
       & psb_get_serial_debug_level, psb_set_serial_debug_level,&
       & psb_clean_errstack, psb_error_handler, &
       & psb_ser_error_handler, psb_par_error_handler, &
       & psb_ser_error_print_stack, psb_par_error_print_stack,&
       & psb_error_print_stack, psb_errmsg, psb_ach_errmsg


  interface psb_error_handler
    subroutine psb_ser_error_handler(err_act)
      import :: psb_ipk_
      integer(psb_ipk_), intent(inout) ::  err_act
    end subroutine psb_ser_error_handler
    subroutine psb_par_error_handler(ictxt,err_act)
      import :: psb_ipk_,psb_mpik_
      integer(psb_ipk_), intent(in) ::  ictxt
      integer(psb_ipk_), intent(in) ::  err_act
    end subroutine psb_par_error_handler
  end interface
 
  interface psb_error
    subroutine psb_serror()
    end subroutine psb_serror
    subroutine psb_perror(ictxt,abrt)
      import :: psb_mpik_
      integer(psb_mpik_), intent(in) ::  ictxt
      logical, intent(in), optional  :: abrt
    end subroutine psb_perror
  end interface


  interface psb_error_print_stack
    subroutine psb_par_error_print_stack(ictxt)
      import :: psb_ipk_,psb_mpik_
      integer(psb_mpik_), intent(in) ::  ictxt
    end subroutine psb_par_error_print_stack
    subroutine psb_ser_error_print_stack()
    end subroutine psb_ser_error_print_stack
  end interface

  interface psb_errcomm
    subroutine psb_errcomm(ictxt, err)
      import :: psb_mpik_, psb_ipk_
      integer(psb_mpik_), intent(in)   :: ictxt
      integer(psb_ipk_), intent(inout):: err
    end subroutine psb_errcomm
  end interface psb_errcomm

  interface psb_errpop
    module procedure psb_errpop, psb_ach_errpop
  end interface

  interface psb_errmsg
    module procedure psb_errmsg, psb_ach_errmsg
  end interface

#if defined(LONG_INTEGERS)
  interface psb_error
    module procedure psb_perror_ipk
  end interface psb_error
  interface psb_errcomm
    module procedure psb_errcomm_ipk
  end interface psb_errcomm
#endif


  private

  type psb_errstack_node

    !  the error code
    integer(psb_ipk_) ::   err_code=0         
    !  the name of the routine generating the error
    character(len=20)        ::   routine=''       
    !  array of integer data to complete the error msg   
    integer(psb_ipk_),dimension(5)     ::   i_err_data=0     
    !     real(psb_dpk_)(dim=10) ::   r_err_data=0.d0   
    !  array of real data to complete the error msg
    !     complex(dim=10)          ::   c_err_data=0.c0    
    !  array of complex data to complete the error msg
    !  array of character data to complete the error msg
    character(len=40)        ::   a_err_data=''      
    !  pointer to the next element in the stack 
    type(psb_errstack_node), pointer :: next         

  end type psb_errstack_node


  type psb_errstack

    !  pointer to the top element of the stack
    type(psb_errstack_node), pointer :: top => null()    
    !  number of entries in the stack
    integer(psb_ipk_) :: n_elems=0        

  end type psb_errstack


  type(psb_errstack), save  :: error_stack         
  integer(psb_ipk_), save   :: error_status    = psb_no_err_    
  integer(psb_ipk_), save   :: verbosity_level = 1 
  integer(psb_ipk_), save   :: err_action      = psb_act_abort_
  integer(psb_ipk_), save   :: debug_level     = 0, debug_unit, serial_debug_level=0

contains

#if defined(LONG_INTEGERS)
  subroutine psb_errcomm_ipk(ictxt, err)
    integer(psb_ipk_), intent(in)   :: ictxt
    integer(psb_ipk_), intent(inout):: err
    integer(psb_mpik_) :: iictxt
    iictxt = ictxt
    call psb_errcomm(iictxt,err)
  end subroutine psb_errcomm_ipk

  subroutine psb_perror_ipk(ictxt)
    integer(psb_ipk_), intent(in)   :: ictxt
    integer(psb_mpik_) :: iictxt
    iictxt = ictxt
    call psb_perror(iictxt)
  end subroutine psb_perror_ipk
#endif
  ! saves action to support error traceback
  ! also changes error action to "return"
  subroutine psb_erractionsave(err_act)
    integer(psb_ipk_), intent(out) :: err_act
    err_act    = err_action
    err_action = psb_act_ret_
  end subroutine psb_erractionsave


  ! return the action to take upon error occurrence
  subroutine psb_get_erraction(err_act)
    integer(psb_ipk_), intent(out) :: err_act
    err_act=err_action
  end subroutine psb_get_erraction

  ! sets the action to take upon error occurrence
  subroutine psb_set_erraction(err_act)
    integer(psb_ipk_), intent(in) :: err_act
    err_action=err_act
  end subroutine psb_set_erraction

  ! sets the action to take upon error occurrence
  subroutine psb_set_erract_return()
    err_action  = psb_act_ret_
  end subroutine psb_set_erract_return
  subroutine psb_set_erract_print()
    err_action  = psb_act_print_
  end subroutine psb_set_erract_print
  subroutine psb_set_erract_abort()
    err_action  = psb_act_abort_
  end subroutine psb_set_erract_abort

  function psb_is_erract_return() result(res)
    logical :: res
    res = (err_action == psb_act_ret_)
  end function psb_is_erract_return
  function psb_is_erract_print() result(res)
    logical :: res
    res = (err_action == psb_act_print_)
  end function psb_is_erract_print
  function psb_is_erract_abort() result(res)
    logical :: res
    res = (err_action == psb_act_abort_)
  end function psb_is_erract_abort


  ! restores error action previously saved with psb_erractionsave
  subroutine psb_erractionrestore(err_act)
    integer(psb_ipk_), intent(in) :: err_act
    err_action = err_act
  end subroutine psb_erractionrestore


  function  psb_get_debug_level()
    integer(psb_ipk_) :: psb_get_debug_level
    psb_get_debug_level = debug_level
  end function psb_get_debug_level

  subroutine psb_set_debug_level(level)
    integer(psb_ipk_), intent(in) :: level
    if (level >= 0) then
      debug_level = level
    else
      debug_level = 0
    end if
  end subroutine psb_set_debug_level

  function  psb_get_debug_unit()
    integer(psb_ipk_) :: psb_get_debug_unit
    psb_get_debug_unit = debug_unit
  end function psb_get_debug_unit

  subroutine psb_set_debug_unit(unit)
    integer(psb_ipk_), intent(in) :: unit
    if ((unit >= 0).or.(unit == psb_err_unit)&
         & .or.(unit == psb_out_unit)) then
      debug_unit = unit
    else
      debug_unit = 0
    end if
  end subroutine psb_set_debug_unit

  function  psb_get_serial_debug_level()
    integer(psb_ipk_) :: psb_get_serial_debug_level
    psb_get_serial_debug_level = serial_debug_level
  end function psb_get_serial_debug_level

  subroutine psb_set_serial_debug_level(level)
    integer(psb_ipk_), intent(in) :: level
    if (level >= 0) then
      serial_debug_level = level
    else
      serial_debug_level = 0
    end if
  end subroutine psb_set_serial_debug_level



  ! sets verbosity of the error message
  subroutine psb_set_errverbosity(v)
    integer(psb_ipk_), intent(in) :: v
    verbosity_level=v
  end subroutine psb_set_errverbosity



  ! returns number of errors
  function psb_get_numerr()
    integer(psb_ipk_) :: psb_get_numerr
    psb_get_numerr = error_stack%n_elems
  end function psb_get_numerr


  ! returns verbosity of the error message
  function psb_get_errverbosity()
    integer(psb_ipk_) :: psb_get_errverbosity

    psb_get_errverbosity=verbosity_level
  end function psb_get_errverbosity



  ! checks the status of the error condition
  function psb_get_errstatus()
    integer(psb_ipk_) :: psb_get_errstatus
    psb_get_errstatus = error_status
  end function psb_get_errstatus

  subroutine  psb_set_errstatus(ircode)
    integer(psb_ipk_) :: ircode
    if ((psb_no_err_<=ircode).and.(ircode <= psb_err_fatal_))&
         &  error_status=ircode
  end subroutine psb_set_errstatus

  function psb_errstatus_fatal() result(res)
    logical :: res
    res = (error_status == psb_err_fatal_)
  end function psb_errstatus_fatal

  function psb_errstatus_warning() result(res)
    logical :: res
    res = (error_status == psb_err_warning_)
  end function psb_errstatus_warning

  function psb_errstatus_ok() result(res)
    logical :: res
    res = (error_status == psb_no_err_)
  end function psb_errstatus_ok

  ! pushes an error on the error stack
  subroutine psb_stackpush(err_c, r_name, i_err, a_err)

    integer(psb_ipk_), intent(in)              ::  err_c
    character(len=*), intent(in)     ::  r_name
    character(len=*), optional       ::  a_err
    integer(psb_ipk_), optional                ::  i_err(5)

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
    nullify(new_node)

  end subroutine psb_stackpush

  ! pushes an error on the error stack
  subroutine psb_errpush(err_c, r_name, i_err, a_err)

    integer(psb_ipk_), intent(in)              ::  err_c
    character(len=*), intent(in)     ::  r_name
    character(len=*), optional       ::  a_err
    integer(psb_ipk_), optional                ::  i_err(5)

    type(psb_errstack_node), pointer     ::  new_node

    call psb_set_errstatus(psb_err_fatal_)
    call psb_stackpush(err_c, r_name, i_err, a_err)

  end subroutine psb_errpush

  ! pushes a warning on the error stack
  subroutine psb_warning_push(err_c, r_name, i_err, a_err)

    integer(psb_ipk_), intent(in)              ::  err_c
    character(len=*), intent(in)     ::  r_name
    character(len=*), optional       ::  a_err
    integer(psb_ipk_), optional                ::  i_err(5)

    type(psb_errstack_node), pointer     ::  new_node


    if (.not.psb_errstatus_fatal())&
         &  call psb_set_errstatus( psb_err_warning_)
    call psb_stackpush(err_c, r_name, i_err, a_err)
  end subroutine psb_warning_push


  ! pops an error from the error stack
  subroutine psb_ach_errpop(achmsg)
    character(len=psb_max_errmsg_len_), allocatable, intent(out) :: achmsg(:)
    integer(psb_ipk_)        ::  err_c
    character(len=20)        ::  r_name
    character(len=40)        ::  a_e_d
    integer(psb_ipk_)        ::  i_e_d(5)

    type(psb_errstack_node), pointer     ::  old_node

    if (error_stack%n_elems > 0) then 
      err_c      =  error_stack%top%err_code
      r_name     =  error_stack%top%routine
      i_e_d      =  error_stack%top%i_err_data
      a_e_d      =  error_stack%top%a_err_data
      call psb_errmsg(achmsg,err_c, r_name, i_e_d, a_e_d)
      old_node   => error_stack%top
      error_stack%top  => old_node%next
      error_stack%n_elems = error_stack%n_elems - 1
      deallocate(old_node)
    end if
    if (error_stack%n_elems == 0) error_status=0
      

  end subroutine psb_ach_errpop

  ! pops an error from the error stack
  subroutine psb_errpop(err_c, r_name, i_e_d, a_e_d)

    integer(psb_ipk_), intent(out)        ::  err_c
    character(len=20), intent(out)        ::  r_name
    character(len=40), intent(out)        ::  a_e_d
    integer(psb_ipk_), intent(out)        ::  i_e_d(5)

    type(psb_errstack_node), pointer     ::  old_node

    if (error_stack%n_elems > 0) then 
      err_c      =  error_stack%top%err_code
      r_name     =  error_stack%top%routine
      i_e_d      =  error_stack%top%i_err_data
      a_e_d      =  error_stack%top%a_err_data
      
      old_node   => error_stack%top
      error_stack%top  => old_node%next
      error_stack%n_elems = error_stack%n_elems - 1
      deallocate(old_node)
    end if
    if (error_stack%n_elems == 0) error_status=psb_no_err_
      

  end subroutine psb_errpop

  ! Clean the error stack
  subroutine psb_clean_errstack()

    integer(psb_ipk_)        ::  err_c
    character(len=20)        ::  r_name
    character(len=40)        ::  a_e_d
    integer(psb_ipk_)        ::  i_e_d(5)

    
    do while (psb_get_numerr() > 0)
      call psb_errpop(err_c, r_name, i_e_d, a_e_d)
    end do
    
  end subroutine psb_clean_errstack

  ! prints the error msg associated to a specific error code
  subroutine psb_ach_errmsg(achmsg,err_c, r_name, i_e_d, a_e_d,me)

    character(len=psb_max_errmsg_len_), allocatable, intent(out) :: achmsg(:)
    integer(psb_ipk_), intent(in)    ::  err_c
    character(len=20), intent(in)    ::  r_name
    character(len=40), intent(in)    ::  a_e_d
    integer(psb_ipk_), intent(in)    ::  i_e_d(5)
    integer(psb_mpik_), optional     ::  me

    character(len=psb_max_errmsg_len_) :: tmpmsg
  
    
    
    if(present(me)) then
      write(tmpmsg,&
           & '("Process: ",i0,".  PSBLAS Error (",i0,") in subroutine: ",a20)')&
           & me,err_c,trim(r_name)
    else
      write(tmpmsg,'("PSBLAS Error (",i0,") in subroutine: ",a)')&
           & err_c,trim(r_name)
    end if

    
    select case (err_c)
    case(:psb_success_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("error on calling perror. err_c must be greater than 0")')
      
    case(psb_err_pivot_too_small_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("pivot too small: ",i0,1x,a)')i_e_d(1),trim(a_e_d)

    case(psb_err_invalid_ovr_num_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid number of ovr:",i0)')i_e_d(1)

    case(psb_err_invalid_input_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid input")')

    case(psb_err_iarg_neg_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," cannot be less than 0")')i_e_d(1)
      write(achmsg(3),'("current value is ",i0)')i_e_d(2)

    case(psb_err_iarg_pos_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," cannot be greater than 0")')i_e_d(1)
      write(achmsg(3),'("current value is ",i0)')i_e_d(2)

    case(psb_err_input_value_invalid_i_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write(achmsg(3),'("current value is ",i0)')i_e_d(2)

    case(psb_err_input_asize_invalid_i_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Size of input array argument n. ",i0," is invalid.")')i_e_d(1)
      write(achmsg(3),'("Current value is ",i0)')i_e_d(2)

    case(psb_err_input_asize_small_i_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Size of input array argument n. ",i0," is too small.")')i_e_d(1)
      write(achmsg(3),'("Current value is ",i0," Should be at least ",i0)') i_e_d(2),i_e_d(3)

    case(psb_err_iarg_invalid_i_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write(achmsg(3),'("current value is ",a)')a_e_d(2:2)

    case(psb_err_iarg_not_gtia_ii_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("input argument n. ",i0," must be equal or greater than input argument n. ",i0)') &
           & i_e_d(1), i_e_d(3)
      write(achmsg(3),'("current values are ",i0," < ",i0)')&
           & i_e_d(2),i_e_d(5)

    case(psb_err_iarg_not_gteia_ii_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," must be greater than or equal to ",i0)')&
           & i_e_d(1),i_e_d(2)
      write(achmsg(3),'("current value is ",i0," < ",i0)')&
           & i_e_d(3), i_e_d(2)

    case(psb_err_iarg_invalid_value_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," in entry # ",i0," has an invalid value")')&
           & i_e_d(1:2)
      write(achmsg(3),'("current value is ",a)')trim(a_e_d)

    case(psb_err_asb_nrc_error_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Impossible error in ASB: nrow>ncol,")')
      write(achmsg(3),'("Actual values are ",i0," > ",i0)')i_e_d(1:2)
      !        ... csr format error ...

    case(psb_err_iarg2_neg_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument ia2(1) is less than 0")')
      write(achmsg(3),'("current value is ",i0)')i_e_d(1)
      !        ... csr format error ...

    case(psb_err_ia2_not_increasing_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("indices in ia2 array are not in  increasing order")')

    case(psb_err_ia1_not_increasing_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("indices in ia1 array are not in increasing order")')
      !        ... csr format error ...

    case(psb_err_ia1_badindices_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("indices in ia1 array are not within problem dimension")')
      write(achmsg(3),'("problem dimension is ",i0)')i_e_d(1)

    case(psb_err_invalid_args_combination_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("invalid combination of input arguments")')

    case(psb_err_invalid_pid_arg_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid process identifier in input array argument n. ",i0,".")')&
           & i_e_d(1)
      write(achmsg(3),'("Current value is ",i0)')i_e_d(2)

    case(psb_err_iarg_n_mbgtian_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," must be greater than input argument n. ",i0)')&
           & i_e_d(1:2)
      write(achmsg(3),'("current values are ",i0," < ",i0)') i_e_d(3:4)

    case(psb_err_dupl_cd_vl)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("there are duplicated entries in vl (input to cdall)")')
      !        ... coo format error ...
      !        ... coo format error ...

    case(psb_err_duplicate_coo)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("there are duplicated elements in coo format")')
      write(achmsg(3),'("and you have chosen psb_dupl_err_ ")')

    case(psb_err_invalid_input_format_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid input format ",a3)')&
           & a_e_d(1:3)

    case(psb_err_unsupported_format_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Format ",a3," not yet supported here")')&
           &a_e_d(1:3)

    case(psb_err_format_unknown_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Format ",a3," is unknown")')&
           & a_e_d(1:3)

    case(psb_err_iarray_outside_bounds_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("indices in input array are not within problem dimension ",2(i0,2x))')&
           &i_e_d(1:2)

    case(psb_err_iarray_outside_process_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("indices in input array are not belonging to the calling process ",i0)')&
           & i_e_d(1)

    case(psb_err_forgot_geall_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("To call this routine you must first call psb_geall on the same matrix")')

    case(psb_err_forgot_spall_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("To call this routine you must first call psb_spall on the same matrix")')

    case(psb_err_wrong_ins_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("Something went wrong before this call to ",a,", probably in cdins/spins")')&
           & trim(r_name)

    case(psb_err_iarg_mbeeiarra_i_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Input argument n. ",i0," must be equal to entry n. ",i0," in array input argument n.",i0)') &
           & i_e_d(1),i_e_d(4),i_e_d(3)
      write(achmsg(3),'("Current values are ",i0," != ",i0)')i_e_d(2), i_e_d(5)

    case(psb_err_mpi_error_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("MPI error:",i0)')i_e_d(1)

    case(psb_err_parm_differs_among_procs_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("Parameter n. ",i0," must be equal on all processes. ",i0)')i_e_d(1)

    case(psb_err_entry_out_of_bounds_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("Entry n. ",i0," out of ",i0," should be between 1 and ",i0," but is ",i0)')&
           & i_e_d(1),i_e_d(3),i_e_d(4),i_e_d(2)

    case(psb_err_inconsistent_index_lists_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Index lists are inconsistent: some indices are orphans")')

    case(psb_err_partfunc_toomuchprocs_)
      allocate(achmsg(4)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("partition function passed as input argument n. ",i0," returns number of processes")')&
           &i_e_d(1)
      write(achmsg(3),&
           & '("greater than No of grid s processes on global point ",i0,". Actual number of grid s ")')&
           &i_e_d(4)
      write(achmsg(4),'("processes is ",i0,", number returned is ",i0)')i_e_d(2),i_e_d(3)

    case(psb_err_partfunc_toofewprocs_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("partition function passed as input argument n. ",i0," returns number of processes")')&
           &i_e_d(1)
      write(achmsg(3),&
           &'("less or equal to 0 on global point ",i0,". Number returned is ",i0)')&
           &i_e_d(3),i_e_d(2)

    case(psb_err_partfunc_wrong_pid_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("partition function passed as input argument n. ",i0," returns wrong processes identifier")')&
           & i_e_d(1)
      write(achmsg(3),&
           & '("on global point ",i0,". Current value returned is : ",i0)')&
           & i_e_d(3),i_e_d(2)

    case(psb_err_no_optional_arg_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("One of the optional arguments  ",a," must be present")')&
           & trim(a_e_d)

    case(psb_err_optional_arg_pair_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("The optional arguments  ",a," must be present together")')&
           & trim(a_e_d)

    case(psb_err_arg_m_required_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Argument M is required when argument PARTS is specified")')

    case(psb_err_missing_override_method_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Base class method ",a," called: the class for ",a," is missing an overriding implementation")')&
           &  trim(r_name), trim(a_e_d)

    case(psb_err_invalid_dynamic_type_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("input argument n. ",i0," has a dynamic type not allowed here.")')&
           & i_e_d(1)
    case(psb_err_rectangular_mat_unsupported_)
      write(achmsg(2),&
           &'("This routine does not support rectangular matrices: ",i0, " /= ",i0)') &
           & i_e_d(1), i_e_d(2)

    case(psb_err_invalid_mat_state_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid state for sparse matrix")')

    case(psb_err_invalid_cd_state_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid state for communication descriptor")')

    case(psb_err_invalid_a_and_cd_state_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid combined state for A and DESC_A")')

    case(psb_err_invalid_vect_state_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid state for vector")')

    case(psb_err_invalid_matrix_sizes_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid combination of matrix sizes")')

    case(1125:1999)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("computational error. code: ",i0)')err_c

    case(psb_err_context_error_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(0,'("Parallel context error. Number of processes=-1")')

    case(psb_err_initerror_neugh_procs_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Initialization error: not enough processes available in the parallel environment")')

    case(psb_err_invalid_matrix_input_state_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid input state for matrix.")')

    case(psb_err_input_no_regen_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Input state for matrix is not adequate for regeneration.")')

    case(2233:2999)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("resource error. code: ",i0)')err_c

    case(3000:3009)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("sparse matrix representation ",a3," not yet implemented")')&
           &a_e_d(1:3)

    case(psb_err_lld_case_not_implemented_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("Case lld not equal matrix_data[N_COL_] is not yet implemented.")')

    case(psb_err_transpose_unsupported_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("transpose option for sparse matrix representation ",a3," not implemented")')&
           & a_e_d(1:3)

    case(psb_err_transpose_c_unsupported_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case trans = C is not yet implemented.")') 

    case(psb_err_transpose_not_n_unsupported_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case trans /= N is not yet implemented.")') 

    case(psb_err_only_unit_diag_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Only unit diagonal so far for triangular matrices. ")') 

    case(3023)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Cases DESCRA(1:1)=S  DESCRA(1:1)=T not yet implemented. ")') 

    case(3024)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Cases DESCRA(1:1)=G not yet implemented. ")') 

    case(psb_err_ja_nix_ia_niy_unsupported_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case ja /= ix or ia/=iy is not yet implemented.")')

    case(psb_err_ix_n1_iy_n1_unsupported_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case ix /= 1 or iy /= 1 is not yet implemented.")')

    case(3050)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case ix /= iy is not yet implemented.")')

    case(3060)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Case ix /= 1 is not yet implemented.")')

    case(3070)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("This operation is only implemented with no overlap.")')

    case(3080)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Decompostion type ",i0," not yet supported.")')&
           & i_e_d(1)

    case(3090)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Insert matrix mode not yet implemented.")')

    case(3100)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Error on index. Element has not been inserted")')
      write(achmsg(3),&
           & '("local index is: ",i0," and global index is:",i0)')&
           & i_e_d(1:2)

    case(psb_err_input_matrix_unassembled_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("Before you call this routine, you must assembly sparse matrix")')

    case(3111)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Before you call this routine, you must initialize the preconditioner")')

    case(3112)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Before you call this routine, you must build the preconditioner")')

    case(3113:3998)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("miscellaneus error. code: ",i0)')err_c

    case(psb_err_missing_aux_lib_)
      allocate(achmsg(3)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           &'("This method requires an external support library.")')
      write(achmsg(3),&
           &'("Fix configure and rebuild the software.")')

    case(psb_err_alloc_dealloc_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Allocation/deallocation error")')

    case(psb_err_internal_error_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Internal error: ",a)') &
           & trim(a_e_d)

    case(psb_err_from_subroutine_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Error from call to subroutine ",a)')&
           & trim(a_e_d)

    case(psb_err_from_subroutine_non_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Error from call to a subroutine ")')

    case(psb_err_from_subroutine_i_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Error ",i0," from call to a subroutine ")')&
           & i_e_d(1)

    case(psb_err_from_subroutine_ai_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Error from call to subroutine ",a," ",i0)')&
           & trim(a_e_d),i_e_d(1)

    case(psb_err_alloc_request_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Error on allocation request for ",i0," items of type ",a)')&
           & i_e_d(1),trim(a_e_d)

    case(4110)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),&
           & '("Error ",i0," from call to an external package in subroutine ",a)')&
           &i_e_d(1),trim(a_e_d)

    case(psb_err_invalid_istop_)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid ISTOP: ",i0)')i_e_d(1)

    case(5002)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid PREC: ",i0)')i_e_d(1)

    case(5003)
      allocate(achmsg(2)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("Invalid PREC: ",a3)')a_e_d(1:3)

    case default
      allocate(achmsg(4)) 
      achmsg(1) = tmpmsg
      write(achmsg(2),'("unknown error (",i0,") in subroutine ",a)')&
           & err_c,trim(r_name)
      write(achmsg(3),'(5(i0,2x))') i_e_d
      write(achmsg(4),'(a)') trim(a_e_d)

    end select

  end subroutine psb_ach_errmsg


  ! prints the error msg associated to a specific error code
  subroutine psb_errmsg(iunit, err_c, r_name, i_e_d, a_e_d,me)
    integer(psb_ipk_), intent(in)  :: iunit
    integer(psb_ipk_), intent(in)  ::  err_c
    character(len=20), intent(in)  ::  r_name
    character(len=40), intent(in)  ::  a_e_d
    integer(psb_ipk_), intent(in)  ::  i_e_d(5)
    integer(psb_mpik_), optional   ::  me

    integer(psb_ipk_) :: i
    character(len=psb_max_errmsg_len_), allocatable :: achmsg(:)

    call psb_ach_errmsg(achmsg,err_c, r_name, i_e_d, a_e_d,me)
    
    do i=1,size(achmsg)
      write(iunit,'(a)') trim(achmsg(i))
    end do

  end subroutine psb_errmsg

end module psb_error_mod
