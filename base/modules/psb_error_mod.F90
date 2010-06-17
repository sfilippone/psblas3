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
  use psb_const_mod
  integer, parameter, public :: psb_act_ret_=0, psb_act_abort_=1, psb_no_err_=0
  integer, parameter, public :: psb_debug_ext_=1, psb_debug_outer_=2
  integer, parameter, public :: psb_debug_comp_=3, psb_debug_inner_=4
  integer, parameter, public :: psb_debug_serial_=8, psb_debug_serial_comp_=9
 
  !
  !     Error handling 
  !
  public psb_errpush, psb_error, psb_get_errstatus,&
       & psb_errpop, psb_errmsg, psb_errcomm, psb_get_numerr, &
       & psb_get_errverbosity, psb_set_errverbosity, &
       & psb_erractionsave, psb_erractionrestore, &
       & psb_get_erraction, psb_set_erraction, &
       & psb_get_debug_level, psb_set_debug_level,&
       & psb_get_debug_unit, psb_set_debug_unit,&
       & psb_get_serial_debug_level, psb_set_serial_debug_level

 
  interface psb_error
    subroutine psb_serror()
    end subroutine psb_serror
    subroutine psb_perror(ictxt)
      integer, intent(in)     ::  ictxt
    end subroutine psb_perror
  end interface

  interface 
    subroutine psb_errcomm(ictxt, err)
      integer, intent(in)   :: ictxt
      integer, intent(inout):: err
    end subroutine psb_errcomm
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
  integer, save            :: error_unit=0

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



  ! sets verbosity of the error message
  subroutine psb_set_errverbosity(v)
    integer, intent(in) :: v
    verbosity_level=v
  end subroutine psb_set_errverbosity



  ! returns number of errors
  function psb_get_numerr()
    integer :: psb_get_numerr

    psb_get_numerr = error_stack%n_elems
  end function psb_get_numerr


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



  ! prints the error msg associated to a specific error code
  subroutine psb_errmsg(err_c, r_name, i_e_d, a_e_d,me)

    integer, intent(in)              ::  err_c
    character(len=20), intent(in)    ::  r_name
    character(len=40), intent(in)    ::  a_e_d
    integer, intent(in)              ::  i_e_d(5)
    integer, optional                ::  me

    if(present(me)) then
      write(error_unit,'("Process: ",i0,".  PSBLAS Error (",i0,") in subroutine: ",a20)')&
           & me,err_c,trim(r_name)
    else
      write(error_unit,'("PSBLAS Error (",i0,") in subroutine: ",a)')err_c,trim(r_name)
    end if


    select case (err_c)
    case(:psb_success_)
      write (error_unit,'("error on calling sperror. err_c must be greater than 0")')
    case(psb_err_pivot_too_small_)
      write (error_unit,'("pivot too small: ",i0,1x,a)')i_e_d(1),trim(a_e_d)
    case(psb_err_invalid_ovr_num_)
      write (error_unit,'("Invalid number of ovr:",i0)')i_e_d(1)
    case(psb_err_invalid_input_)
      write (error_unit,'("Invalid input")')

    case(psb_err_iarg_neg_)
      write (error_unit,'("input argument n. ",i0," cannot be less than 0")')i_e_d(1)
      write (error_unit,'("current value is ",i0)')i_e_d(2)

    case(psb_err_iarg_pos_)
      write (error_unit,'("input argument n. ",i0," cannot be greater than 0")')i_e_d(1)
      write (error_unit,'("current value is ",i0)')i_e_d(2)
    case(psb_err_input_value_invalid_i_)
      write (error_unit,'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write (error_unit,'("current value is ",i0)')i_e_d(2)
    case(psb_err_input_asize_invalid_i_)
      write (error_unit,'("Size of input array argument n. ",i0," is invalid.")')i_e_d(1)
      write (error_unit,'("Current value is ",i0)')i_e_d(2)
    case(psb_err_iarg_invalid_i_)
      write (error_unit,'("input argument n. ",i0," has an invalid value")')i_e_d(1)
      write (error_unit,'("current value is ",a)')a_e_d(2:2)
    case(psb_err_iarg_not_gtia_ii_)
      write (error_unit,'("input argument n. ",i0," must be equal or greater than input argument n. ",i0)') &
           & i_e_d(1), i_e_d(3)
      write (error_unit,'("current values are ",i0," < ",i0)')&
           & i_e_d(2),i_e_d(5)
    case(psb_err_iarg_not_gteia_ii_)
      write (error_unit,'("input argument n. ",i0," must be greater than or equal to ",i0)')&
           & i_e_d(1),i_e_d(2)
      write (error_unit,'("current value is ",i0," < ",i0)')&
           & i_e_d(3), i_e_d(2)
    case(psb_err_iarg_invalid_value_)
      write (error_unit,'("input argument n. ",i0," in entry # ",i0," has an invalid value")')&
           & i_e_d(1:2)
      write (error_unit,'("current value is ",a)')trim(a_e_d)
    case(psb_err_asb_nrc_error_)
      write (error_unit,'("Impossible error in ASB: nrow>ncol,")')
      write (error_unit,'("Actual values are ",i0," > ",i0)')i_e_d(1:2)
      !        ... csr format error ...
    case(psb_err_iarg2_neg_)
      write (error_unit,'("input argument ia2(1) is less than 0")')
      write (error_unit,'("current value is ",i0)')i_e_d(1)
      !        ... csr format error ...
    case(psb_err_ia2_not_increasing_)
      write (error_unit,'("indices in ia2 array are not in  increasing order")')
    case(psb_err_ia1_not_increasing_)
      write (error_unit,'("indices in ia1 array are not in increasing order")')
      !        ... csr format error ...
    case(psb_err_ia1_badindices_)
      write (error_unit,'("indices in ia1 array are not within problem dimension")')
      write (error_unit,'("problem dimension is ",i0)')i_e_d(1)
    case(psb_err_invalid_args_combination_)
      write (error_unit,'("invalid combination of input arguments")')
    case(psb_err_invalid_pid_arg_)
      write (error_unit,'("Invalid process identifier in input array argument n. ",i0,".")')&
           & i_e_d(1)
      write (error_unit,'("Current value is ",i0)')i_e_d(2)
    case(psb_err_iarg_n_mbgtian_)
      write (error_unit,'("input argument n. ",i0," must be greater than input argument n. ",i0)')&
           & i_e_d(1:2)
      write (error_unit,'("current values are ",i0," < ",i0)') i_e_d(3:4)
    case(psb_err_dupl_cd_vl)
      write (error_unit,'("there are duplicated entries in vl (input to cdall)")')
      !        ... coo format error ...
      !        ... coo format error ...
    case(psb_err_duplicate_coo)
      write (error_unit,'("there are duplicated elements in coo format")')
      write (error_unit,'("and you have chosen psb_dupl_err_ ")')
    case(psb_err_invalid_input_format_)
      write (error_unit,'("Invalid input format ",a3)')&
           & a_e_d(1:3)
    case(psb_err_unsupported_format_)
      write (error_unit,'("Format ",a3," not yet supported here")')&
           &a_e_d(1:3)
    case(psb_err_format_unknown_)
      write (error_unit,'("Format ",a3," is unknown")')&
           & a_e_d(1:3)
    case(psb_err_iarray_outside_bounds_)
      write (error_unit,'("indices in input array are not within problem dimension ",2(i0,2x))')&
           &i_e_d(1:2)
    case(psb_err_iarray_outside_process_)
      write (error_unit,'("indices in input array are not belonging to the calling process ",i0)')&
           & i_e_d(1)
    case(psb_err_forgot_geall_)
      write (error_unit,'("To call this routine you must first call psb_geall on the same matrix")')
    case(psb_err_forgot_spall_)
      write (error_unit,'("To call this routine you must first call psb_spall on the same matrix")')
    case(psb_err_wrong_ins_)
      write (0,'("Something went wrong before this call to ",a,", probably in cdins/spins")')&
           & trim(r_name)
    case(psb_err_iarg_mbeeiarra_i_)
      write (error_unit,&
           & '("Input argument n. ",i0," must be equal to entry n. ",i0," in array input argument n.",i0)') &
           & i_e_d(1),i_e_d(4),i_e_d(3)
      write (error_unit,'("Current values are ",i0," != ",i0)')i_e_d(2), i_e_d(5)
    case(psb_err_mpi_error_)
      write (error_unit,'("MPI error:",i0)')i_e_d(1)
    case(psb_err_parm_differs_among_procs_)
      write (error_unit,'("Parameter n. ",i0," must be equal on all processes. ",i0)')i_e_d(1)
    case(psb_err_entry_out_of_bounds_)
      write (error_unit,'("Entry n. ",i0," out of ",i0," should be between 1 and ",i0," but is ",i0)')&
           & i_e_d(1),i_e_d(3),i_e_d(4),i_e_d(2)
    case(psb_err_inconsistent_index_lists_)
      write (error_unit,'("Index lists are inconsistent: some indices are orphans")')
    case(psb_err_partfunc_toomuchprocs_)
      write (error_unit,&
           &'("partition function passed as input argument n. ",i0," returns number of processes")')&
           &i_e_d(1)
      write (error_unit,&
           & '("greater than No of grid s processes on global point ",i0,". Actual number of grid s ")')&
           &i_e_d(4)
      write (error_unit,'("processes is ",i0,", number returned is ",i0)')i_e_d(2),i_e_d(3)
    case(psb_err_partfunc_toofewprocs_)
      write (error_unit,'("partition function passed as input argument n. ",i0," returns number of processes")')&
           &i_e_d(1)
      write (error_unit,'("less or equal to 0 on global point ",i0,". Number returned is ",i0)')&
           &i_e_d(3),i_e_d(2)
    case(psb_err_partfunc_wrong_pid_)
      write (error_unit,&
           &'("partition function passed as input argument n. ",i0," returns wrong processes identifier")')&
           & i_e_d(1)
      write (error_unit,'("on global point ",i0,". Current value returned is : ",i0)')&
           & i_e_d(3),i_e_d(2)
    case(psb_err_no_optional_arg_)
      write (error_unit,'("One of the optional arguments  ",a," must be present")')&
           & trim(a_e_d)
    case(psb_err_arg_m_required_)
      write (error_unit,'("Argument M is required when argument PARTS is specified")')
    case(psb_err_spmat_invalid_state_)
      write (error_unit,&
           & '("Sparse Matrix and descriptors are in an invalid state for this subroutine call: ",i0)')&
           &i_e_d(1)
    case(psb_err_missing_override_method_)
      write (error_unit,&
           & '("Base class method ",a," called: the class for ",a," is missing an overriding implementation")')&
           &  trim(r_name), trim(a_e_d)
    case (psb_err_invalid_cd_state_)
      write (error_unit,'("Invalid state for communication descriptor")')
    case (psb_err_invalid_a_and_cd_state_)
      write (error_unit,'("Invalid combined state for A and DESC_A")')
    case(1124:1999)
      write (error_unit,'("computational error. code: ",i0)')err_c
    case(psb_err_context_error_)
      write (0,'("Parallel context error. Number of processes=-1")')
    case(psb_err_initerror_neugh_procs_)
      write (error_unit,&
           & '("Initialization error: not enough processes available in the parallel environment")')
    case(psb_err_invalid_matrix_input_state_)
      write (error_unit,'("Invalid input state for matrix.")')
    case(psb_err_input_no_regen_)
      write (error_unit,'("Input state for matrix is not adequate for regeneration.")')
    case (2233:2999)
      write(error_unit,'("resource error. code: ",i0)')err_c
    case(3000:3009)
      write (error_unit,&
           & '("sparse matrix representation ",a3," not yet implemented")')&
           &a_e_d(1:3)
    case(psb_err_lld_case_not_implemented_)
      write (error_unit,&
           &'("Case lld not equal matrix_data[N_COL_] is not yet implemented.")')
    case(psb_err_transpose_unsupported_)
      write (error_unit,&
           & '("transpose option for sparse matrix representation ",a3," not implemented")')&
           & a_e_d(1:3)
    case(psb_err_transpose_c_unsupported_)
      write (error_unit,'("Case trans = C is not yet implemented.")') 
    case(psb_err_transpose_not_n_unsupported_)
      write (error_unit,'("Case trans /= N is not yet implemented.")') 
    case(psb_err_only_unit_diag_)
      write (error_unit,'("Only unit diagonal so far for triangular matrices. ")') 
    case(3023)
      write (error_unit,'("Cases DESCRA(1:1)=S  DESCRA(1:1)=T not yet implemented. ")') 
    case(3024)
      write (error_unit,'("Cases DESCRA(1:1)=G not yet implemented. ")') 
    case(psb_err_ja_nix_ia_niy_unsupported_)
      write (error_unit,'("Case ja /= ix or ia/=iy is not yet implemented.")')
    case(psb_err_ix_n1_iy_n1_unsupported_)
      write (error_unit,'("Case ix /= 1 or iy /= 1 is not yet implemented.")')
    case(3050)
      write (error_unit,'("Case ix /= iy is not yet implemented.")')
    case(3060)
      write (error_unit,'("Case ix /= 1 is not yet implemented.")')
    case(3070)
      write (error_unit,&
           & '("This operation is only implemented with no overlap.")')
    case(3080)
      write (error_unit,&
           & '("Decompostion type ",i0," not yet supported.")')&
           & i_e_d(1)
    case(3090)
      write (error_unit,'("Insert matrix mode not yet implemented.")')
    case(3100)
      write (error_unit,&
           & '("Error on index. Element has not been inserted")')
      write (error_unit,&
           & '("local index is: ",i0," and global index is:",i0)')&
           & i_e_d(1:2)
    case(psb_err_input_matrix_unassembled_)
      write (error_unit,&
           &'("Before you call this routine, you must assembly sparse matrix")')
    case(3111)
      write (error_unit,&
           & '("Before you call this routine, you must initialize the preconditioner")')
    case(3112)
      write (error_unit,&
           & '("Before you call this routine, you must build the preconditioner")')
    case(3113:3999)
      write(error_unit,'("miscellaneus error. code: ",i0)')err_c
    case(psb_err_alloc_dealloc_)
      write(error_unit,'("Allocation/deallocation error")')
    case(psb_err_internal_error_)
      write(error_unit,'("Internal error: ",a)') &
           & trim(a_e_d)
    case(psb_err_from_subroutine_)
      write (error_unit,'("Error from call to subroutine ",a)')&
           & trim(a_e_d)
    case(psb_err_from_subroutine_non_)
      write (error_unit,'("Error from call to a subroutine ")')
    case(psb_err_from_subroutine_i_)
      write (error_unit,'("Error ",i0," from call to a subroutine ")')&
           & i_e_d(1)
    case(psb_err_from_subroutine_ai_)
      write (error_unit,'("Error from call to subroutine ",a," ",i0)')&
           & trim(a_e_d),i_e_d(1)
    case(psb_err_alloc_request_)
      write (error_unit,&
           & '("Error on allocation request for ",i0," items of type ",a)')&
           & i_e_d(1),trim(a_e_d)
    case(4110)
      write (error_unit,&
           & '("Error ",i0," from call to an external package in subroutine ",a)')&
           &i_e_d(1),trim(a_e_d)
    case (psb_err_invalid_istop_)
      write (error_unit,'("Invalid ISTOP: ",i0)')i_e_d(1)
    case (5002)
      write (error_unit,'("Invalid PREC: ",i0)')i_e_d(1)
    case (5003)
      write (error_unit,'("Invalid PREC: ",a3)')a_e_d(1:3)
    case default
      write(error_unit,'("unknown error (",i0,") in subroutine ",a)')&
           & err_c,trim(r_name)
      write(error_unit,'(5(i0,2x))') i_e_d
      write(error_unit,'(a)') trim(a_e_d)

    end select

  end subroutine psb_errmsg



end module psb_error_mod
