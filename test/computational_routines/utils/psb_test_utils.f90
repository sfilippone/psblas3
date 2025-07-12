!> @file psb_test_utils.f90ù
!!
!! @brief PSBLAS Test Utilities Module:
!!        This module provides utility functions for testing PSBLAS computational routines.
!!        It includes functions for initializing test environments, generating input vectors,
!!        checking results, and logging test outcomes. It also implements validation criteria
!!        for single and double precision computations.
!!
!! @date 2023-10-01
!! @author Luca Pepé Sciarria, Simone Staccone
!! @university Tor Vergata, Rome, Italy
!! @version 1.0
!!
module psb_test_utils
    use psb_base_mod
    use psb_util_mod
    use, intrinsic :: ieee_arithmetic

    implicit none

    ! Define the enumeration values to represent testing criteria
    integer, parameter :: VALUE = 1
    integer, parameter :: GAMMA = 2

    ! Define test metadata struct
    type :: test_info_
        integer(psb_ipk_)   :: current_test = 1         !> The test that is currently beeing run  
        integer(psb_ipk_)   :: total_tests = 1          !> The number of the total tests to run
        integer(psb_ipk_)   :: success = 0              !> the number of tests that succeded
        integer(psb_ipk_)   :: failure = 0              !> The number of tests that failed
        integer(psb_ipk_)   :: output_unit = 6          !> The output file handles (stdout by default)
        character(len=32)   :: kernel_name = "unknown"  !> The PSBLAS kernel that is beeing tested
        integer(psb_ipk_)   :: threshold_type = VALUE   !> The criteria used to pass a test (VAL,...)        
        real(psb_dpk_)      :: threshold = 0.0          !> The threashold value used for acceptance
        type(psb_ctxt_type) :: ctxt                     !> The PSBLAS context (Used for MPI communications)
        integer(psb_ipk_)   :: my_rank = 0              !> The rank of the current process in the MPI communicator
    end type test_info_


contains

    include 'psb_test_log.inc'

    !> @brief Function to initialize the test environment.
    !!        It is used to set the output unit and the kernel name.
    !!
    !! @param test_info The test information structure to be initialized.
    !!
    subroutine psb_test_init(test_info)
        type(test_info_), intent(inout) :: test_info
        integer(psb_ipk_)               :: output_unit, info
        ! MPI variables
        integer(psb_ipk_)               :: my_rank, np


        call psb_init(test_info%ctxt)
        call psb_info(test_info%ctxt,test_info%my_rank,np)

        ! Check if the kernel name is set
        if (trim(test_info%kernel_name) == "unknown") then
            write(*, '(A)') "Error: Kernel name is not set. Please set the kernel name before running the test."
            call psb_exit(test_info%ctxt)
        end if

        ! Set the output unit based on the number of processes


        ! Set the output unit to stdout by default
        if(np == 1) then
            open(newunit=output_unit, file=trim(test_info%kernel_name)//'_test.log', &
            & status='replace', action='write', iostat=info)
        else
            open(newunit=output_unit, file=trim(test_info%kernel_name)//'_test.log', &
            & status='old', action='write', position='append', iostat=info)
        end if

        ! Check if the file was opened successfully
        if (info /= 0) then
           write(*, '(A,I0)') "Error opening log file for kernel ", test_info%kernel_name
           write(*, '(A,I0)') "I/O Status Code:", info
           write(*, '(A)') "Please check if the file is accessible and writable."
           call psb_exit(test_info%ctxt)
        end if

        test_info%output_unit = output_unit

        write(test_info%output_unit,'(A,A)') 'Welcome to PSBLAS version: ',psb_version_string_
        write(test_info%output_unit,'(A)') 'This is the psb_gedot_test sample program'
        write(test_info%output_unit,'(A,I0)') 'Number of processes used in this computation: ', np
        write(test_info%output_unit,'(A)') ''

    end subroutine

    !> @brief Function to finalize the test environment, it is used to close the output unit 
    !!        and finalize the PSBLAS context.
    !!
    !! @param test_info The test information structure to be finalized.
    !!
    subroutine psb_test_exit(test_info)
        type(test_info_), intent(inout) :: test_info
        integer(psb_ipk_)               :: info

        ! Finalize test
        if(test_info%my_rank == psb_root_) then
            write(*,'(A)') "[INFO]    Tests completed successfully!"
            write(*,'(A,I0)') "          Passed: ", test_info%success
            write(*,'(A,I0)') "          Failed: ", test_info%failure
            write(*,'(A,I0)') "          Total:  ", test_info%total_tests
            write(*,'(A)') "[INFO]    Check " // trim(test_info%kernel_name) // "_test.log for a full description"
            write(test_info%output_unit, *) ""

            ! Close the output unit
            if (test_info%output_unit /= 0) then
                close(test_info%output_unit, iostat=info)
                if (info /= 0) then
                    write(*, '(A,I0,A,I0)') "[ERROR]    Error closing log file for kernel ", &
                    & test_info%kernel_name, " I/O Status Code ", info
                    write(*,'(A)') ''
                end if
            end if
        end if

        ! Finalize PSBLAS context
        call psb_exit(test_info%ctxt)
    end subroutine


    !> @brief Function to validate the test information structure.
    !!        It sets the threshold value based on the threshold type.
    !!
    !! @param result_single The single precision result to be validated.
    !! @param result_double The double precision result to be validated.
    !! @param test_info The test information structure to be validated.
    !! @param arr_size The size of the array to be used for validation.
    !!
    function psb_test_validate(result_single, result_double, test_info, arr_size) result(pass)
        type(test_info_)    :: test_info
        real(psb_spk_)      :: result_single
        real(psb_dpk_)      :: result_double
        integer(psb_ipk_)   :: arr_size,int_digits
        real(psb_dpk_)      :: gamma_n, unit_roundoff, delta, rel_err
        logical             :: pass   

        unit_roundoff = 5.96D-08 !! 1.11D-16

        call shift_decimal_single(result_single,int_digits)
        call shift_decimal_double(result_double,int_digits)

        delta = abs(result_double - real(result_single,psb_dpk_))


        if(test_info%threshold_type == VALUE) then
            if(delta < test_info%threshold) then
                pass = .true.
            else
                pass = .false.
            end if
        else if(test_info%threshold_type == GAMMA) then
            gamma_n = (arr_size * unit_roundoff) / (1 - unit_roundoff * arr_size)
            test_info%threshold = gamma_n
            rel_err = delta / real(result_single,psb_dpk_)

            !! Handle case when result_single is zero
            if ((ieee_is_nan(rel_err)).and.(result_single == 0.0)) then
                rel_err = 0.0
            end if
            
            if( rel_err < gamma_n) then
                pass = .true.
            else
                pass = .false.
            end if
        else
            write(test_info%output_unit,'(A,I0)') "Error: Invalid threshold type: ", test_info%threshold_type
            write(*,*)
            call psb_test_exit(test_info)
        end if

    end function psb_test_validate

    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    !! @param arr_size The size of the vectors to be generated.
    !!
    subroutine psb_test_generate_input_vectors(arr_size)
        integer(psb_ipk_)                           :: arr_size
        real(psb_dpk_), allocatable                 :: x(:), y(:)
        integer(psb_ipk_)                           :: i, info
        logical                                     :: exists


        ! Check if output directory exists
        inquire(file='vectors/', exist=exists)
        if (.not.exists) then
            call system('mkdir vectors/')
        end if

        allocate(x(arr_size))
        allocate(y(arr_size))
        
        call random_init(repeatable=.true.,image_distinct=.true.) 
        call random_number(x)
        call random_number(y)

        ! Write only positive in x_1
        call mm_array_write(x,"Positive vector",info,filename="vectors/x1.mtx")
        call mm_array_write(y,"Positive vector",info,filename="vectors/y1.mtx")

        ! Write only negative in x_2
        do i=1,arr_size 
            x(i) = -x(i)
        end do  

        do i=1,arr_size
            y(i) = -y(i)
        end do

        call mm_array_write(x,"Negative vector",info,filename="vectors/x2.mtx")
        call mm_array_write(y,"Negative vector",info,filename="vectors/y2.mtx")


        ! Since numbers are less than one and always positive, we have to generate negative ones subtractiong 50
        do i=1,arr_size
            x(i) = -x(i) ! Make the values positive again  
            x(i) = x(i) - 0.5
        end do   

        do i=1,arr_size
            y(i) = -y(i) ! Make the values positive again  
            y(i) = y(i) - 0.5
        end do

        ! Write random in x_3
        call mm_array_write(x,"Random vector",info,filename="vectors/x3.mtx")
        call mm_array_write(y,"Random vector",info,filename="vectors/y3.mtx")

        ! Write zero in x_4
        do i=1,arr_size
            x(i) = 0
        end do   

        do i=1,arr_size
            y(i) = 0
        end do

        call mm_array_write(x,"Null vector",info,filename="vectors/x4.mtx")
        call mm_array_write(y,"Null vector",info,filename="vectors/y4.mtx")

        deallocate(x)
        deallocate(y)

    end subroutine

    !> @brief Subroutine to shift the decimal point of a single precision number
    !!        and count the number of digits in the integer part.
    !!
    !! @param num The single precision number whose decimal point is to be shifted.
    !! @param int_digits The integer to store the number of digits in the integer part.
    !!
    subroutine shift_decimal_double(num, int_digits)
        real(psb_dpk_),intent(inout)    :: num
        integer(psb_ipk_), intent(out)  :: int_digits
        integer(psb_ipk_)               :: n_digits
        character(len=20)               :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(num))
      
        ! Count number of digits
        int_digits = int(floor(log10(abs(num)))) + 1
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        num = abs(num) / 10.0**n_digits
      
    end subroutine

    !> @brief Subroutine to shift the decimal point of a single precision number
    !!        and count the number of digits in the integer part.
    !!
    !! @param num The single precision number whose decimal point is to be shifted.
    !! @param int_digits The integer to store the number of digits in the integer part.
    !!
    subroutine shift_decimal_single(num, int_digits)
        real(psb_spk_),intent(inout)    :: num
        integer(psb_ipk_), intent(out)  :: int_digits
        integer(psb_ipk_)               :: n_digits
        character(len=20)               :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(num))
      
        ! Count number of digits
        int_digits = int(floor(log10(abs(num)))) + 1
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        num = abs(num) / 10.0**n_digits
      
    end subroutine

    !> @brief Subroutine to check the results of a single and double precision computation.
    !!        It compares the results and logs the outcome.
    !!
    !! @param result_single The single precision result to be checked.
    !! @param result_double The double precision result to be checked.
    !! @param test_info The test information structure containing the threshold and logging details.
    !! @param arr_size The size of the array used in the computation.
    !!
    subroutine psb_test_single_double_check(result_single, result_double, test_info, arr_size)
        integer(psb_ipk_)               :: int_digits, arr_size 
        real(psb_spk_), intent(inout)   :: result_single
        real(psb_dpk_), intent(inout)   :: result_double
        real(psb_dpk_)                  :: delta
        type(test_info_), intent(inout) :: test_info


        call psb_test_progress_bar(test_info)

        if(psb_test_validate(result_single, result_double, test_info, arr_size)) then 
            call psb_test_log_passed(test_info)
            test_info%success = test_info%success + 1  
        else
            call psb_test_log_failed(test_info)
            write(psb_out_unit,'(A,F20.10)') "Single precision result:  ", result_single
            write(psb_out_unit,'(A,F20.10)') "Double precision result:  ", result_double
            write(psb_out_unit,'(A,F20.10)') "Computed delta:           ", delta
            write(psb_out_unit,'(A,F20.10)') "Threshold used:           ", test_info%threshold
            test_info%failure = test_info%failure + 1
        end if
        
    end subroutine
!! 
    !! subroutine psb_parallel_check()
    !! 
    !! end subroutine
!! 
    !! subroutine psb_threshold_check()
    !! 
    !! end subroutine

end module psb_test_utils