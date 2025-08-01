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
    integer, parameter :: DEFAULT = -1  
    integer, parameter :: VALUE = 1
    integer, parameter :: GAMMA = 2

    ! Define test metadata struct
    type :: psb_test_info
        integer(psb_ipk_)   :: current_test = 1         !> The test that is currently beeing run  
        integer(psb_ipk_)   :: total_tests = 1          !> The number of the total tests to run
        integer(psb_ipk_)   :: success = 0              !> the number of tests that succeded
        integer(psb_ipk_)   :: failure = 0              !> The number of tests that failed
        integer(psb_ipk_)   :: output_unit = 6          !> The output file handles (stdout by default)
        character(len=32)   :: kernel_name = "default"  !> The PSBLAS kernel that is beeing tested
        integer(psb_ipk_)   :: threshold_type = DEFAULT !> The criteria used to pass a test (VAL,...)        
        real(psb_dpk_)      :: threshold = 0.0          !> The threashold value used for acceptance
        type(psb_ctxt_type) :: ctxt                     !> The PSBLAS context (Used for MPI communications)
        integer(psb_ipk_)   :: my_rank = 0              !> The rank of the current process in the MPI communicator
        integer(psb_ipk_)   :: np = 1                   !> The number of processes in the MPI communicator
        integer(psb_ipk_)   :: mat_dist = DEFAULT       !> The distribution of the matrix (default is block row distribution)
    end type psb_test_info


contains

    include 'psb_test_log.inc'
    include "psb_test_env.inc"
 

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

    !> @brief Subroutine to save the result of a single precision computation
    !!        to a file in the results directory.
    !!
    !! @param result_single The single precision result to be saved.
    !! @param test_info The test information structure containing the current test count.
    !!
    subroutine psb_test_save_result(result_single, test_info)
        type(psb_test_info), intent(inout)  :: test_info
        real(psb_spk_), intent(in)          :: result_single
        integer(psb_ipk_)                   :: info, unit
        character(len=32)                   :: filename
        logical                             :: exists

        ! Check if results directory exists
        inquire(file='results/', exist=exists)
        if (.not.exists) then
            call system('mkdir results/')
        end if

        ! Set the filename based on the test count
        write(filename, '(A,I0,A)') 'results/result_', test_info%current_test, '.txt'

        ! Open the file for writing
        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=info)

        ! Check if the file was opened successfully
        if (info /= 0) then
            write(*, '(A,I0)') "Error opening result file: ", info
            return
        end if

        ! Write the result to the file
        write(unit, '(F20.10)') result_single

        ! Close the file
        close(unit, iostat=info)
    end subroutine

    !> @brief Subroutine to save the result of a single precision vector computation
    !!        to a file in the results directory.
    !! @param result_single The single precision vector result to be saved.
    !! @param test_info The test information structure containing the current test count.
    !! 
    subroutine psb_test_save_vector_result(result_single, test_info)
        type(psb_test_info), intent(inout)      :: test_info
        real(psb_spk_), allocatable,intent(in)  :: result_single(:)
        integer(psb_ipk_)                       :: info, unit
        character(len=32)                       :: filename
        logical                                 :: exists

        ! Check if results directory exists
        inquire(file='results/', exist=exists)
        if (.not.exists) then
            call system('mkdir results/')
        end if

        ! Set the filename based on the test count
        write(filename, '(A,I0,A)') 'results/result_', test_info%current_test, '.txt'

        ! Open the file for writing
        open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=info)

        ! Check if the file was opened successfully
        if (info /= 0) then
            write(*, '(A,I0)') "Error opening result file: ", info
            return
        end if

        ! Close the file
        close(unit, iostat=info)

        ! Write the result to the file
        call mm_array_write(result_single,"",info,filename=filename)
        if (info /= 0) then
            write(*, '(A,I0)') "Error writing result file: ", info
            return
        end if


    end subroutine


    !> @brief Subroutine to shift the decimal point of a single precision number
    !!        and count the number of digits in the integer part.
    !!
    !! @param num The single precision number whose decimal point is to be shifted.
    !!
    !! @return shifted_num The single precision number with the decimal point shifted.
    !!
    function shift_decimal_double(num) result(shifted_num)
        real(psb_dpk_)      :: num, shifted_num
        integer(psb_ipk_)   :: n_digits
        character(len=20)   :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(num))
      
        ! Count number of digits
        ! int_digits = int(floor(log10(abs(num)))) + 1
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        shifted_num = abs(num) / 10.0**n_digits
      
    end function shift_decimal_double

    !> @brief Function to shift the decimal point of a single precision number.
    !!
    !! @param num The single precision number whose decimal point is to be shifted.
    !!
    !! @return shifted_num The single precision number with the decimal point shifted.
    !!
    function shift_decimal_single(num) result(shifted_num)
        real(psb_spk_)      :: num, shifted_num
        integer(psb_ipk_)   :: n_digits
        character(len=20)   :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(num))
      
        ! Count number of digits
        ! int_digits = int(floor(log10(abs(num)))) + 1
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        shifted_num = abs(num) / 10.0**n_digits
      
    end function shift_decimal_single


    !> @brief Function to validate the test information structure.
    !!        It sets the threshold value based on the threshold type.
    !!
    !! @param result_single The single precision result to be validated.
    !! @param result_double The double precision result to be validated.
    !! @param test_info The test information structure to be validated.
    !! @param arr_size The size of the array to be used for validation.
    !!
    subroutine psb_test_validate(result_single, result_double, test_info, arr_size, pass) 
        type(psb_test_info), intent(inout)  :: test_info
        real(psb_spk_), intent(in)          :: result_single
        real(psb_dpk_), intent(in)          :: result_double
        integer(psb_ipk_), intent(in)       :: arr_size 
        logical, intent(inout)              :: pass

        real(psb_spk_)                      :: local_single
        real(psb_dpk_)                      :: local_double
        integer(psb_ipk_)                   :: n
        real(psb_dpk_)                      :: gamma_n, unit_roundoff, delta, rel_err

        unit_roundoff = 5.96D-08 !! 1.11D-16

        delta = abs(result_double - real(result_single,psb_dpk_))
        rel_err = delta / abs(real(result_single,psb_dpk_))
        n = (arr_size / test_info%np) + (test_info%np - 1)


        ! write(psb_out_unit,'(A,F20.10)') "Computed delta:           ", delta

        if(test_info%threshold_type == VALUE) then
            ! Lower down values in order to match threshold for absolute error check
            local_single = shift_decimal_single(result_single)
            local_double = shift_decimal_double(result_double)

            delta = abs(result_double - real(result_single,psb_dpk_))

            if(delta < test_info%threshold) then
                pass = .true.
            else
                pass = .false.
            end if
        else if(test_info%threshold_type == GAMMA) then

            if(n * unit_roundoff >= 1) then
                write(test_info%output_unit,'(A)') "Error: Invalid GAMMA computation, n * U is greater than 1"
                write(*,*)
                call psb_test_exit(test_info)
            end if

            gamma_n = (n * unit_roundoff) / (1.0D0 - real(n * unit_roundoff, psb_dpk_) )
            test_info%threshold = gamma_n

            !! Handle case when result_single is zero
            if ((ieee_is_nan(rel_err)).and.(result_single == 0.0D0)) then
                rel_err = 0.0D0
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

    end subroutine


    !> @brief Subroutine to check the results of a single and double precision computation.
    !!        It compares the results and logs the outcome.
    !!
    !! @param result_single The single precision result to be checked.
    !! @param result_double The double precision result to be checked.
    !! @param test_info The test information structure containing the threshold and logging details.
    !! @param arr_size The size of the array used in the computation.
    !!
    subroutine psb_test_single_double_scalar_check(result_single, result_double, test_info, arr_size) 
        type(psb_test_info), intent(inout)  :: test_info
        real(psb_spk_), intent(inout)       :: result_single
        real(psb_dpk_), intent(inout)       :: result_double
        real(psb_dpk_)                      :: delta
        integer(psb_ipk_)                   :: int_digits, arr_size
        logical                             :: pass 
        character(len=64)                   :: out_string

        out_string = "Double precision check: "


        call psb_test_progress_bar(test_info)
        call psb_test_validate(result_single, result_double, test_info, arr_size, pass)

        if(pass .eqv. .true.) then 
            call psb_test_log_passed(test_info, out_string)
            test_info%success = test_info%success + 1  
        else
            call psb_test_log_failed(test_info, out_string)
            test_info%failure = test_info%failure + 1
            write(psb_out_unit,'(A,F20.10)') "Threshold used:           ", test_info%threshold
        end if
        write(psb_out_unit,'(A,F20.10)') "Single precision result:  ", result_single
        write(psb_out_unit,'(A,F20.10)') "Double precision result:  ", result_double
    end subroutine

    !> @brief Subroutine to check the results of a single and double precision vector computation.
    !!        It compares the results element-wise and logs the outcome.
    !! @param result_single The single precision vector result to be checked.
    !! @param result_double The double precision vector result to be checked.
    !! @param test_info The test information structure containing the threshold and logging details.
    !! @param arr_size The size of the array used in the computation.
    !!
    subroutine psb_test_single_double_vector_check(result_single, result_double, test_info, arr_size) 
        type(psb_test_info), intent(inout)          :: test_info
        real(psb_spk_), allocatable, intent(inout)  :: result_single(:)
        real(psb_dpk_), allocatable, intent(inout)  :: result_double(:)
        real(psb_dpk_)                              :: delta
        integer(psb_ipk_)                           :: int_digits, arr_size, i
        logical                                     :: pass 
        character(len=64)                           :: out_string

        out_string = "Double precision check: "
        pass = .true.

        call psb_test_progress_bar(test_info)
        do i = 1, size(result_single)
            call psb_test_validate(result_single(i), result_double(i), test_info, arr_size, pass)
            if(pass .eqv. .false. ) exit
        end do

        if(pass .eqv. .true.) then 
            call psb_test_log_passed(test_info, out_string)
            test_info%success = test_info%success + 1  
        else
            call psb_test_log_failed(test_info, out_string)
            test_info%failure = test_info%failure + 1
            write(psb_out_unit,'(A,I0)') "Comparison error occurred at index:  ", i
            write(psb_out_unit,'(A,F20.10)') "Single precision result:  ", result_single(i)
            write(psb_out_unit,'(A,F20.10)') "Double precision result:  ", result_double(i)        
        end if

        ! write(psb_out_unit,'(A,F20.10)') "Threshold used:           ", test_info%threshold
    end subroutine


    !> @brief Subroutine to check the global and local results of a single precision computation.
    !!        It compares the global result with the local result and logs the outcome.
    !!
    !! @param global_result_single The global single precision result to be checked.
    !! @param result_single The local single precision result to be checked.
    !! @param test_info The test information structure containing the logging details.
    !!
    subroutine psb_test_check_global_local(global_result_single, result_single, test_info)
        real(psb_spk_), intent(in)  :: global_result_single,result_single
        type(psb_test_info)         :: test_info
        character(len=64)           :: out_string

        out_string = "Global vs Local check: "

        call psb_test_progress_bar(test_info)
        
        ! Check if the global result is equal to the local result
        if (global_result_single == result_single) then
            call psb_test_log_passed(test_info, out_string)
            test_info%success = test_info%success + 1  
        else
            call psb_test_log_failed(test_info, out_string)
            test_info%failure = test_info%failure + 1
        end if
        write(test_info%output_unit,'(A,F20.10)') "Global single precision result: ", global_result_single
        write(test_info%output_unit,'(A,F20.10)') "Local single precision result:  ", result_single
    end subroutine


    !> @brief Subroutine to check the result of a single precision computation
    !!        against a previously saved result from a single process test.
    !!
    !! @param result_single The single precision result to be checked.
    !! @param test_info The test information structure containing the current test count.
    !!
    subroutine psb_test_process_check(result_single, test_info)
        real(psb_spk_), intent(inout)       :: result_single
        type(psb_test_info), intent(inout)  :: test_info
        real(psb_spk_)                      :: saved_result, local_saved, local_single
        integer(psb_ipk_)                   :: unit, info, file_size
        character(len=32)                   :: filename
        logical                             :: exists
        character(len=64)                   :: out_string

        out_string = "Multiprocess check: "
        call psb_test_progress_bar(test_info)

        ! Set the filename based on the test count
        write(filename, '(A,I0,A)') 'results/result_', test_info%current_test, '.txt'

        ! Check if the file exists
        inquire(file=trim(filename), exist=exists)
        if (.not.exists) then
            write(test_info%output_unit, '(A)') "Error: Result file does not exist."
            write(test_info%output_unit, '(A)') "Please ensure the single process test is run first to generate the result file."
            call psb_test_exit(test_info)
        end if

        ! Open the file for reading
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=info)

        ! Check if the file was opened successfully
        if (info /= 0) then
            write(*, '(A,I0)') "Error opening result file: ", info
            call psb_test_exit(test_info)
        end if

        ! Read the saved result
        read(unit, '(F20.10)') saved_result

        ! Close the file
        close(unit, iostat=info)

        local_saved = shift_decimal_single(saved_result)
        local_single = shift_decimal_single(result_single)

        ! Compare the saved result with the new result_single
        if (abs(local_saved - local_single) <= test_info%threshold) then
            call psb_test_log_passed(test_info, out_string)
            test_info%success = test_info%success + 1 
        else
            call psb_test_log_failed(test_info, out_string)
            test_info%failure = test_info%failure + 1
            write(test_info%output_unit, '(A,F20.10)') "Delta: ", abs(saved_result - result_single)
        end if

        write(test_info%output_unit, '(A,F20.10)') "Multi-process result: ", result_single
        write(test_info%output_unit, '(A,F20.10)') "Single process result: ", saved_result

    end subroutine


    subroutine psb_test_process_vector_check(result_single, test_info)
        real(psb_spk_), allocatable, intent(inout)  :: result_single(:)
        type(psb_test_info), intent(inout)          :: test_info
        real(psb_spk_), allocatable                 :: saved_result(:)
        integer(psb_ipk_)                           :: unit, info, file_size, int_digits, i
        character(len=32)                           :: filename
        logical                                     :: exists, pass
        character(len=64)                           :: out_string

        out_string = "Multiprocess check: "
        pass = .true.
        call psb_test_progress_bar(test_info)


        ! Set the filename based on the test count
        write(filename, '(A,I0,A)') 'results/result_', test_info%current_test, '.txt'

        ! Check if the file exists
        inquire(file=trim(filename), exist=exists)
        if (.not.exists) then
            write(test_info%output_unit, '(A)') "Error: Result file does not exist."
            write(test_info%output_unit, '(A)') "Please ensure the single process test is run first to generate the result file."
            call psb_test_exit(test_info)
        end if

        ! Open the file for reading
        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=info)

        ! Check if the file was opened successfully
        if (info /= 0) then
            write(*, '(A,I0)') "Error opening result file: ", info
            call psb_test_exit(test_info)
        end if

        ! Close the file
        close(unit, iostat=info)

        ! Read the saved result
        call mm_array_read(saved_result, info, filename=trim(filename))
        if (info /= 0) then
            write(*, '(A,I0)') "Error reading result file: ", info
            call psb_test_exit(test_info)
        end if

        do i = 1, size(result_single)
            ! call shift_decimal_single(saved_result(i),int_digits)
            ! call shift_decimal_single(result_single(i),int_digits)

            ! Compare the saved result with the new result_single
            if (abs(saved_result(i) - result_single(i)) > test_info%threshold) then
                pass = .false. 

                write(psb_out_unit,'(A,I0)') "Comparison error occurred at index:  ", i
                write(test_info%output_unit, '(F20.10,F20.10,A,L,A,L)') &
                & saved_result(i) - result_single(i), result_single(i) - saved_result(i), " ", & 
                & saved_result(i) - result_single(i) == 0, " ", result_single(i) - saved_result(i) == 0 
                write(test_info%output_unit, '(A,F20.10)') "Multi-process result: ", result_single(i)
                write(test_info%output_unit, '(A,F20.10)') "Single process result: ", saved_result(i)
                exit
            end if
        end do

        if(pass .eqv. .true.) then 
            call psb_test_log_passed(test_info, out_string)
            test_info%success = test_info%success + 1  
        else
            call psb_test_log_failed(test_info, out_string)
            test_info%failure = test_info%failure + 1
            


        end if




    end subroutine



end module psb_test_utils