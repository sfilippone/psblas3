!! Utils module containing all tools usefull to do test validation and result checking
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!!
module psb_test_utils
    use psb_base_mod
    use psb_util_mod

    implicit none

    ! Define the enumeration values to represent testing criteria
    integer, parameter :: VALUE = 1

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
    end type test_info_


contains

    include 'psb_test_log.inc'

    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    subroutine generate_vectors(arr_size)
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


    subroutine psb_test_single_double_check(result_single, result_double, test_info)
        integer(psb_ipk_)               :: int_digits 
        real(psb_spk_), intent(inout)   :: result_single
        real(psb_dpk_), intent(inout)   :: result_double
        real(psb_dpk_)                  :: delta
        type(test_info_), intent(inout) :: test_info


        call psb_test_progress_bar(test_info)

        call shift_decimal_single(result_single,int_digits)
        call shift_decimal_double(result_double,int_digits)

        delta = abs(result_double - result_single)

        if(delta < test_info%threshold) then 
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