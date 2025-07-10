!! Logger of the PSBLAS test suite. THis module should expose all the features 
!! usefull for the PSBLAS test developers, in order to generate an automatic log of all the tests.
!! This has been implemented as a .inc file so that it is possible to use only
!! use psb_test_utils in other files instead of multiple includes.
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!!


!> @brief Subroutine used to foramt correctly the message in the log about a test that has paseed
!!
!! @param test_info is a data structure containing all the metadata usefull for test output
!!
subroutine psb_test_log_passed(test_info)
    type(test_info_), intent(in)    :: test_info
    ! time stats variables
    character(len=8)                :: date  
    character(len=10)               :: time    
    character(len=5)                :: zones 
    integer                         :: values(8)

    call date_and_time(date, time, zones, values)

    write(test_info%output_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
    & A,A,A,I0,A,I0,T110,A)') &
    & values(1), values(2), values(3), values(5), values(6), values(7), &
    & "Generation ", trim(test_info%kernel_name), " single precision ", & 
    & test_info%current_test , "/", test_info%total_tests, "[OK]"

end subroutine

!> @brief Subroutine used to foramt correctly the message in the log about a test that has failed
!!
!! @param test_info is a data structure containing all the metadata usefull for test output
!!
subroutine psb_test_log_failed(test_info)
    type(test_info_), intent(in)    :: test_info
    ! time stats variables
    character(len=8)                :: date  
    character(len=10)               :: time
    character(len=5)                :: zones 
    integer(psb_ipk_)               :: values(8)

    call date_and_time(date, time, zones, values)

    write(test_info%output_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
    & A,A,A,I0,A,I0,T110,A)') &
    & values(1), values(2), values(3), values(5), values(6), values(7), &
    & "Generation ", trim(test_info%kernel_name), " single precision ", & 
    & test_info%current_test , "/", test_info%total_tests, "[FAIL]"

end subroutine

!> @brief Subroutine used to show the advancement of the test excecution as a progress bar.
!!        No output unit is specified since it is intended to be shown on the stdout.
!!
!! @param test_info is a data structure containing all the metadata usefull for test output
!!
subroutine psb_test_progress_bar(test_info)
    type(test_info_), intent(in)    :: test_info 
    integer(psb_ipk_)               :: percent, bar_width, num_hashes, num_spaces
    character(len=100)              :: bar

    bar_width = 40 

    ! Calculate percentage
    percent = int(real(test_info%current_test) / real(test_info%total_tests) * 100.0)

    ! Calculate how many '#' and '-' characters
    num_hashes = int(real(bar_width) * percent / 100.0)
    num_spaces = bar_width - num_hashes

    bar = '[' // repeat('#', num_hashes) // repeat('-', num_spaces) // ']'
    write(*,'(A,A,I3,A)', advance='no') char(13), 'Progress: '//bar//' ', percent, '%'
    call flush(6)
    
    ! Update stdout to get it back in normal mode when reaching the 100%
    if(test_info%current_test == test_info%total_tests) then
        write(*,'(A)', advance='no') char(13)//repeat(' ', 120)
        write(*,*)
        call flush(6)
    end if
end subroutine


