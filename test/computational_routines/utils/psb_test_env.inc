!> @brief Function to initialize the test environment.
!!        It is used to set the output unit and the kernel name.
!!
!! @param test_info The test information structure to be initialized.
!!
subroutine psb_test_init(test_info)
    type(psb_test_info), intent(inout)  :: test_info
    integer(psb_ipk_)                   :: output_unit, info
    ! MPI variables
    integer(psb_ipk_)                   :: my_rank


    call psb_init(test_info%ctxt)
    call psb_info(test_info%ctxt,test_info%my_rank,test_info%np)

    ! Check if the kernel name is set
    if (trim(test_info%kernel_name) == "default") then
        write(*, '(A)') "Error: Kernel name is not set. Please set the kernel name before running the test."
        call psb_exit(test_info%ctxt)
    end if

    ! Set the output unit to stdout by default
    if(test_info%np == 1) then
        open(newunit=output_unit, file=trim(test_info%kernel_name)//'_test.log', &
        & status='replace', action='write', iostat=info)
    else
        open(newunit=output_unit, file=trim(test_info%kernel_name)//'_test.log', &
        & status='old', action='write', position='append', iostat=info)
    end if

    ! Check if the file was opened successfully
    if (info /= 0) then
       write(*, '(A,I0)') "Error opening log file for kernel ", test_info%kernel_name
       write(*, '(A,I0)') "I-O Status Code: ", info
       write(*, '(A)') "Please check if the file is accessible and writable."
       call psb_test_exit(test_info)
    end if

    test_info%output_unit = output_unit

    write(test_info%output_unit,'(A,A)') 'Welcome to PSBLAS version: ',psb_version_string_
    write(test_info%output_unit,'(A)') 'This is the psb_gedot_test sample program'
    write(test_info%output_unit,'(A,I0)') 'Number of processes used in this computation: ', test_info%np
    write(test_info%output_unit,'(A)') ''

end subroutine

!> @brief Function to finalize the test environment, it is used to close the output unit 
!!        and finalize the PSBLAS context.
!!
!! @param test_info The test information structure to be finalized.
!!
subroutine psb_test_exit(test_info)
    type(psb_test_info), intent(inout)  :: test_info
    integer(psb_ipk_)                   :: info

    ! Finalize test
    if(test_info%my_rank == psb_root_) then
        write(test_info%output_unit, *) ""

        ! Check if the kernel name is set to default
        if(test_info%kernel_name == "default") then
            write(test_info%output_unit,'(A)') "Warning: Kernel name is not set. Using default kernel name 'default'."
        else
            write(test_info%output_unit,'(A)') "Kernel name set to: " // trim(test_info%kernel_name)
        end if

        ! Check the threshold type and value
        if(test_info%mat_dist == DEFAULT) then
            write(test_info%output_unit,'(A)') "Matrix distribution set to default (block row distribution)."
        else
            write(test_info%output_unit,'(A,I0)') "Matrix distribution set to: ", test_info%mat_dist
        end if
        
        ! Check the threshold type and value
        if(test_info%threshold_type == DEFAULT) then
            write(test_info%output_unit,'(A,G20.10)') "Threshold type is set to default&
            & using single precision IEEE unit roundoff as threshould ", 5.96D-08 
        else if(test_info%threshold_type == VALUE) then
            write(test_info%output_unit,'(A,G20.10)') "Threshold type is set to VALUE, so absolute error will be checked&
            & using threshold: ", test_info%threshold
        else if(test_info%threshold_type == GAMMA) then
            write(test_info%output_unit,'(A, G20.10)') "Threshold type is set to GAMMA, so relative error will& 
            & be checked using threshold: ", test_info%threshold 
        end if
        write(test_info%output_unit,'(A)') ''
            
        write(*,'(A)') "[INFO]    Tests completed successfully!"
        write(*,'(A,I0)') "          Passed: ", test_info%success
        write(*,'(A,I0)') "          Failed: ", test_info%failure
        write(*,'(A,I0)') "          Total:  ", test_info%total_tests
        write(*,'(A)') "[INFO]    Check " // trim(test_info%kernel_name) // "_test.log for a full description"

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