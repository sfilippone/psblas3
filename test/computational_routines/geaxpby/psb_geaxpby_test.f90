!> Test program for y = alpha * x + betha * y  psb_geaxbpy routine
!! Check the README.md to see all details about the tests.
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!! 
!! psb_geaxpby(alpha, x, beta, y, desc_a, info)
!!
!! Type: Synchronous.
!! 
!!          ======================================
!!          | Data type |      Precision         |
!!          ======================================
!!          | psb_spk_  | Short Precision Real   |
!!          | psb_dpk_  | Long Precision Real    |
!!          | psb_cpk_  | Short Precision Complex|
!!          | psb_zpk_  | Long Precision Complex |
!!          ======================================
!!                  Table 1: Data types
!! 
!!      ROUTINE PARAMETERS
!!
!! Input:
!!
!! alpha    Description: the scalar α.
!!          Scope: global
!!          Type: required
!!          Intent: in
!!          Specified as: a number of the data type indicated in Table 1.
!!
!! x        Description: the local portion of global dense matrix x.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: a rank one or two array or an object of type psb_T_vect_type
!!          containing numbers of type specified in Table 1. The rank of x must be
!!          the same of y.
!!
!! beta     Description: the scalar β.
!!          Scope: global
!!          Type: required
!!          Intent: in.
!!          Specified as: a number of the data type indicated in Table 1.
!! 
!! y        Description: the local portion of the global dense matrix y.
!!          Scope: local
!!          Type: required
!!          Intent: inout
!!          Specified as: a rank one or two array or an object of type psb_T_vect_type
!!          containing numbers of the type indicated in Table 1. The rank of y must
!!          be the same of x.
!!
!! desc_a   Description: contains data structures for communications.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: an object of type psb desc type.
!!
!! Output:
!!
!! y        Description: the local portion of the global dense matrix y.
!!          Scope: local
!!          Type: required
!!          Intent: inout
!!          Specified as: a rank one or two array or an object of type psb_T_vect_type
!!          containing numbers of the type indicated in Table 1. The rank of y must
!!          be the same of x.
!!
!! info     Description: Error code.
!!          Scope: local
!!          Type: required
!!          Intent: out.
!!          Specified as: An integer value; 0 means no error has been detected.
!!





program main
    use psb_base_mod
    use psb_util_mod 
    use psb_test_utils

    implicit none

    ! MPI variables
    integer(psb_ipk_)               :: my_rank, np

    ! Communicator variable
    type(psb_ctxt_type)             :: ctxt

    ! parameters array
    character(len=64)               :: x(4),y(4)  
    real(psb_dpk_)                  :: alpha(3), beta(3)
    integer(psb_ipk_)               :: arr_size  
    integer(psb_ipk_)               :: tests_number, count

    ! cycle indexes variables
    integer(psb_ipk_)               :: i,j,k,h,l
    integer(psb_ipk_)               :: info, ret, unit

    ! time stats variables
    character(len=8)                :: date    ! YYYYMMDD
    character(len=10)               :: time    ! HHMMSS.sss
    character(len=5)                :: zones    ! Time zone
    integer                         :: values(8)

    ! others
    character(len=:), allocatable   :: output_file_name
    type(test_info_)                :: test_info

    ! Initialize parameters
    x(1) = "vectors/x1.mtx"
    x(2) = "vectors/x2.mtx"
    x(3) = "vectors/x3.mtx"
    x(4) = "vectors/x4.mtx"

    y(1) = "vectors/y1.mtx"
    y(2) = "vectors/y2.mtx"
    y(3) = "vectors/y3.mtx"
    y(4) = "vectors/y4.mtx"

    alpha(1) = done
    alpha(2) = -done
    alpha(3) = dzero

    beta(1) = done
    beta(2) = -done
    beta(3) = dzero

    arr_size = 10000
    tests_number = size(x) * size(y) * size(alpha) * size(beta)
    count = 0
    
    !! Initialize test metadata
    test_info%total_tests = size(x) * size(y)
    test_info%threshold_type = GAMMA
    test_info%threshold = 0.0
    test_info%kernel_name = "psb_geaxpby"

    call psb_test_init(test_info)

    if(test_info%my_rank == psb_root_) then
        psb_out_unit = test_info%output_unit
        call psb_test_generate_input_vectors(arr_size) 
    end if

    call psb_bcast(test_info%ctxt,test_info%output_unit)
    call psb_barrier(test_info%ctxt)

    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Starting test excecution ..."
    !call psb_init(ctxt)
    !call psb_info(ctxt,my_rank,np)
    ! call generate_vectors(arr_size) 

    !call psb_bcast(ctxt,psb_out_unit)
    !call psb_barrier(ctxt)


    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Starting single precision computation..."

    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    
                    call psb_geaxpby_kernel(x_file=x(i), y_file=y(j), alpha = real(alpha(k),psb_spk_),& 
                    & beta = real(beta(h),psb_spk_), arr_size = arr_size, ctxt = test_info%ctxt, ret = ret, &
                    & output_file_name = output_file_name)

                    if(test_info%my_rank == psb_root_) then
                        count = count + 1
                        call date_and_time(date, time, zones, values)

                        if(ret /= -1) then 
                            ! Success formatted output
                            write(psb_out_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
                            & A,A,A,I0,A,I0,T110,A)') &
                            & values(1), values(2), values(3), values(5), values(6), values(7), &
                            & "Generation geaxpby single precision result file ", & 
                            & output_file_name , ' ', count , "/", tests_number, "[OK]"
                        else
                            ! Fail formatted output
                            write(psb_out_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
                            & A,A,A,I0,A,I0,T110,A)') &
                            & values(1), values(2), values(3), values(5), values(6), values(7), &
                            & "Generation geaxpby single precision result file ", & 
                            & output_file_name , ' ', count , "/", tests_number, "[FAIL]"
                            goto 9998
                        end if
                    end if
                    call psb_barrier(test_info%ctxt)
                end do
            end do
        end do
    end do

    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Single precision computation completed succesfully!"


    if(test_info%my_rank == psb_root_) then
        write(psb_out_unit, *) ''
        count = 0
    end if


    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Starting double precision check..."



    call psb_barrier(test_info%ctxt)

    ! Here double precision comparison should be done
    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    call psb_geaxpby_check(x_file=x(i), y_file=y(j), alpha = alpha(k), beta = beta(h), & 
                    & arr_size = arr_size, ctxt = test_info%ctxt, ret = ret, output_file_name = output_file_name)
                    
                    if(my_rank == psb_root_) then
                        count = count + 1
                        call date_and_time(date, time, zones, values)

                        if(ret == 0) then 
                            ! Success formatted output
                            write(psb_out_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
                            & A,A,A,I0,A,I0,T110,A)') &
                            & values(1), values(2), values(3), values(5), values(6), values(7), &
                            & "Double precision check on file ", & 
                            & output_file_name , ' ', count , "/", tests_number, "[OK]"
                        else
                            ! Fail formatted output
                            write(psb_out_unit,'("[", I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"] ",& 
                            & A,A,A,I0,A,I0,T110,A)') &
                            & values(1), values(2), values(3), values(5), values(6), values(7), &
                            & "Double precision check on file ", & 
                            & output_file_name , ' ', count , "/", tests_number, "[FAIL]"
                            write(psb_out_unit,'(A,I0)') "[ERROR]   Error at element ", abs(ret)
                            write(psb_out_unit,'(A,F10.8)') "Alpha:", alpha(k)
                            write(psb_out_unit,'(A,F15.8)') "Beta: ", beta(h)

                            goto 9999
                        end if
                    end if
                    call psb_barrier(test_info%ctxt)
                end do
            end do
        end do
    end do

    if(test_info%my_rank == psb_root_) then 
        write(*,'(A)') "[INFO]    Duble precision check completed succesfully!"
        close(unit)
    end if
     
    call psb_exit(test_info%ctxt)
    return 

    9998 continue
    if(test_info%my_rank == psb_root_) then
        close(unit) 
        write(*,'(A,I0,A,I0,A)') "[ERROR]   Error in geaxpby single precision computation ", & 
        & count, "/", tests_number, " see log file for details"
    end if

    9999 continue
    if(test_info%my_rank == psb_root_) then
        close(unit) 
        write(*,'(A,I0,A,I0,A)') "[ERROR]   Error in geaxpby double precision check ", &
        & count, "/", tests_number, " see log file for details"
    end if
    
    call psb_test_exit(test_info)
    return



    contains

    !> @brief Function to excecute psb_geaxpby in single precision and
    !!        save the results on file
    !!
    subroutine psb_geaxpby_kernel(x_file, y_file, alpha, beta, arr_size, ctxt, ret, output_file_name)
        ! implicit none 

        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
        real(psb_spk_), intent(in)                  :: alpha, beta
        integer(psb_ipk_), intent(in)               :: arr_size
        type(psb_ctxt_type), intent(in)             :: ctxt

        ! output parameters
        integer(psb_ipk_), intent(out)              :: ret
        character(len=:), allocatable, intent(out)  :: output_file_name       

        ! vectors
        type(psb_s_vect_type)           	        :: x, y

        ! matrix descriptor data structure
        type(psb_desc_type)             	        :: desc_a

        ! communication context
        integer(psb_ipk_)               	        :: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_spk_), allocatable     	        :: x_global(:), y_global(:)
        integer(psb_ipk_)               	        :: i

        ! others
        logical                                     :: exists



        info = psb_success_

        call psb_info(ctxt,my_rank,np)

        if (my_rank < 0) then 
            ! This should not happen, but just in case
            call psb_error(ctxt) 
        endif

        ! Generate random array for b using always the same seed
        if(my_rank == psb_root_) then
            allocate(x_global(arr_size))
            allocate(y_global(arr_size))
            call mm_array_read(x_global,info,filename=x_file)
            call mm_array_read(y_global,info,filename=y_file)
        end if

        ! Allocate descriptor as if it was a block rows distribution
        call psb_cdall(ctxt, desc_a, info,nl=arr_size/np)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating desc_a data structure"
            goto 9999
        end if

        call psb_cdasb(desc_a, info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error assembling desc_a data structure"
            goto 9999
        end if


        call psb_geall(x,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating x data structure"
            goto 9999
        end if


        ! Populate x class using data from x_global vector
        call psb_scatter(x_global,x,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate x data structure"
            goto 9999
        end if


        call psb_geall(y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating y data structure"
            goto 9999
        end if

        ! Populate y class using data from y_global vector
        call psb_scatter(y_global,y,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate y data structure"
            goto 9999
        end if


        ! y = alpha * x + beta * y
        call psb_geaxpby(alpha,x,beta,y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_geaxpby routine"
            goto 9999
        end if

        ! Make the root process be the one that saves everything on file
        if(np == 1) then 
            ! Check if output directory exists
            inquire(file='serial/', exist=exists)
            if (.not.exists) then
                call system('mkdir serial/')
            end if
            output_file_name = "serial/"
        else 
            ! Check if output directory exists
            inquire(file='parallel/', exist=exists)
            if (.not.exists) then
                call system('mkdir parallel/')
            end if
            output_file_name = "parallel/"
        end if

        output_file_name = output_file_name // "sol_" // x_file(9:10) // "_" // y_file(9:10)

        if(alpha == sone) then 
            output_file_name = output_file_name // "_a1"
        else if(alpha == -sone) then 
            output_file_name = output_file_name // "_a2"
        else if(alpha == szero) then
            output_file_name = output_file_name // "_a3"
        end if

        if(beta == sone) then 
            output_file_name = output_file_name // "_b1.mtx"
        else if(beta == -sone) then 
            output_file_name = output_file_name // "_b2.mtx"
        else if(beta == szero) then
            output_file_name = output_file_name // "_b3.mtx"
        end if

        ! gather the result combining all the partial ones
        call psb_gather(y_global, y, desc_a, info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error gathering global vector x to write on file"
            goto 9999
        end if

        ! Save result to output file
        if(my_rank == psb_root_) then
            call mm_array_write(y_global,"Result vector",info,filename=output_file_name)
        end if

        ! Deallocate
        call psb_gefree(x, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in vector x free routine"
            goto 9999
        end if

        call psb_gefree(y, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in vector y free routine"
            goto 9999
        end if

        call psb_cdfree(desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in matrix descriptor free routine"
            goto 9999
        end if

        if(my_rank == 0) then
            deallocate(x_global)
            deallocate(y_global)
        end if

        return


        ! Error handling
        9999 ret = -1 
        stop

    end subroutine



    !> @brief Function to excecute psb_geaxpby in double precision and
    !!        compare the results with the ones on file
    !!
    subroutine psb_geaxpby_check(x_file, y_file, alpha, beta, arr_size, ctxt, ret, output_file_name)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
        real(psb_dpk_), intent(in)                  :: alpha, beta
        integer(psb_ipk_), intent(in)               :: arr_size
        type(psb_ctxt_type), intent(in)             :: ctxt

        ! output parameters
        integer(psb_ipk_), intent(out)              :: ret
        character(len=:), allocatable, intent(out)  :: output_file_name      
        ! vectors
        type(psb_d_vect_type)           	        :: x, y
        type(psb_s_vect_type)           	        :: y_check

        ! matrix descriptor data structure
        type(psb_desc_type)             	        :: desc_a

        ! communication context
        integer(psb_ipk_)               	        :: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_dpk_), allocatable     	        :: x_global(:), y_global(:) 
        integer(psb_ipk_)               	        :: i

        ! others
        logical                                     :: exists
 



        info = psb_success_

        call psb_info(ctxt,my_rank,np)

        if (my_rank < 0) then 
            ! This should not happen, but just in case
            call psb_error(ctxt) 
        endif

        ! Generate random array for b using always the same seed
        if(my_rank == psb_root_) then
            allocate(x_global(arr_size))
            allocate(y_global(arr_size))
            call mm_array_read(x_global,info,filename=x_file)
            call mm_array_read(y_global,info,filename=y_file)
        end if

        ! Allocate descriptor as if it was a block rows distribution
        call psb_cdall(ctxt, desc_a, info,nl=10000/np)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating desc_a data structure"
            goto 9999
        end if

        call psb_cdasb(desc_a, info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error assembling desc_a data structure"
            goto 9999
        end if


        call psb_geall(x,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating x data structure"
            goto 9999
        end if


        ! Populate x class using data from x_global vector
        call psb_scatter(x_global,x,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate x data structure"
            goto 9999
        end if


        call psb_geall(y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating y data structure"
            goto 9999
        end if

        ! Populate y class using data from y_global vector
        call psb_scatter(y_global,y,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate y data structure"
            goto 9999
        end if


        call psb_geall(y_check,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating y_check data structure"
            goto 9999
        end if


        ! y = alpha * x + beta * y
        call psb_geaxpby(alpha,x,beta,y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_geaxpby routine"
            goto 9999
        end if


        ! gather the result combining all the partial ones
        call psb_gather(y_global, y, desc_a, info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error gathering global vector y used for comparison"
            goto 9999
        end if

        if(my_rank == psb_root_) then 
            ! Make the root process be the one that saves everything on file
            if(np == 1) then 
                ! Check if output directory exists
                inquire(file='serial/', exist=exists)
                if(.not.exists) then
                    write(psb_out_unit,'(A)') "Error in psb_geaxpby_check routine, no single precision result is saved on file"
                    goto 9999
                end if
                output_file_name = "serial/"
            else 
                ! Check if output directory exists
                inquire(file='parallel/', exist=exists)
                if(.not.exists) then
                    write(psb_out_unit,'(A)') "Error in psb_geaxpby_check routine, no single precision result is saved on file"
                    goto 9999
                end if
                output_file_name = "parallel/"
            end if
    
            output_file_name = output_file_name // "sol_" // x_file(9:10) // "_" // y_file(9:10)
    
            if(alpha == done) then 
                output_file_name = output_file_name // "_a1"
            else if(alpha == -done) then 
                output_file_name = output_file_name // "_a2"
            else if(alpha == dzero) then
                output_file_name = output_file_name // "_a3"
            end if
    
            if(beta == done) then 
                output_file_name = output_file_name // "_b1.mtx"
            else if(beta == -done) then 
                output_file_name = output_file_name // "_b2.mtx"
            else if(beta == dzero) then
                output_file_name = output_file_name // "_b3.mtx"
            end if
    

            ! Read single precision result from file
            call mm_array_read(y_check,info,filename=output_file_name)
            if(info /= psb_success_) then
                write(psb_out_unit,'(A)') "Error in mm_array_read for y_check data structure"
                goto 9999
            end if
            
            ! 5.96e-08 is 2^-24 (Single precision unit roundoff)
            ! 1.19e-07 is 2^-23 (Single precision unit interval)
            do i=1, arr_size
                !! write(*, *) abs(y_global(i) - y_check%v%v(i)) > 5.96D-08, &
                !! & y_global(i), y_check%v%v(i), output_file_name
                if(abs(y_global(i) - y_check%v%v(i)) > 1.19e-07) then
                   ret = -i
                   write(psb_out_unit, '(A,F10.8)') "Y computed in double precision: ", y_global(i) 
                   write(psb_out_unit, '(A,F10.8)') "Y read from single precision file: ", y_check%v%v(i) 
                   write(psb_out_unit, '(A,F10.8)') "Diff: ", abs(y_global(i) - y_check%v%v(i)) 
                   exit
                end if
            end do
        end if 
        
        ! Deallocate
        call psb_gefree(x, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in vector x free routine"
            goto 9999
        end if

        call psb_gefree(y, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in vector y free routine"
            goto 9999
        end if

        call psb_gefree(y_check, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in vector y_check free routine"
            goto 9999
        end if


        call psb_cdfree(desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in matrix descriptor free routine"
            goto 9999
        end if

        if(my_rank == 0) then
            deallocate(x_global)
            deallocate(y_global)
        end if

        ret = 0

        return


        ! Error handling
        9999 ret = -1 
        stop

    end subroutine





end program main
