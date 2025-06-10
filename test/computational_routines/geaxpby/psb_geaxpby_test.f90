!> Test program for y = alpha * x + betha * y  psb_geaxbpy routine
!! Check the README.md to see all details about the tests.
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!! 
!! Type: Synchronous.
!! 
!! ROUTINE PARAMETERS
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
!!          Specified as: a rank one or two array or an object of type psb T vect type
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
!!          Specified as: a rank one or two array or an object of type psb T vect type
!!          containing numbers of the type indicated in Table 1. The rank of y must
!!          be the same of x.
!!
!! desc_a   Description: contains data structures for communications.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: an object of type psb desc type.
!!
!! info     Description: Error code.
!!          Scope: local
!!          Type: required
!!          Intent: out.
!!          Specified as: An integer value; 0 means no error has been detected.
!!
module psb_geaxpby_test

    contains

    !> @brief Function to excecute psb_geaxpby in single precision and
    !!        save the results on file
    !!
    subroutine psb_geaxpby_kernel(x_file, y_file, alpha, beta, arr_size, ctxt, ret)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        ! input parameters
        character(len = *), intent(in)      :: x_file, y_file
        real(psb_spk_), intent(in)          :: alpha, beta
        integer(psb_ipk_), intent(in)       :: arr_size
        type(psb_ctxt_type), intent(in)     :: ctxt

        ! output parameters
        integer(psb_ipk_), intent(out)       :: ret

        ! vectors
        type(psb_s_vect_type)           	:: x, y

        ! matrix descriptor data structure
        type(psb_desc_type)             	:: desc_a

        ! communication context
        integer(psb_ipk_)               	:: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_spk_), allocatable     	:: x_global(:), y_global(:)
        integer(psb_ipk_)               	:: i

        ! others
        logical                             :: exists
        character(len=:), allocatable       :: output_file_name       



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
        if(np > 1) then 
            call psb_gather(x_global, x, desc_a, info)
            if(info /= psb_success_) then
                write(psb_out_unit,'(A)') "Error gathering global vector x to write on file"
                goto 9999
            end if
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
    subroutine psb_geaxpby_check(x_file, y_file, alpha, beta, arr_size, ctxt, ret)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        ! input parameters
        character(len = *), intent(in)      :: x_file, y_file
        real(psb_dpk_), intent(in)          :: alpha, beta
        integer(psb_ipk_), intent(in)       :: arr_size
        type(psb_ctxt_type), intent(in)     :: ctxt

        ! output parameters
        integer(psb_ipk_), intent(out)      :: ret

        ! vectors
        type(psb_d_vect_type)           	:: x, y
        type(psb_s_vect_type)           	:: y_check

        ! matrix descriptor data structure
        type(psb_desc_type)             	:: desc_a

        ! communication context
        integer(psb_ipk_)               	:: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_dpk_), allocatable     	:: x_global(:), y_global(:)
        integer(psb_ipk_)               	:: i

        ! others
        logical                             :: exists
        character(len=:), allocatable       :: output_file_name       



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

            ! 5.96e-08 is 2^-24 (Single precision unit roundoff)
            ! 1.19e-07 is 2^-23 (Single precision unit interval)
            do i=1, arr_size
                ! write(psb_out_unit, *) abs(y%v%v(i) - real(y%v%v(i),psb_spk_)), 1.19D-07, &
                ! & abs(y%v%v(i) - real(y%v%v(i),psb_spk_)) > 1.19D-07
                if(abs(y%v%v(i) - y_check%v%v(i)) > 1.19D-07) then
                   ret = -1 
                end if
            end do
        end if 

        call psb_barrier(ctxt)
        
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



    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    subroutine generate_vectors(arr_size)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        integer(psb_ipk_), intent(in)               :: arr_size
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

end module psb_geaxpby_test