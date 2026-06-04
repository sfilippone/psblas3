!> Test program for y = x^T * y or y = x^H * y psb_gedot routine
!! Check the README.md to see all details about the tests.
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!! 
!! psb_gedot(x, y, desc_a, info [,global])
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
!! x        Description: the local portion of global dense matrix x.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: a rank one or two array or an object of type psb_T_vect_type
!!          containing numbers of type specified in Table 1. The rank of x must be
!!          the same of y.
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
!! global   Descritption: Specifies whether the computation should include the global 
!!          reduction across all processes.
!!          Scope: global
!!          Type: optional
!!          Intent: in
!!          Specified as: a logical scalar. 
!!          Default: global=.true.
!!
!! Output:
!!
!! Function value   the dot product of vectors x and y.
!!                  Scope: global unless the optional variable global=.false. 
!1                  has been specified
!!                  Specified as: a number of the data type indicated in Table 1.
!!
!! info     Description: Error code.
!!          Scope: local
!!          Type: required
!!          Intent: out
!!          Specified as: An integer value; 0 means no error has been detected.
!!
!!
!! NOTES
!! 
!! 1.       The computation of a global result requires a global communication, which
!!          entails a significant overhead. It may be necessary and/or advisable to
!!          compute multiple dot products at the same time; in this case, it is possible
!!          to improve the runtime efficiency by using the following scheme:
!!
!!              vres(1) = psb_gedot(x1,y1,desc_a,info,global=.false.)
!!              vres(2) = psb_gedot(x2,y2,desc_a,info,global=.false.)
!!              vres(3) = psb_gedot(x3,y3,desc_a,info,global=.false.)
!!              call psb_sum(ctxt,vres(1:3))
!!
!!          In this way the global communication, which for small sizes is a latency-
!!          bound operation, is invoked only once.
!!
module psb_gedot_test
    use psb_base_mod
    use psb_util_mod 

    contains

    !> @brief Function to excecute psb_geaxpby in single precision and
    !!        save the results on file
    !!
    subroutine psb_gedot_kernel(x_file, y_file, arr_size, ctxt, ret, output_file_name)

        implicit none 

        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
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
        real(psb_spk_)                              :: result(1)


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


        ! y = x^T * y
        result(1) = psb_gedot(x,y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedot routine"
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

        output_file_name = output_file_name // "sol_" // x_file(9:10) // "_" // y_file(9:10) // ".mtx"

        ! Save result to output file
        if(my_rank == psb_root_) then
            call mm_array_write(result,"Result of the scalar product computation",info,filename=output_file_name)
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

        ret = 0
        return


        ! Error handling
        9999 ret = -1 
        stop

    end subroutine



    !> @brief Function to excecute psb_geaxpby in double precision and
    !!        compare the results with the ones on file
    !!
    subroutine psb_gedot_check(x_file, y_file, arr_size, ctxt, ret, output_file_name)

        implicit none 

        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
        integer(psb_ipk_), intent(in)               :: arr_size
        type(psb_ctxt_type), intent(in)             :: ctxt

        ! output parameters
        integer(psb_ipk_), intent(out)              :: ret
        character(len=:), allocatable, intent(out)  :: output_file_name      
        ! vectors
        type(psb_d_vect_type)           	        :: x, y
        type(psb_s_vect_type)           	        :: result_check

        ! matrix descriptor data structure
        type(psb_desc_type)             	        :: desc_a

        ! communication context
        integer(psb_ipk_)               	        :: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_dpk_), allocatable     	        :: x_global(:), y_global(:)
        integer(psb_ipk_)               	        :: i

        ! others
        logical                                     :: exists
        real(psb_dpk_)                              :: result(1)


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


        call psb_geall(result_check,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating y_check data structure"
            goto 9999
        end if


        ! y = x^T * y
        result(1) = psb_gedot(x,y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedot routine"
            goto 9999
        end if

        if(my_rank == psb_root_) then 
            ! Make the root process be the one that saves everything on file
            if(np == 1) then 
                ! Check if output directory exists
                inquire(file='serial/', exist=exists)
                if(.not.exists) then
                    write(psb_out_unit,'(A)') "Error in psb_gedot_check routine, no single precision result is saved on file"
                    goto 9999
                end if
                output_file_name = "serial/"
            else 
                ! Check if output directory exists
                inquire(file='parallel/', exist=exists)
                if(.not.exists) then
                    write(psb_out_unit,'(A)') "Error in psb_gedot_check routine, no single precision result is saved on file"
                    goto 9999
                end if
                output_file_name = "parallel/"
            end if
    
            output_file_name = output_file_name // "sol_" // x_file(9:10) // "_" // y_file(9:10) // ".mtx"

            ! Read single precision result from file
            call mm_array_read(result_check,info,filename=output_file_name)
            if(info /= psb_success_) then
                write(psb_out_unit,'(A)') "Error in mm_array_read for y_check data structure"
                goto 9999
            end if

            ! 5.96e-08 is 2^-24 (Single precision unit roundoff)
            ! 1.19e-07 is 2^-23 (Single precision unit interval)
            !! call shift_decimal_double(result(1))
            !! call shift_decimal_single(result_check%v%v(1))

            ! write(*,*) result(1),result_check%v%v(1), (arr_size * 1.19D-07) / (done-arr_size * 1.19D-07)
            if(abs(result(1) - result_check%v%v(1)) > (arr_size * 1.19D-07) / (done-arr_size * 1.19D-07)) then
               ret = -1 
               return
            end if
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

        call psb_gefree(result_check, desc_a,info)
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

        ret = 0
        return


        ! Error handling
        9999 ret = -1 
        stop

    end subroutine

    subroutine shift_decimal_double(n)
        
        implicit none

        real(psb_dpk_),intent(inout)    :: n
        integer                         :: n_digits
        character(len=20)               :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(n))
      
        ! Count number of digits
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        n = abs(n) / 10.0**n_digits
      
    end subroutine


    subroutine shift_decimal_single(n)
        
        implicit none

        real(psb_spk_),intent(inout)    :: n
        integer                         :: n_digits
        character(len=20)               :: int_str
      
      
        ! Convert the absolute value of the integer part to string
        write(int_str, '(I0)') int(abs(n))
      
        ! Count number of digits
        n_digits = len_trim(adjustl(int_str))
      
        ! Shift the decimal point
        n = abs(n) / 10.0**n_digits
      
    end subroutine


    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    subroutine generate_vectors(arr_size)


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

end module psb_gedot_test