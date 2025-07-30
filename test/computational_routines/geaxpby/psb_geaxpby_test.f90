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

    ! others
    character(len=:), allocatable   :: output_file_name
    type(psb_test_info)             :: test_info
    
    ! result vectors
    real(psb_spk_), allocatable     :: y_single(:)
    real(psb_dpk_), allocatable     :: y_double(:)


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

    !! Initialize test metadata
    test_info%total_tests = size(x) * size(y) * size(alpha) * size(beta)
    test_info%threshold_type = VALUE
    test_info%threshold = 1.0D-06
    test_info%kernel_name = "psb_geaxpby"

    call psb_test_init(test_info)

    if(test_info%my_rank == psb_root_) then
        psb_out_unit = test_info%output_unit
        call psb_test_generate_input_vectors(arr_size) 
        allocate(y_single(arr_size))
        allocate(y_double(arr_size))
    end if

    call psb_bcast(test_info%ctxt,test_info%output_unit)
    call psb_barrier(test_info%ctxt)

    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Starting test excecution ..."

    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)

                    call psb_geaxpby_real_kernel(x(i), y(j), alpha(k),beta(h), arr_size, test_info%ctxt, y_single, y_double)

                    if(test_info%my_rank == psb_root_) then
                        if(test_info%np > 1) then 
                            ! If the program is being run on multiple processes, we need to
                            ! check the result on the root process with the one computed only using 
                            ! a single process
                            call psb_test_process_vector_check(y_single, test_info)
                        else
                            call psb_test_single_double_vector_check(y_single,y_double,test_info, arr_size) 
                            
                            ! If the program is being run on a single process, we can save the result directly
                            call psb_test_save_vector_result(y_single, test_info)
                        end if

                        test_info%current_test = test_info%current_test + 1                        
                    end if
                    call psb_barrier(test_info%ctxt)
                end do
            end do
        end do
    end do
    
    call psb_test_exit(test_info)

contains

    !> @brief This subroutine implements the psb_geaxpby kernel for real vectors 
    !!        performing the operation y = alpha * x + beta * y. It reads input from files, 
    !!        performs the operation, and outputs the result.
    !!
    !! @param x_file: file name of the input vector x
    !! @param y_file: file name of the input vector y
    !! @param alpha: scalar alpha
    !! @param beta: scalar beta
    !! @param arr_size: size of the vectors
    !! @param ctxt: communication context
    !! @param y_single_global: single precision output vector
    !! @param y_double_global: double precision output vector
    !!
    subroutine psb_geaxpby_real_kernel(x_file, y_file, alpha, beta, arr_size, ctxt, y_single_global, y_double_global)
        ! implicit none 

        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
        real(psb_dpk_), intent(in)                  :: alpha, beta
        integer(psb_ipk_), intent(in)               :: arr_size
        type(psb_ctxt_type), intent(in)             :: ctxt

        ! vectors
        type(psb_s_vect_type)                       :: x_single, y_single
        type(psb_d_vect_type)                       :: x_double, y_double 

        ! matrix descriptor data structure
        type(psb_desc_type)                         :: desc_a

        ! communication context
        integer(psb_ipk_)                           :: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_spk_), allocatable                 :: x_single_global(:) 
        real(psb_spk_), allocatable, intent(inout)  :: y_single_global(:)
        real(psb_dpk_), allocatable                 :: x_double_global(:) 
        real(psb_dpk_), allocatable, intent(inout)  :: y_double_global(:)
        integer(psb_ipk_)                           :: i, nl

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
            allocate(x_single_global(arr_size))
            allocate(x_double_global(arr_size))

            call mm_array_read(x_single_global,info,filename=x_file)
            call mm_array_read(y_single_global,info,filename=y_file)
            call mm_array_read(x_double_global,info,filename=x_file)
            call mm_array_read(y_double_global,info,filename=y_file)
        end if

        ! Allocate descriptor as if it was a block rows distribution
        nl = (arr_size)/np + mod(arr_size,np)

        call psb_bcast(test_info%ctxt,test_info%output_unit)
        call psb_barrier(test_info%ctxt)


        call psb_cdall(ctxt, desc_a, info,nl=nl)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating desc_a data structure"
            goto 9999
        end if

        call psb_cdasb(desc_a, info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error assembling desc_a data structure"
            goto 9999
        end if


        call psb_geall(x_single,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating single precision x data structure"
            goto 9999
        end if

        call psb_geall(x_double,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating double precision x data structure"
            goto 9999
        end if

        ! Populate x class using data from x_global vector
        call psb_scatter(x_single_global,x_single,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate single precision x data structure"
            return 
        end if

        call psb_scatter(x_double_global,x_double,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate double precision x data structure"
            return 
        end if

        call psb_geall(y_single,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating single precision y data structure"
            goto 9999
        end if

        call psb_geall(y_double,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error allocating double precision y data structure"
            goto 9999
        end if

        ! Populate y class using data from y_global vector
        call psb_scatter(y_single_global,y_single,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate single precision y data structure"
            goto 9999
        end if

        call psb_scatter(y_double_global,y_double,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate double precision y data structure"
            goto 9999
        end if


        ! y = alpha * x + beta * y
        call psb_geaxpby(real(alpha,psb_spk_),x_single,real(beta,psb_spk_),y_single,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_geaxpby routine"
            goto 9999
        end if

        ! y = alpha * x + beta * y
        call psb_geaxpby(alpha,x_double,beta,y_double,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_geaxpby routine"
            goto 9999
        end if

        ! Gather final result
        call psb_gather(y_single_global, y_single, desc_a, info, psb_root_) 
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gather routine"
            goto 9999
        end if
    
        call psb_gather(y_double_global, y_double, desc_a, info, psb_root_) 
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gather routine"
            goto 9999
        end if


        ! Deallocate
        9999 call psb_gefree(x_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector x free routine"
        end if

        call psb_gefree(y_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector y free routine"
        end if

        call psb_gefree(x_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector x free routine"
        end if

        call psb_gefree(y_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector y free routine"
        end if


        call psb_cdfree(desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in matrix descriptor free routine"
        end if

        if(my_rank == 0) then
            deallocate(x_single_global)
            deallocate(x_double_global)
        end if

        return


    end subroutine



end program main
