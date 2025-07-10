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
program main
    use psb_base_mod
    use psb_util_mod 
    use psb_test_utils

    implicit none



    ! Communicator variable
    type(psb_ctxt_type)             :: ctxt

    ! parameters array
    character(len=64)               :: x(4),y(4)  
    integer(psb_ipk_)               :: arr_size  
    integer(psb_ipk_)               :: tests_number, count

    ! cycle indexes variables
    integer(psb_ipk_)               :: i,j,k,h,l
    integer(psb_ipk_)               :: info, unit

    ! results
    real(psb_spk_)                  :: result_single
    real(psb_dpk_)                  :: result_double
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

    arr_size = 100000

    !! Initialize test metadata
    test_info%total_tests = size(x) * size(y)
    test_info%threshold_type = VALUE
    test_info%threshold = 0.0
    test_info%kernel_name = "psb_gedot"

    call psb_test_init(test_info)

    if(test_info%my_rank == psb_root_) then
        psb_out_unit = test_info%output_unit
        call psb_test_generate_input_vectors(arr_size) 
    end if

    call psb_bcast(test_info%ctxt,test_info%output_unit)
    call psb_barrier(test_info%ctxt)


    if(test_info%my_rank == psb_root_) write(*,'(A)') "[INFO]    Starting test excecution ..."

    ! Iterate over test parameters
    do i=1,size(x)
        do j=1,size(y)
            call psb_gedot_real_kernel(x(i), y(j), arr_size, test_info%ctxt, result_single, result_double)
            
            if(test_info%my_rank == psb_root_) then
                call psb_test_single_double_check(result_single,result_double,test_info)
                test_info%current_test = test_info%current_test + 1 
            end if
            call psb_barrier(test_info%ctxt)            
        end do
    end do
     
    call psb_test_exit(test_info)


contains

    !> @brief Function to excecute psb_gedot in single precision real
    !!        vector and compare with the same computation in double 
    !!        precision
    !!
    !! @param x_file file name of the first vector
    !! @param y_file file name of the second vector
    !! @param arr_size size of the vectors
    !! @param ctxt communication context
    !! @param result_single result of the single precision computation
    !! @param result_double result of the double precision computation
    !! 
    subroutine psb_gedot_real_kernel(x_file, y_file, arr_size, ctxt, result_single, result_double)
        ! input parameters
        character(len = *), intent(in)  :: x_file, y_file
        integer(psb_ipk_), intent(in)   :: arr_size
        type(psb_ctxt_type), intent(in) :: ctxt

        ! output parameters
        real(psb_spk_), intent(out)     :: result_single
        real(psb_dpk_), intent(out)     :: result_double

        ! vectors
        type(psb_s_vect_type)           :: x_single, y_single
        type(psb_d_vect_type)           :: x_double, y_double

        ! matrix descriptor data structure
        type(psb_desc_type)             :: desc_a

        ! communication context
        integer(psb_ipk_)               :: my_rank, np, info, err_act

        ! variables outside PSLBALS data structures
        real(psb_spk_), allocatable     :: x_single_global(:), y_single_global(:)
        real(psb_dpk_), allocatable     :: x_double_global(:), y_double_global(:)
        integer(psb_ipk_)               :: i

        ! others
        logical                         :: exists

        info = psb_success_

        call psb_info(ctxt,my_rank,np)

        if (my_rank < 0) then 
            ! This should not happen, but just in case
            call psb_error(ctxt) 
        endif

        ! Generate random array for b using always the same seed
        if(my_rank == psb_root_) then
            allocate(x_single_global(arr_size))
            allocate(y_single_global(arr_size))
            allocate(x_double_global(arr_size))
            allocate(y_double_global(arr_size))

            call mm_array_read(x_single_global,info,filename=x_file)
            call mm_array_read(y_single_global,info,filename=y_file)
            call mm_array_read(x_double_global,info,filename=x_file)
            call mm_array_read(y_double_global,info,filename=y_file)
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
            goto 9999
        end if

        call psb_scatter(x_double_global,x_double,desc_a,info,root=psb_root_)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_scatter to populate double precision x data structure"
            goto 9999
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

        ! y = x^T * y
        result_single = psb_gedot(x_single,y_single,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedot routine in single precision"
            goto 9999
        end if


        result_double = psb_gedot(x_double,y_double,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedot routine in double precision"
            goto 9999
        end if

        ! Deallocate
        9999 call psb_gefree(x_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector x free routine"
            goto 9999
        end if

        call psb_gefree(y_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector y free routine"
            goto 9999
        end if

        call psb_gefree(x_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector x free routine"
            goto 9999
        end if

        call psb_gefree(y_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector y free routine"
            goto 9999
        end if


        call psb_cdfree(desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in matrix descriptor free routine"
            goto 9999
        end if

        if(my_rank == 0) then
            deallocate(x_single_global)
            deallocate(y_single_global)
            deallocate(x_double_global)
            deallocate(y_double_global)
        end if

        return

    end subroutine

end program main