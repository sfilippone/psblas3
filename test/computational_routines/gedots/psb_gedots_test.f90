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
!!          Specified as: an object of type psb_desc_type.
!!
!! Output:
!!
!! res      Description: is the dot product of vectors x and y
!!          Scope: global
!!          Intent: out
!!          Specified as: a number or a rank-one array of the data type indicated in Table 1
!!
!! info     Description: Error code.
!!          Scope: local
!!          Type: required
!!          Intent: out
!!          Specified as: An integer value; 0 means no error has been detected.
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
    real(psb_spk_), allocatable     :: result_single_vector(:)
    real(psb_dpk_), allocatable     :: result_double_vector(:)
    type(psb_test_info)             :: test_info



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
    test_info%threshold_type = GAMMA
    test_info%threshold = 1.0D-06
    test_info%kernel_name = "psb_gedots"

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
            call psb_gedots_real_kernel(x(i), y(j), arr_size, test_info%ctxt,result_single, result_double)

            if(test_info%my_rank == psb_root_) then                
                if(test_info%np > 1) then 
                    ! If the program is being run on multiple processes, we need to
                    ! check the result on the root process with the one computed only using 
                    ! a single process
                    call psb_test_process_check(result_single, test_info)
                else
                    call psb_test_single_double_scalar_check(result_single,result_double,test_info, arr_size)
                    
                    ! If the program is being run on a single process, we can save the result directly
                    call psb_test_save_result(result_single, test_info)
                end if
                test_info%current_test = test_info%current_test + 1
            end if

            call psb_barrier(test_info%ctxt)
        end do
    end do


    ! Test using matrices
    x(1) = "../matrix/1138_bus.mtx"
    x(2) = "../matrix/crystm03.mtx"
    x(3) = "../matrix/qc2534.mtx"
    x(4) = "../matrix/rdb5000.mtx"

    y(1) = "../matrix/1138_bus.mtx"
    y(2) = "../matrix/crystm03.mtx"
    y(3) = "../matrix/qc2534.mtx"
    y(4) = "../matrix/rdb5000.mtx"



    if(test_info%my_rank == psb_root_) then
        allocate(result_single_vector(arr_size))
        allocate(result_double_vector(arr_size))
    end if

    ! Iterate over test parameters
    !! do i=1,size(x)
    !!     do j=1,size(y)
    !!         call psb_gedots_real_matrix_kernel(x(i), y(j), arr_size, test_info%ctxt,result_single_vector, result_double_vector)
!! 
    !!         if(test_info%my_rank == psb_root_) then
    !!             if(test_info%np > 1) then 
    !!                 ! If the program is being run on multiple processes, we need to
    !!                 ! check the result on the root process with the one computed only using 
    !!                 ! a single process
    !!                 call psb_test_process_vector_check(result_single_vector, test_info)
    !!             else
    !!                 call psb_test_single_double_vector_check(result_single_vector,result_double_vector,test_info, arr_size)
    !!                 
    !!                 ! If the program is being run on a single process, we can save the result directly
    !!                 call psb_test_save_vector_result(result_single_vector, test_info)
    !!             end if
    !!             test_info%current_test = test_info%current_test + 1
    !!         end if
!! 
    !!         call psb_barrier(test_info%ctxt)
    !!     end do
    !! end do

    if(test_info%my_rank == psb_root_) then
        deallocate(result_single_vector)
        deallocate(result_double_vector)
    end if
     
    call psb_test_exit(test_info)

contains


    subroutine psb_gedots_real_matrix_kernel(x_file, y_file, arr_size, ctxt, result_single, result_double)
        ! input parameters
        character(len = *), intent(in)              :: x_file, y_file
        integer(psb_ipk_), intent(in)               :: arr_size
        type(psb_ctxt_type), intent(in)             :: ctxt

        ! output parameters
        real(psb_spk_), allocatable, intent(inout)  :: result_single(:)
        real(psb_dpk_), allocatable, intent(inout)  :: result_double(:)
        
        ! matrices
        type(psb_lsspmat_type)                      :: x_single_aux, y_single_aux
        type(psb_ldspmat_type)                      :: x_double_aux, y_double_aux

        ! sparse matrices
        type(psb_sspmat_type)                       :: x_single, y_single
        type(psb_dspmat_type)                       :: x_double, y_double
        
        ! matrix parameters
        integer(psb_ipk_)                           :: x_rows, x_cols, x_nnz
        integer(psb_ipk_)                           :: y_rows, y_cols, y_nnz
        integer(psb_ipk_)                           :: nl

        ! matrix descriptor data structure
        type(psb_desc_type)                         :: desc_a

        ! communication context
        integer(psb_ipk_)                           :: my_rank, np, info




                ! Allocate dense arrays for the local part only
        integer(psb_ipk_) :: num_local_rows_x, num_local_cols_x, num_local_rows_y, num_local_cols_y
        real(psb_spk_), allocatable :: x_dense_single(:,:), y_dense_single(:,:)
        real(psb_dpk_), allocatable :: x_dense_double(:,:), y_dense_double(:,:)
        ! For single precision x
        integer(psb_ipk_), allocatable :: row_x(:), col_x(:)
        real(psb_spk_), allocatable   :: val_x(:)
        integer(psb_ipk_)             :: nnz_x, idx
        type(psb_s_coo_sparse_mat)        :: x_coo_single
        ! For single precision y
        integer(psb_ipk_), allocatable :: row_y(:), col_y(:)
        real(psb_spk_), allocatable   :: val_y(:)
        integer(psb_ipk_)             :: nnz_y
        type(psb_s_coo_sparse_mat)        :: y_coo_single

        ! For double precision x
        integer(psb_ipk_), allocatable :: row_xd(:), col_xd(:)
        real(psb_dpk_), allocatable   :: val_xd(:)
        integer(psb_ipk_)             :: nnz_xd
        type(psb_d_coo_sparse_mat)        :: x_coo_d
        ! For double precision y
        integer(psb_ipk_), allocatable :: row_yd(:), col_yd(:)
        real(psb_dpk_), allocatable   :: val_yd(:)
        integer(psb_ipk_)             :: nnz_yd
        type(psb_d_coo_sparse_mat)        :: y_coo_d



        class(psb_s_coo_sparse_mat), pointer :: a_ptr_s
        class(psb_d_coo_sparse_mat), pointer :: a_ptr_d






        call psb_info(ctxt,my_rank,np)
        
        if (my_rank < 0) then 
            ! This should not happen, but just in case
            call psb_error(ctxt) 
        endif

        call mm_mat_read(a=x_single_aux,info=info, iunit=17, filename=x_file)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error reading single precision matrix x"
            return
        end if
        
        call mm_mat_read(a=x_double_aux,info=info, iunit=17, filename=x_file)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error reading double precision matrix x"
            return
        end if

        call mm_mat_read(a=y_single_aux,info=info, iunit=17, filename=y_file)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error reading single precision matrix y"
            return
        end if

        call mm_mat_read(a=y_double_aux,info=info, iunit=17, filename=y_file)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error reading double precision matrix y"
            return
        end if


        x_rows    = x_single_aux%get_nrows()
        x_cols    = x_single_aux%get_ncols()
        x_nnz     = x_single_aux%get_nzeros()

        y_rows    = y_single_aux%get_nrows()
        y_cols    = y_single_aux%get_ncols()
        y_nnz     = y_single_aux%get_nzeros()


        ! Allocate descriptor as if it was a block rows distribution
        nl = (arr_size)/np + mod(arr_size,np)


        ! part_block it's a macro defined in psb_blockpart_mod to identify BLOCK ROWS distribution
        call psb_matdist(x_single_aux, x_single, ctxt,desc_a,info,fmt="COO",parts=part_block) 
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_matdist for single precision matrix x"
            goto 9999
        end if

        call psb_matdist(x_double_aux, x_double, ctxt,desc_a,info,fmt="COO",parts=part_block)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_matdist for double precision matrix x"
            goto 9999
        end if

        call psb_matdist(y_single_aux, y_single, ctxt,desc_a,info,fmt="COO",parts=part_block)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_matdist for single precision matrix y"
            goto 9999
        end if  

        call psb_matdist(y_double_aux, y_double, ctxt,desc_a,info,fmt="COO",parts=part_block)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_matdist for double precision matrix y"
            goto 9999
        end if



        

        ! --------------------------------------------------------------------- !



        ! Get local dimensions from the descriptor
        num_local_rows_x = x_single%get_nrows()
        num_local_cols_x = x_single%get_ncols()
        num_local_rows_y = y_single%get_nrows()
        num_local_cols_y = y_single%get_ncols()

        allocate(x_dense_single(num_local_rows_x, num_local_cols_x))
        allocate(y_dense_single(num_local_rows_y, num_local_cols_y))
        allocate(x_dense_double(num_local_rows_x, num_local_cols_x))
        allocate(y_dense_double(num_local_rows_y, num_local_cols_y))

        ! Initialize to zero
        x_dense_single = 0.0_psb_spk_
        y_dense_single = 0.0_psb_spk_
        x_dense_double = 0.0_psb_dpk_
        y_dense_double = 0.0_psb_dpk_

        ! Explicitly extract COO values and insert them in the dense arrays




        select type(a_ptr_s => x_single%a)
        type is (psb_s_coo_sparse_mat)
            x_coo_single = a_ptr_s
        class default
            write(psb_out_unit,'(A)') "Error: x_single%a is not of type psb_s_coo_sparse_mat"
            return
        end select

        nnz_x = x_coo_single%get_nzeros()
        allocate(row_x(nnz_x), col_x(nnz_x), val_x(nnz_x))
        row_x = x_coo_single%ia
        col_x = x_coo_single%ja
        val_x = x_coo_single%val

        do idx = 1, nnz_x
            x_dense_single(row_x(idx), col_x(idx)) = val_x(idx)
        end do



        select type(a_ptr_s => y_single%a)
        type is (psb_s_coo_sparse_mat)
            y_coo_single = a_ptr_s
        class default
            write(psb_out_unit,'(A)') "Error: x_single%a is not of type psb_s_coo_sparse_mat"
            return
        end select        
        
        nnz_y = y_coo_single%get_nzeros()
        allocate(row_y(nnz_y), col_y(nnz_y), val_y(nnz_y))
        row_y = y_coo_single%ia
        col_y = y_coo_single%ja
        val_y = y_coo_single%val
        do idx = 1, nnz_y
            y_dense_single(row_y(idx), col_y(idx)) = val_y(idx)
        end do



        select type(a_ptr_d => x_double%a)
        type is (psb_d_coo_sparse_mat)
            x_coo_d = a_ptr_d
        class default
            write(psb_out_unit,'(A)') "Error: x_single%a is not of type psb_s_coo_sparse_mat"
            return
        end select      


        nnz_xd = x_coo_d%get_nzeros()
        allocate(row_xd(nnz_xd), col_xd(nnz_xd), val_xd(nnz_xd))
        row_xd = x_coo_d%ia
        col_xd = x_coo_d%ja
        val_xd = x_coo_d%val
        
        do idx = 1, nnz_xd
            x_dense_double(row_xd(idx), col_xd(idx)) = val_xd(idx)
        end do





        select type(a_ptr_d => y_double%a)
        type is (psb_d_coo_sparse_mat)
            y_coo_d = a_ptr_d
        class default
            write(psb_out_unit,'(A)') "Error: x_single%a is not of type psb_s_coo_sparse_mat"
            return
        end select      


        nnz_yd = y_coo_d%get_nzeros()
        allocate(row_yd(nnz_yd), col_yd(nnz_yd), val_yd(nnz_yd))
        row_yd = y_coo_d%ia
        col_yd = y_coo_d%ja
        val_yd = y_coo_d%val
        do idx = 1, nnz_yd
            y_dense_double(row_yd(idx), col_yd(idx)) = val_yd(idx)
        end do




        ! --------------------------------------------------------------------- !

        ! y = x^T * y
        call psb_gedots(result_single,x_dense_single,y_dense_single,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedots routine in single precision"
            goto 9999
        end if
!! 
        !! ! y = x^T * y
        call psb_gedots(result_double,x_dense_double,y_dense_double,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedots routine in double precision"
            goto 9999
        end if








        9999 call psb_spfree(x_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector x free routine"
        end if

        call psb_spfree(y_single, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in single precision vector y free routine"
        end if

        call psb_spfree(x_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector x free routine"
        end if

        call psb_spfree(y_double, desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in double precision vector y free routine"
        end if

        call psb_cdfree(desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in matrix descriptor free routine"
        end if


    end subroutine


    !> @brief Function to excecute psb_gedots in single precision real
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
    subroutine psb_gedots_real_kernel(x_file, y_file, arr_size, ctxt, result_single, result_double)
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
        integer(psb_ipk_)               :: i, nl

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
        nl = (arr_size)/np + mod(arr_size,np)

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
        call psb_gedots(result_single,x_single%v%v,y_single%v%v,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedots routine in single precision"
            goto 9999
        end if


        call psb_gedots(result_double,x_double%v%v,y_double%v%v,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,'(A)') "Error in psb_gedots routine in double precision"
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
            deallocate(y_single_global)
            deallocate(x_double_global)
            deallocate(y_double_global)
        end if

    end subroutine

end program main