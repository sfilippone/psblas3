!> Test program for y = AX spsb_pmm routine
!! Check the README.md to see all details about the tests.
!!
!! Author: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
program psb_spmm_test
    use psb_base_mod
    use psb_util_mod

    implicit none 

    ! input parameters
    character(len=256)              :: mtx_file, file_name, name

    ! Testing parameters
    character(len=256),dimension(2) :: mtx_files

    ! sparse matrices
    type(psb_sspmat_type)           :: a
    type(psb_lsspmat_type)          :: aux_a

    ! vectors
    type(psb_s_vect_type)           :: x, y

    ! matrix descriptor data structure
    type(psb_desc_type)             :: desc_a

    ! communication context
    type(psb_ctxt_type)             :: ctxt
    integer(psb_ipk_)               :: my_rank, np, info, err_act

    ! matrix parameters
    integer(psb_ipk_)               :: m, n, nnz
    integer(psb_ipk_)               :: nr, nt ! In BLOCK ROWS distributin, the number of rows 


    real(psb_spk_), allocatable     :: x_global(:), y_global(:)
    integer(psb_ipk_)               :: i

    name     = "psb_spmm_test"    ! Name of the program to output in case of error
    info     = psb_success_
    mtx_file = "matrix/1138_bus.mtx"
    mtx_files(1) = "matrix/1138_bus.mtx" 


    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    if (my_rank < 0) then 
        ! This should not happen, but just in case
        call psb_error(ctxt) 
    endif

    if (my_rank == psb_root_) then 
        write(psb_out_unit,*) 'Welcome to PSBLAS version: ',psb_version_string_
        write(psb_out_unit,*) 'This is the ',trim(name),' sample program'
    end if


    call mm_mat_read(aux_a,info,filename=mtx_file)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error while reading matric ", mtx_file
        goto 9999        
    end if 

    
    ! part_block it's a macro defined in psb_blockpart_mod to identify BLOCK ROWS distribution
    call psb_matdist(aux_a, a, ctxt,desc_a,info,fmt="COO",parts=part_block) 
     

    m   = aux_a%get_nrows()
    n   = aux_a%get_ncols()
    nnz = aux_a%get_nzeros()

    call psb_bcast(ctxt,m)
    call psb_bcast(ctxt,n)
    call psb_bcast(ctxt,nnz)

    if(my_rank == psb_root_) then 
        write(psb_out_unit,*) "Matrix stats"
        write(psb_out_unit,*) "ROWS:", m
        write(psb_out_unit,*) "COLS:", n
        write(psb_out_unit,*) "NNZ: ", nnz
    end if

    ! Generate random array for b using always the same seed
    if(my_rank == psb_root_) then
        allocate(x_global(n))
        allocate(y_global(n))
        call generate_vectors(n) ! True for x
        call mm_array_read(x_global,info,filename="vectors/x1.mtx")
        call mm_array_read(y_global,info,filename="vectors/y1.mtx")
    end if

    call psb_geall(x,desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error allocating x data structure"
        goto 9999
    end if

    ! Populate x class using data from x_global vector
    call psb_scatter(x_global,x,desc_a,info,root=psb_root_)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in psb_scatter to populate x data structure"
        goto 9999
    end if

    
    call psb_geall(y,desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error allocating y data structure"
        goto 9999
    end if

    ! Populate y class using data from y_global vector
    call psb_scatter(y_global,y,desc_a,info,root=psb_root_)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in psb_scatter to populate y data structure"
        goto 9999
    end if


    ! y = alpha * A * x + betha * y
    call psb_spmm(sone,a,x,sone,y,desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in psb_spmm routine"
        goto 9999
    end if

    ! Save result to output file
    call mm_array_write(y,"Result vector",info,filename="sol_m1_x1_y1.mtx")
    
    ! Deallocate
    call psb_gefree(x, desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in vector x free routine"
        goto 9999
    end if

    call psb_gefree(y, desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in vector y free routine"
        goto 9999
    end if

    call psb_spfree(a, desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in matrix A free routine"
        goto 9999
    end if
    
    call psb_cdfree(desc_a,info)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error in matrix descriptor free routine"
        goto 9999
    end if

    deallocate(x_global)
    deallocate(y_global)
    
    call psb_exit(ctxt)
    stop


    ! Error handling
    9999 call psb_error(ctxt)
    
    call psb_errpush(info,name)
    call psb_error_handler(ctxt,err_act)
    call psb_exit(ctxt)
    stop


contains

    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    !! @param n     The size of the vector.
    subroutine generate_vectors(n)
        use psb_base_mod

        implicit none 

        integer(psb_ipk_), intent(in)               :: n
        integer(psb_ipk_), allocatable              :: x(:), y(:)
        real(psb_spk_), allocatable                 :: v(:)
        integer(psb_ipk_)                           :: i

        allocate(x(n))
        allocate(y(n))
        allocate(v(2*n))

        
        call random_init(repeatable=.true.,image_distinct=.true.) 
        call random_number(v)

        do i = 1,n 
            x(i) = int(v(i) * 100)
        end do

        do i = 1,n
            y(i) = int(v(i+n) * 100)
        end do


        ! Write only positive in x_1
        call mm_array_write(x,"Positive vector",info,filename="vectors/x1.mtx")
        call mm_array_write(y,"Positive vector",info,filename="vectors/y1.mtx")

        ! Write only negative in x_2
        do i=1,n  
            x(i) = -x(i)
            y(i) = -y(i)
        end do  

        call mm_array_write(x,"Negative vector",info,filename="vectors/x2.mtx")
        call mm_array_write(y,"Negative vector",info,filename="vectors/y2.mtx")


        ! Since numbers are less than one and always positive, we have to generate negative ones subtractiong 50
        do i=1,n
            x(i) = -x(i) ! Make the values positive again  
            x(i) = x(i) - 50
            y(i) = -y(i) ! Make the values positive again  
            y(i) = y(i) - 50
        end do   

        ! Write random in x_3
        call mm_array_write(x,"Random vector",info,filename="vectors/x3.mtx")
        call mm_array_write(y,"Random vector",info,filename="vectors/y3.mtx")

        ! Write zero in x_4
        do i=1,n
            x(i) = 0
            y(i) = 0
        end do   

        call mm_array_write(x,"Null vector",info,filename="vectors/x4.mtx")
        call mm_array_write(y,"Null vector",info,filename="vectors/y4.mtx")

        deallocate(x)
        deallocate(y)
        deallocate(v)

    end subroutine

end program  

