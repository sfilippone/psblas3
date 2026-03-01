!> Test program for y = AX psb_spmm routine
!! Check the README.md to see all details about the tests.
!!
!! Author: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
module psb_spmm_test
    
    contains

    subroutine psb_spmm_kernel(mtx_file,x_file, y_file, alpha, beta, ctxt)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        ! input parameters
        character(len = *), intent(in)      :: mtx_file, x_file, y_file
        real(psb_spk_), intent(in)          :: alpha, beta

        character(len=:), allocatable       :: output_file_name       

        ! sparse matrices
        type(psb_sspmat_type)           	:: a
        type(psb_lsspmat_type)          	:: aux_a

        ! vectors
        type(psb_s_vect_type)           	:: x, y

        ! matrix descriptor data structure
        type(psb_desc_type)             	:: desc_a

        ! communication context
        type(psb_ctxt_type), intent(in)     :: ctxt
        integer(psb_ipk_)               	:: my_rank, np, info, err_act

        ! matrix parameters
        integer(psb_ipk_)               	:: rows, cols, nnz
        integer(psb_ipk_)               	:: nr, nt ! In BLOCK ROWS distributin, the number of rows 

        ! variables outside PSLBALS data structures
        real(psb_spk_), allocatable     	:: x_global(:), y_global(:)
        integer(psb_ipk_)               	:: i

        info = psb_success_

        call psb_info(ctxt,my_rank,np)

        if (my_rank < 0) then 
            ! This should not happen, but just in case
            call psb_error(ctxt) 
        endif

        call mm_mat_read(aux_a,info,filename=mtx_file)
        if(info /= psb_success_) then
            write(psb_out_unit,*) "Error while reading matrix ", mtx_file
            goto 9999        
        end if 


        ! part_block it's a macro defined in psb_blockpart_mod to identify BLOCK ROWS distribution
        call psb_matdist(aux_a, a, ctxt,desc_a,info,fmt="COO",parts=part_block) 


        rows    = aux_a%get_nrows()
        cols    = aux_a%get_ncols()
        nnz     = aux_a%get_nzeros()

        call psb_bcast(ctxt,rows)
        call psb_bcast(ctxt,cols)
        call psb_bcast(ctxt,nnz)

        ! Generate random array for b using always the same seed
        if(my_rank == psb_root_) then
            allocate(x_global(cols))
            allocate(y_global(rows))
            call mm_array_read(x_global,info,filename=x_file)
            call mm_array_read(y_global,info,filename=y_file)
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


        ! y = alpha * A * x + beta * y
        call psb_spmm(alpha,a,x,beta,y,desc_a,info)
        if(info /= psb_success_) then
            write(psb_out_unit,*) "Error in psb_spmm routine"
            goto 9999
        end if

        ! Make the root process be the one that saves everything on file
        if(np == 1) then 
            output_file_name = "serial/"
        else 
            output_file_name = "parallel/"
        end if

        ! Insert matrix number in output file naming convention
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

        ! Save result to output file
        call mm_array_write(y,"Result vector",info,filename=output_file_name)

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

        if(my_rank == 0) then
            deallocate(x_global)
            deallocate(y_global)
        end if

        return


        ! Error handling
        9999 call psb_error(ctxt)

        call psb_error_handler(ctxt,err_act)
        stop
    end subroutine

    !> @brief Function to randomly generate x and y vectors 
    !!        and save them on multiple files based on their
    !!        coefficients values.
    !!
    subroutine generate_vectors(rows, cols)
        use psb_base_mod
        use psb_util_mod

        implicit none 

        integer(psb_ipk_), intent(in)               :: rows, cols
        real(psb_spk_), allocatable                 :: x(:), y(:)
        integer(psb_ipk_)                           :: i, info

        allocate(x(rows))
        allocate(y(cols))
        
        call random_init(repeatable=.true.,image_distinct=.true.) 
        call random_number(x)
        call random_number(y)

        ! Write only positive in x_1
        call mm_array_write(x,"Positive vector",info,filename="vectors/x1.mtx")
        call mm_array_write(y,"Positive vector",info,filename="vectors/y1.mtx")

        ! Write only negative in x_2
        do i=1,rows 
            x(i) = -x(i)
        end do  

        do i=1,cols
            y(i) = -y(i)
        end do

        call mm_array_write(x,"Negative vector",info,filename="vectors/x2.mtx")
        call mm_array_write(y,"Negative vector",info,filename="vectors/y2.mtx")


        ! Since numbers are less than one and always positive, we have to generate negative ones subtractiong 50
        do i=1,rows
            x(i) = -x(i) ! Make the values positive again  
            x(i) = x(i) - 0.5
        end do   

        do i=1,cols
            y(i) = -y(i) ! Make the values positive again  
            y(i) = y(i) - 0.5
        end do

        ! Write random in x_3
        call mm_array_write(x,"Random vector",info,filename="vectors/x3.mtx")
        call mm_array_write(y,"Random vector",info,filename="vectors/y3.mtx")

        ! Write zero in x_4
        do i=1,rows
            x(i) = 0
        end do   

        do i=1,cols
            y(i) = 0
        end do

        call mm_array_write(x,"Null vector",info,filename="vectors/x4.mtx")
        call mm_array_write(y,"Null vector",info,filename="vectors/y4.mtx")

        deallocate(x)
        deallocate(y)

    end subroutine


    subroutine read_matrix_market_size(filename, rows, cols)
        implicit none
        character(len=*), intent(in)    :: filename
        integer, intent(out)            :: rows, cols
        integer                         :: ret, nnz
        character(len=256)              :: line
        logical                         :: found
        integer                         :: unit
      
        ! Open the file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ret)
        if (ret /= 0) then
           print *, 'Error opening file: ', filename
           stop
        end if
      
        found = .false.
      
        ! Skip comment lines (starting with %)
        do
           read(unit, '(A)', iostat=ret) line
           if (ret /= 0) exit
           if (line(1:1) /= '%') then
              read(line, *) rows, cols, nnz
              found = .true.
              exit
           end if
        end do
      
        if (.not. found) then
           print *, 'Error: header not found in Matrix Market file.'
           stop
        end if
      
        close(unit)
      end subroutine read_matrix_market_size

end module psb_spmm_test

