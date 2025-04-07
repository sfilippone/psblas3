! Test program for y = AX spmm routine
! Author: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
program psb_spmm_test
    use psb_base_mod
    use psb_prec_mod
    use psb_util_mod

    implicit none 

    ! input parameters
    character(len=40)           :: mtx_file, file_name, name

    ! sparse matrices
    type(psb_sspmat_type)       :: a
    type(psb_lsspmat_type)      :: aux_a

    ! vectors
    type(psb_s_vect_type)       :: x, y

    ! matrix descriptor data structure
    type(psb_desc_type)         :: desc_a

    ! communication context
    type(psb_ctxt_type)         :: ctxt
    integer(psb_ipk_)           :: my_rank, np, info, err_act

    ! matrix parameters
    integer(psb_ipk_)           :: m, n, nnz
    integer(psb_ipk_)           :: nr, nt ! In BLOCK ROWS distributin, the number of rows 


    real(psb_spk_), allocatable :: x_global(:)
    integer(psb_ipk_)           :: i

    name     = "psb_spmm_test"    ! Name of the program to output in case of error
    info     = psb_success_
    mtx_file = "matrix/1138_bus.mtx"


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
        call random_init(repeatable=.true.,image_distinct=.false.) 
        allocate(x_global(n))
        call random_number(x_global)
    end if

    ! nt = (m+np-1)/np
    ! nr = max(0,min(nt,m-(my_rank*nt)))

    ! Check if distribution metadata is correct
    ! nt = nr
    ! call psb_sum(ctxt,nt)
    ! if (nt /= m) then
    !   write(psb_err_unit,*) my_rank, 'Initialization error ',nr,nt,m
    !   info = -1
    !   call psb_barrier(ctxt)
    !   call psb_abort(ctxt)
    !   return
    ! end if

    ! call psb_cdall(ctxt,desc_a,info, nl=nr)
    ! if(info /= psb_success_) then
    !     write(psb_out_unit,*) "Error in dexcriptor allocator routine using BLOCK ROWS distribution"
    !     goto 9999
    ! end if

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
    
    call y%zero()


    ! y = alpha * A * x + betha * y
    call psb_spmm(sone,a,x,sone,y,desc_a,info)

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

    if(my_rank == psb_root_) then
        deallocate(x_global)
    end if
    
    call psb_exit(ctxt)
    stop


    ! Error handling
    9999 call psb_error(ctxt)
    
    call psb_errpush(info,name)
    call psb_error_handler(ctxt,err_act)
    call psb_exit(ctxt)
    stop

end program  