program main
    use psb_spmm_test
    use psb_base_mod

    implicit none

    ! matrix stats variables
    integer(psb_ipk_)       :: rows, cols

    ! MPI variables
    integer(psb_ipk_)       :: my_rank, np

    ! parameters array
    character(len=64)       :: x(4),y(4)  
    real(psb_ipk_)          :: alpha(3), beta(3)  

    ! cycle indexes variables
    integer(psb_ipk_)       :: i,j,k,h,l

    ! Communicator variable
    type(psb_ctxt_type)     :: ctxt

    ! Initialize parameters
    x(1) = "vectors/x1.mtx"
    x(2) = "vectors/x2.mtx"
    x(3) = "vectors/x3.mtx"
    x(4) = "vectors/x4.mtx"

    y(1) = "vectors/y1.mtx"
    y(2) = "vectors/y2.mtx"
    y(3) = "vectors/y3.mtx"
    y(4) = "vectors/y4.mtx"

    alpha(1) = sone
    alpha(2) = -sone
    alpha(3) = szero

    beta(1) = sone
    beta(2) = -sone
    beta(3) = szero

    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    if(my_rank == psb_root_) then
        write(psb_out_unit,*) 'Welcome to PSBLAS version: ',psb_version_string_
        write(psb_out_unit,*) 'This is the psb_spmm_test sample program'

        call read_matrix_market_size("matrix/1138_bus.mtx", rows, cols)

        call generate_vectors(rows,cols) 
    end if

    call psb_barrier(ctxt)



    !! 1138_bus matrix (sparse)
    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    call psb_spmm_kernel(mtx_file="matrix/1138_bus.mtx",x_file=x(i), y_file=y(j), &
                        & alpha = alpha(k), beta = beta(h), ctxt = ctxt)
                end do
            end do
        end do
    end do


    call psb_exit(ctxt)

    return
end program main