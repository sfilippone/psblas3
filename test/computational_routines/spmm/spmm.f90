program main
    use psb_spmm_test
    use psb_base_mod
    use psb_util_mod

    implicit none

    ! matrix stats variables
    integer(psb_ipk_)       :: rows, cols

    ! MPI variables
    integer(psb_ipk_)       :: my_rank, np

    ! parameters array
    character(len=64)       :: x(4),y(4)  
    real(psb_spk_)          :: alpha(3), beta(3)  

    ! cycle indexes variables
    integer(psb_ipk_)       :: i,j,k,h,l
    character(len=256)      :: matrix_file
    logical                 :: matrix_exists
    integer(psb_ipk_)       :: info
    integer(psb_ipk_)       :: tests_number, count, last_percent

    ! sparse matrices
    type(psb_sspmat_type)           :: a
    type(psb_lsspmat_type)          :: aux_a

    ! matrix descriptor data structure
    type(psb_desc_type)             :: desc_a

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

    tests_number = size(x) * size(y) * size(alpha) * size(beta)
    count = 0
    last_percent = -1

    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    matrix_file = "matrix/1138_bus.mtx"
    inquire(file=matrix_file, exist=matrix_exists)
    if (.not.matrix_exists) then
        matrix_file = "../../comm/spmv/Geo_1438.mtx"
        inquire(file=matrix_file, exist=matrix_exists)
    end if
    if (.not.matrix_exists) then
        if (my_rank == psb_root_) then
            write(psb_out_unit,*) 'Matrix file not found. Expected matrix/1138_bus.mtx'
        end if
        call psb_abort(ctxt)
    end if

    if(my_rank == psb_root_) then
        write(psb_out_unit,*) 'Welcome to PSBLAS version: ',psb_version_string_
        write(psb_out_unit,*) 'This is the psb_spmm_test sample program'

        call read_matrix_market_size(matrix_file, rows, cols)

        call generate_vectors(rows,cols) 
    end if

    call psb_barrier(ctxt)

    !! Read and distribute matrix once
    call mm_mat_read(aux_a,info,filename=matrix_file)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error while reading matrix ", matrix_file
        call psb_abort(ctxt)
    end if

    rows = aux_a%get_nrows()
    cols = aux_a%get_ncols()

    call psb_matdist(aux_a, a, ctxt, desc_a, info, fmt="COO", parts=part_block)
    if(info /= psb_success_) then
        write(psb_out_unit,*) "Error while distributing matrix"
        call psb_abort(ctxt)
    end if

    !! 1138_bus matrix (sparse) - reuse distributed matrix across all parameter combinations
    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    call psb_spmm_kernel(a=a, desc_a=desc_a, rows=rows, cols=cols, &
                        & x_file=x(i), y_file=y(j), alpha=alpha(k), beta=beta(h), ctxt=ctxt)
                    if (my_rank == psb_root_) then
                        count = count + 1
                        call print_progress(count, tests_number, last_percent, "spmm testcases")
                    end if
                end do
            end do
        end do
    end do

    !! Deallocate matrix structures
    call psb_spfree(a, desc_a, info)
    call psb_cdfree(desc_a, info)


    call psb_exit(ctxt)

    return
contains

    subroutine print_progress(current, total, last_percent, label)
        implicit none
        integer(psb_ipk_), intent(in)    :: current, total
        integer(psb_ipk_), intent(inout) :: last_percent
        character(len=*), intent(in)     :: label
        integer(psb_ipk_)                :: percent, filled, width, i

        if (total <= 0) return
        percent = int(real(current) / real(total) * 100.0)
        if (percent == last_percent) return
        last_percent = percent

        width = 30
        filled = int(real(percent) / 100.0 * width)

        write(*,'(A)',advance='no') "[INFO]    Progress " // trim(label) // ": ["
        do i = 1, filled
            write(*,'(A)',advance='no') "#"
        end do
        do i = filled + 1, width
            write(*,'(A)',advance='no') "-"
        end do
        write(*,'(A,I3,A,I0,A,I0,A)',advance='no') "] ", percent, "% (", current, "/", total, ")"
        write(*,'(A)',advance='no') char(13)
        call flush(6)
        if (percent == 100) write(*,'(A)') ""
    end subroutine print_progress

end program main