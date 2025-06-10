program main
    use psb_geaxpby_test
    use psb_base_mod

    implicit none

    ! MPI variables
    integer(psb_ipk_)       :: my_rank, np

    ! Communicator variable
    type(psb_ctxt_type)     :: ctxt

    ! parameters array
    character(len=64)       :: x(4),y(4)  
    real(psb_dpk_)          :: alpha(3), beta(3)
    integer(psb_ipk_)       :: arr_size  
    integer(psb_ipk_)       :: tests_number, count

    ! cycle indexes variables
    integer(psb_ipk_)       :: i,j,k,h,l
    integer(psb_ipk_)       :: info, ret, unit

    ! Setup logger output
    open(unit, file='psblas_geaxpby_test.log', status='replace', action='write', iostat=info)
    if (info /= 0) then
       print *, 'Error opening output file.'
       stop
    end if
    
    ! Set psb_out_unit to redirect PSBLAS output
    psb_out_unit = unit

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
    tests_number = size(x) * size(y) * size(alpha) * size(beta)
    count = 0

    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    if(my_rank == psb_root_) then
        write(psb_out_unit,'(A,A)') 'Welcome to PSBLAS version: ',psb_version_string_
        write(psb_out_unit,'(A)') 'This is the psb_geaxpby_test sample program'
        write(psb_out_unit,'(A)') ''

        call generate_vectors(arr_size) 
    end if

    call psb_barrier(ctxt)

    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    call psb_geaxpby_kernel(x_file=x(i), y_file=y(j), alpha = real(alpha(k),psb_spk_),& 
                    & beta = real(beta(h),psb_spk_), arr_size = arr_size, ctxt = ctxt, ret = ret)
                    if(my_rank == psb_root_) then
                        count = count + 1
                        if(ret /= -1) then 
                            write(psb_out_unit, '(A,I0,A,I0,A,A)') & 
                            & "Generation geaxpby single precision result file ", count , "/", tests_number, CHAR(9), "[OK]"
                        else
                            write(psb_out_unit, '(A,I0,A,I0,A,A)') & 
                            & "Generation geaxpby single precision result file ", count , "/", tests_number, CHAR(9), "[FAIL]"
                        end if
                    end if
                    call psb_barrier(ctxt)
                end do
            end do
        end do
    end do

    if(my_rank == psb_root_) then
        write(psb_out_unit, *) ''
        count = 0
    end if

    call psb_barrier(ctxt)

    ! Here double precision comparison should be done
    do i=1,size(x)
        do j=1,size(y)
            do k=1,size(alpha)
                do h=1,size(beta)
                    call psb_geaxpby_check(x_file=x(i), y_file=y(j), alpha = alpha(k), beta = beta(h), & 
                    & arr_size = arr_size, ctxt = ctxt, ret = ret)
                    if(my_rank == psb_root_) then
                        count = count + 1
                        if(ret == 0) then 
                            write(psb_out_unit, '(A,I0,A,I0,A,A)') & 
                            & "Double precision check ", count , "/", tests_number, CHAR(9), "[OK]"
                        else
                            write(psb_out_unit, '(A,I0,A,I0,A,A)') & 
                            & "Double precision check ", count , "/", tests_number, CHAR(9), "[FAIL]"
                            goto 9999
                        end if
                    end if
                    call psb_barrier(ctxt)
                end do
            end do
        end do
    end do


    9999 call psb_exit(ctxt)

    return
end program main