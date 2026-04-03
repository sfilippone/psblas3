program main
    use psb_spmv_overlap_test
    use psb_base_mod

    implicit none

    integer(psb_ipk_)       :: my_rank, np
    integer(psb_ipk_)       :: k, h
    type(psb_ctxt_type)     :: ctxt



    call psb_init(ctxt)
    call psb_info(ctxt, my_rank, np)

    if (my_rank == psb_root_) then
        write(psb_out_unit,*) 'Welcome to PSBLAS version: ', psb_version_string_
        write(psb_out_unit,*) 'This is the psb_spmv_overlap_test sample program'
    end if

    call psb_barrier(ctxt)


    call psb_spmv_overlap_kernel(ctxt)


    call psb_exit(ctxt)
end program main
