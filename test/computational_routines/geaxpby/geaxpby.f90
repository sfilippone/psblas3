program main
    use psb_geaxpby_test
    use psb_base_mod

    implicit none


    ! MPI variables
    integer(psb_ipk_)       :: my_rank, np

    ! Communicator variable
    type(psb_ctxt_type)     :: ctxt


    call psb_init(ctxt)
    call psb_info(ctxt,my_rank,np)

    if(my_rank == psb_root_) then
        write(psb_out_unit,*) 'Welcome to PSBLAS version: ',psb_version_string_
        write(psb_out_unit,*) 'This is the psb_geaxpby_test sample program'
    end if

    call psb_barrier(ctxt)




    call psb_exit(ctxt)

    return
end program main