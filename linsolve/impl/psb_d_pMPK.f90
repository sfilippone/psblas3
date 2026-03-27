subroutine psb_d_pMPK_packd(spmat, prec, vec_in, mvec_out, s, desc, info, base_type, alpha, beta, gamma)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    type(psb_dspmat_type), intent(in)           :: spmat
    class(psb_dprec_type), intent(inout)        :: prec
    type(psb_d_vect_type), intent(inout)        :: vec_in
    type(psb_d_multivect_type), intent(inout)   :: mvec_out
    integer(psb_ipk_), intent(in)               :: s
    type(psb_desc_type), intent(in)             :: desc
    integer(psb_ipk_), intent(out)              :: info
    character, optional, intent(in)             :: base_type
    real(psb_dpk_), optional, intent(in)        :: alpha, beta, gamma

    integer(psb_ipk_) :: err_act
    real(psb_dpk_)    :: gamma_
    logical           :: save_r
    character         :: base_type_
    character(len=20) :: name = "psb_d_pMPK"

    call psb_erractionsave(err_act)
    info = psb_success_

    !Check s value
    if(s <= 0)  then
        info = psb_err_iarg_invalid_value_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Dimension checks. mvec_out = n*(2s + 1) or n*2s
    if((mvec_out%get_ncols() /= 2*s) .and. (mvec_out%get_ncols() /= 2*s + 1)) then
        info = psb_err_invalid_mvect_size_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Select if save r as the first vector of Q using is dimension
    save_r = (mvec_out%get_ncols() == 2*s + 1);

    !Check type base
    if(present(base_type)) then 
        base_type_ = psb_toupper(base_type)
    else
        base_type_ = "M"
    endif

    select case(base_type_)
        case ("M")
            call psb_d_pMPK_packd_monomial() 
        
        case ("C")
            ! Check presence of Chebyshev parameters alfa and beta (error or default values?)
            if((.not. present(alpha)) .or. (.not. present(beta))) then
                info = psb_err_invalid_args_combination_
                call psb_errpush(info, name)
                goto 9999
            endif
            !Chebyshev parameters gamma defaulted to 1 if not present
            if(present(gamma)) then
                gamma_ = gamma
            else
                gamma_ = done
            endif

            call psb_d_pMPK_packd_chebyshev()
        
        case default
            info = psb_err_invalid_input_
            call psb_errpush(info, name)
            goto 9999
    end select

9999 call psb_error_handler(err_act)
    return

contains
    subroutine psb_d_pMPK_packd_monomial()
        implicit none
        character(len=20) :: name = "psb_d_pMPK_monomial"
        integer(psb_ipk_) :: i, idx_Z, idx_Q

        ! Copy r in the first column of Q
        idx_Q = s + 1
        call psb_geaxpby(done, vec_in, dzero, mvec_out, idx_Q, desc, info)

        ! First iteration
        idx_Z = 1
        call prec%apply(mvec_out, idx_Q, mvec_out, idx_Z, desc, info)
        idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
        call psb_spmm(done, spmat, mvec_out, idx_Z, dzero, mvec_out, idx_Q, desc, info)
        
        if(s == 1) return

        do i = 2, s
            idx_Z = idx_Z + 1
            call prec%apply(mvec_out, idx_Q, mvec_out, idx_Z, desc, info)
            
            idx_Q = idx_Q + 1
            call psb_spmm(done, spmat, mvec_out, idx_Z, dzero, mvec_out, idx_Q, desc, info)
        end do
    end subroutine psb_d_pMPK_packd_monomial

    subroutine psb_d_pMPK_packd_chebyshev()
        implicit none
        character(len=20) :: name = "psb_d_pMPK_chebyshev"
        integer(psb_ipk_) :: i, idx_Z, idx_Q, ind_tmp
        type(psb_d_multivect_type) :: mvec_tmp

        call psb_geall(mvec_tmp, desc, info, n = 3)    ! I need only tre stored temp vectors
        call psb_geasb(mvec_tmp, desc, info)
        
        ! Inizialize first column of temp multivector
        ind_tmp = 1
        call psb_geaxpby(done, vec_in, dzero, mvec_tmp, ind_tmp, desc, info)

        ! First iteration
        idx_Z = 1
        call prec%apply(mvec_tmp, ind_tmp, mvec_out, idx_Z, desc, info)
        idx_Q = s + 1
        call psb_spmm(done, spmat, mvec_out, idx_Z, dzero, mvec_out, idx_Q, desc, info)

        !Check first early exit
        if(s == 1) goto 9998

        ! First update of temp multivector + second iteration
        ind_tmp = ind_tmp + 1
        call psb_geaxpby(alpha, mvec_out, idx_Q, -beta, mvec_tmp, ind_tmp - 1, dzero, mvec_tmp, ind_tmp, desc, info)
        idx_Z = idx_Z + 1
        call prec%apply(mvec_tmp, ind_tmp, mvec_out, idx_Z, desc, info)
        idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
        call psb_spmm(done, spmat, mvec_out, idx_Z, dzero, mvec_out, idx_Q, desc, info)

        !Check second early exit
        if(s == 2) goto 9998

        ! Loop for s > 2
        do i = 3, s
            ind_tmp = mod(ind_tmp, 3) + 1
            call psb_geaxpby(2*alpha, mvec_out, idx_Q, -2*beta, mvec_tmp, mod(ind_tmp - 1, 3) + 1, &
                                & -gamma_, mvec_tmp, mod(ind_tmp - 2, 3) + 1, mvec_tmp, ind_tmp, desc, info)

            idx_Z = idx_Z + 1
            call prec%apply(mvec_tmp, ind_tmp, mvec_out, idx_Z, desc, info)

            idx_Q = idx_Q + 1
            call psb_spmm(done, spmat, mvec_out, idx_Z, dzero, mvec_out, idx_Q, desc, info)
        end do

    9998 call psb_gefree(mvec_tmp, desc, info)
        return
    end subroutine psb_d_pMPK_packd_chebyshev
end subroutine psb_d_pMPK_packd    

subroutine psb_d_pMPK_split(spmat, prec, vec_in, Z, Q, s, desc, info, base_type, alpha, beta, gamma)
    use psb_base_mod
    use psb_prec_mod
    implicit none

    type(psb_dspmat_type), intent(in)           :: spmat
    class(psb_dprec_type), intent(inout)        :: prec
    type(psb_d_vect_type), intent(inout)        :: vec_in
    type(psb_d_multivect_type), intent(inout)   :: Z, Q
    integer(psb_ipk_), intent(in)               :: s
    type(psb_desc_type), intent(in)             :: desc
    integer(psb_ipk_), intent(out)              :: info
    character, optional, intent(in)             :: base_type
    real(psb_dpk_), optional, intent(in)        :: alpha, beta, gamma

    integer(psb_ipk_) :: err_act
    real(psb_dpk_)    :: gamma_
    logical           :: save_r
    character         :: base_type_
    character(len=20) :: name = "psb_d_pMPK"

    call psb_erractionsave(err_act)
    info = psb_success_

    !Check s value
    if(s <= 0)  then
        info = psb_err_iarg_invalid_value_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Dimension checks. Z = n*s Q = n*(s+1) or n*s
    if((Z%get_ncols() /= s) .or. ((Q%get_ncols() /= s) .and. (Q%get_ncols() /= s + 1))) then
        info = psb_err_invalid_mvect_size_
        call psb_errpush(info, name)
        goto 9999
    endif

    !Select if save r as the first vector of Q using is dimension
    save_r = (Q%get_ncols() == s + 1);

    !Check type base
    if(present(base_type)) then 
        base_type_ = psb_toupper(base_type)
    else
        base_type_ = "M"
    endif

    select case(base_type_)
        case ("M")
            call psb_d_pMPK_split_monomial() 
        
        case ("C")
            ! Check presence of Chebyshev parameters alfa and beta (error or default values?)
            if((.not. present(alpha)) .or. (.not. present(beta))) then
                info = psb_err_invalid_args_combination_
                call psb_errpush(info, name)
                goto 9999
            endif

            !Chebyshev parameters gamma defaulted to 1 if not present
            if(present(gamma)) then
                gamma_ = gamma
            else
                gamma_ = done
            endif

            call psb_d_pMPK_split_chebyshev()
        
        case default
            info = psb_err_invalid_input_
            call psb_errpush(info, name)
            goto 9999
    end select

9999 call psb_error_handler(err_act)
    return

contains
    subroutine psb_d_pMPK_split_monomial()
        implicit none
        character(len=20) :: name = "psb_d_pMPK_monomial"
        integer(psb_ipk_) :: i, idx_Z, idx_Q

        ! Copy r in the first column of Q
        idx_Q = 1
        call psb_geaxpby(done, vec_in, dzero, Q, idx_Q, desc, info)


        ! First iteration
        idx_Z = 1
        call prec%apply(Q, idx_Q, Z, idx_Z, desc, info)
        idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
        call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info)
        
        if(s == 1) return

        do i = 2, s
            idx_Z = idx_Z + 1
            call prec%apply(Q, idx_Q, Z, idx_Z, desc, info)
            
            idx_Q = idx_Q + 1
            call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info)
        end do
    end subroutine psb_d_pMPK_split_monomial

    subroutine psb_d_pMPK_split_chebyshev()
        implicit none
        character(len=20) :: name = "psb_d_pMPK_chebyshev"
        integer(psb_ipk_) :: i, idx_Z, idx_Q, ind_tmp
        type(psb_d_multivect_type) :: mvec_tmp

        call psb_geall(mvec_tmp, desc, info, n = 3)    ! I need only tre stored temp vectors
        call psb_geasb(mvec_tmp, desc, info)
        
        ! Copy r in the first column of Q
        idx_Q = 1
        call psb_geaxpby(done, vec_in, dzero, Q, idx_Q, desc, info)
        
        ! Inizialize first column of temp multivector
        ind_tmp = 1
        call psb_geaxpby(done, vec_in, dzero, mvec_tmp, ind_tmp, desc, info)

        ! First iteration
        idx_Z = 1 
        call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info)
        idx_Q = merge(idx_Q + 1, idx_Q, save_r) !If not save_r overwrite it with the first vector
        call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info)

        !Check first early exit
        if(s == 1) goto 9998

        ! First update of temp multivector + second iteration
        ind_tmp = ind_tmp + 1
        call psb_geaxpby(alpha, Q, idx_Q, -beta, mvec_tmp, ind_tmp - 1, dzero, mvec_tmp, ind_tmp, desc, info)
        idx_Z = idx_Z + 1
        call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info)
        idx_Q = idx_Q + 1
        call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info)

        !Check second early exit
        if(s == 2) goto 9998

        ! Loop for s > 2
        do i = 3, s
            ind_tmp = mod(ind_tmp, 3) + 1
            call psb_geaxpby(2*alpha, Q, idx_Q, -2*beta, mvec_tmp, mod(ind_tmp - 1, 3) + 1, &
                                & -gamma_, mvec_tmp, mod(ind_tmp - 2, 3) + 1, mvec_tmp, ind_tmp, desc, info)

            idx_Z = idx_Z + 1
            call prec%apply(mvec_tmp, ind_tmp, Z, idx_Z, desc, info)

            idx_Q = idx_Q + 1
            call psb_spmm(done, spmat, Z, idx_Z, dzero, Q, idx_Q, desc, info)
        end do
    
    9998 call psb_gefree(mvec_tmp, desc, info)
        return
    end subroutine psb_d_pMPK_split_chebyshev
end subroutine psb_d_pMPK_split