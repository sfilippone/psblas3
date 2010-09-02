!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!
! File: psb_krylov_mod.f90
!  Interfaces for Krylov subspace iterative methods.
!
Module psb_c_inner_krylov_mod

  use psb_base_inner_krylov_mod

  interface psb_init_conv
    module procedure psb_c_init_conv
  end interface

  interface psb_check_conv
    module procedure psb_c_check_conv
  end interface


contains
  subroutine psb_c_init_conv(methdname,stopc,trace,itmax,a,b,eps,desc_a,stopdat,info)
    use psb_sparse_mod
    implicit none 
    character(len=*), intent(in)      :: methdname
    integer, intent(in)               :: stopc, trace, itmax
    type(psb_c_sparse_mat), intent(in) :: a
    complex(psb_spk_), intent(in)   :: b(:)
    real(psb_spk_), intent(in)      :: eps
    type(psb_desc_type), intent(in)   :: desc_a
    type(psb_itconv_type)             :: stopdat
    integer, intent(out)              :: info

    integer                           :: ictxt, me, np, err_act
    character(len=20)                 :: name

    info = psb_success_
    name = 'psb_init_conv'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)

    stopdat%controls(:) = 0
    stopdat%values(:)   = 0.0d0

    stopdat%controls(psb_ik_stopc_) = stopc
    stopdat%controls(psb_ik_trace_) = trace
    stopdat%controls(psb_ik_itmax_) = itmax

    select case(stopdat%controls(psb_ik_stopc_))
    case (1) 
      stopdat%values(psb_ik_ani_) = psb_spnrmi(a,desc_a,info)
      if (info == psb_success_) stopdat%values(psb_ik_bni_) = psb_geamax(b,desc_a,info)

    case (2) 
      stopdat%values(psb_ik_bn2_) = psb_genrm2(b,desc_a,info)

    case default
      info=psb_err_invalid_istop_
      call psb_errpush(info,name,i_err=(/stopc,0,0,0,0/))
      goto 9999      
    end select
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,a_err="Init conv check data")
      goto 9999
    end if

    stopdat%values(psb_ik_eps_)    = eps
    stopdat%values(psb_ik_errnum_) = dzero
    stopdat%values(psb_ik_errden_) = done

    if ((stopdat%controls(psb_ik_trace_) > 0).and. (me == 0))&
         &  call log_header(methdname) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine psb_c_init_conv


  function psb_c_check_conv(methdname,it,x,r,desc_a,stopdat,info)
    use psb_sparse_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    complex(psb_spk_), intent(in)   :: x(:), r(:)
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    logical                         :: psb_c_check_conv
    integer, intent(out)            :: info

    integer                         :: ictxt, me, np, err_act
    character(len=20)               :: name

    info = psb_success_
    name = 'psb_check_conv'
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)
    psb_c_check_conv = .false. 

    select case(stopdat%controls(psb_ik_stopc_)) 
    case(1)
      stopdat%values(psb_ik_rni_) = psb_geamax(r,desc_a,info)
      if (info == psb_success_) stopdat%values(psb_ik_xni_) = psb_geamax(x,desc_a,info)
      stopdat%values(psb_ik_errnum_) = stopdat%values(psb_ik_rni_)
      stopdat%values(psb_ik_errden_) = &
           & (stopdat%values(psb_ik_ani_)*stopdat%values(psb_ik_xni_)+stopdat%values(psb_ik_bni_))
    case(2)
      stopdat%values(psb_ik_rn2_)  = psb_genrm2(r,desc_a,info)
      stopdat%values(psb_ik_errnum_) = stopdat%values(psb_ik_rn2_)
      stopdat%values(psb_ik_errden_) = stopdat%values(psb_ik_bn2_)

    case default
      info=psb_err_internal_error_
      call psb_errpush(info,name,a_err="Control data in stopdat messed up!")
      goto 9999      
    end select
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    if (stopdat%values(psb_ik_errden_) == dzero) then 
      psb_c_check_conv = (stopdat%values(psb_ik_errnum_) <= stopdat%values(psb_ik_eps_))
    else
      psb_c_check_conv = &
           & (stopdat%values(psb_ik_errnum_) <= stopdat%values(psb_ik_eps_)*stopdat%values(psb_ik_errden_))
    end if

    psb_c_check_conv = (psb_c_check_conv.or.(stopdat%controls(psb_ik_itmax_) <= it))

    if ( (stopdat%controls(psb_ik_trace_) > 0).and.&
         & ((mod(it,stopdat%controls(psb_ik_trace_)) == 0).or.psb_c_check_conv)) then 
      call log_conv(methdname,me,it,1,stopdat%values(psb_ik_errnum_),&
           & stopdat%values(psb_ik_errden_),stopdat%values(psb_ik_eps_))
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end function psb_c_check_conv

end module psb_c_inner_krylov_mod
