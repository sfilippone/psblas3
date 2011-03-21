!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define PREC_DATA,           !!
!!      structure for preconditioning.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module psb_prec_const_mod

  use psb_base_mod, only : psb_dpk_, psb_spk_, psb_long_int_k_,&
       & psb_err_unit, psb_inp_unit, psb_out_unit

  integer, parameter :: psb_min_prec_=0, psb_noprec_=0, psb_diag_=1, &
       & psb_bjac_=2, psb_max_prec_=2

  ! Entries in iprcparm: preconditioner type, factorization type,
  ! prolongation type, restriction type, renumbering algorithm,
  ! number of overlap layers, pointer to SuperLU factors, 
  ! levels of fill in for ILU(N), 
  integer, parameter :: psb_p_type_=1, psb_f_type_=2
  integer, parameter :: psb_ilu_fill_in_=8
  !Renumbering. SEE BELOW
  integer, parameter :: psb_renum_none_=0, psb_renum_glb_=1, psb_renum_gps_=2
  integer, parameter :: psb_ifpsz=10
  ! Entries in rprcparm: ILU(E) epsilon, smoother omega
  integer, parameter :: psb_fact_eps_=1
  integer, parameter :: psb_rfpsz=4
  ! Factorization types: none, ILU(N), ILU(E)
  integer, parameter :: psb_f_none_=0,psb_f_ilu_n_=1
  ! Fields for sparse matrices ensembles: 
  integer, parameter :: psb_l_pr_=1, psb_u_pr_=2, psb_bp_ilu_avsz=2
  integer, parameter :: psb_max_avsz=psb_bp_ilu_avsz


  interface psb_check_def
    module procedure psb_icheck_def, psb_scheck_def, psb_dcheck_def
  end interface


contains


  function pr_to_str(iprec)

    integer, intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(psb_noprec_)
      pr_to_str='NOPREC'
    case(psb_diag_)         
      pr_to_str='DIAG'
    case(psb_bjac_)         
      pr_to_str='BJAC'
    case default
      pr_to_str='???'
    end select

  end function pr_to_str


  function is_legal_prec(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_prec

    is_legal_prec = ((ip>=psb_noprec_).and.(ip<=psb_bjac_))
    return
  end function is_legal_prec
  function is_legal_ml_fact(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact

    is_legal_ml_fact = (ip == psb_f_ilu_n_)
    return
  end function is_legal_ml_fact
  function is_legal_ml_eps(ip)
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_ml_eps

    is_legal_ml_eps = (ip>=0.0d0)
    return
  end function is_legal_ml_eps


  subroutine psb_icheck_def(ip,name,id,is_legal)
    integer, intent(inout) :: ip
    integer, intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        integer, intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(psb_err_unit,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine psb_icheck_def

  subroutine psb_scheck_def(ip,name,id,is_legal)
    real(psb_spk_), intent(inout) :: ip
    real(psb_spk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        use psb_const_mod
        real(psb_spk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(psb_err_unit,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine psb_scheck_def

  subroutine psb_dcheck_def(ip,name,id,is_legal)
    real(psb_dpk_), intent(inout) :: ip
    real(psb_dpk_), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        use psb_const_mod
        real(psb_dpk_), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(psb_err_unit,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine psb_dcheck_def



end module psb_prec_const_mod
