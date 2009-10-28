!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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

module psb_prec_type

  ! Reduces size of .mod file.
  use psb_base_mod, only : psb_dpk_, psb_spk_, psb_long_int_k_,&
       & psb_desc_type, psb_sizeof, psb_free, psb_cdfree,&
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus,&
       & psb_s_sparse_mat,  psb_d_sparse_mat,  psb_c_sparse_mat,  psb_z_sparse_mat

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


  type psb_s_base_prec_type
  contains
    procedure, pass(prec) :: apply     => s_base_apply
    procedure, pass(prec) :: precbld   => s_base_precbld
    procedure, pass(prec) :: s_base_precseti
    procedure, pass(prec) :: s_base_precsetr
    procedure, pass(prec) :: s_base_precsetc
    procedure, pass(prec) :: sizeof  => s_base_sizeof
    generic, public       :: precset => s_base_precseti, s_base_precsetr, s_base_precsetc
    procedure, pass(prec) :: precinit  => s_base_precinit
    procedure, pass(prec) :: precfree  => s_base_precfree
    procedure, pass(prec) :: precdescr => s_base_precdescr
  end type psb_s_base_prec_type
  
  type psb_sprec_type
    class(psb_s_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: s_apply2v
    procedure, pass(prec)               :: s_apply1v
    generic, public                     :: apply => s_apply2v, s_apply1v
  end type psb_sprec_type

  type psb_d_base_prec_type
  contains
    procedure, pass(prec) :: apply => d_base_apply
    procedure, pass(prec) :: precbld   => d_base_precbld
    procedure, pass(prec) :: d_base_precseti
    procedure, pass(prec) :: d_base_precsetr
    procedure, pass(prec) :: d_base_precsetc
    procedure, pass(prec) :: sizeof  => d_base_sizeof
    generic, public       :: precset => d_base_precseti, d_base_precsetr, d_base_precsetc
    procedure, pass(prec) :: precinit  => d_base_precinit
    procedure, pass(prec) :: precfree  => d_base_precfree
    procedure, pass(prec) :: precdescr => d_base_precdescr
  end type psb_d_base_prec_type
  
  type psb_dprec_type
    class(psb_d_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: d_apply2v
    procedure, pass(prec)               :: d_apply1v
    generic, public                     :: apply => d_apply2v, d_apply1v
  end type psb_dprec_type


  type psb_c_base_prec_type
  contains
    procedure, pass(prec) :: apply     => c_base_apply
    procedure, pass(prec) :: precbld   => c_base_precbld
    procedure, pass(prec) :: c_base_precseti
    procedure, pass(prec) :: c_base_precsetr
    procedure, pass(prec) :: c_base_precsetc
    procedure, pass(prec) :: sizeof  => c_base_sizeof
    generic, public       :: precset => c_base_precseti, c_base_precsetr, c_base_precsetc
    procedure, pass(prec) :: precinit  => c_base_precinit
    procedure, pass(prec) :: precfree  => c_base_precfree
    procedure, pass(prec) :: precdescr => c_base_precdescr
  end type psb_c_base_prec_type
  
  type psb_cprec_type
    class(psb_c_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: c_apply2v
    procedure, pass(prec)               :: c_apply1v
    generic, public                     :: apply => c_apply2v, c_apply1v
  end type psb_cprec_type

  type psb_z_base_prec_type
  contains
    procedure, pass(prec) :: apply => z_base_apply
    procedure, pass(prec) :: precbld   => z_base_precbld
    procedure, pass(prec) :: z_base_precseti
    procedure, pass(prec) :: z_base_precsetr
    procedure, pass(prec) :: z_base_precsetc
    procedure, pass(prec) :: sizeof  => z_base_sizeof
    generic, public       :: precset => z_base_precseti, z_base_precsetr, z_base_precsetc
    procedure, pass(prec) :: precinit  => z_base_precinit
    procedure, pass(prec) :: precfree  => z_base_precfree
    procedure, pass(prec) :: precdescr => z_base_precdescr
  end type psb_z_base_prec_type
  
  type psb_zprec_type
    class(psb_z_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: z_apply2v
    procedure, pass(prec)               :: z_apply1v
    generic, public                     :: apply => z_apply2v, z_apply1v
  end type psb_zprec_type


  character(len=15), parameter, private :: &
       &  fact_names(0:2)=(/'None          ','ILU(n)        ',&
       &  'ILU(eps)      '/)

  interface psb_precfree
    module procedure psb_s_precfree, psb_d_precfree,&
         & psb_c_precfree, psb_z_precfree
  end interface

  interface psb_nullify_prec
    module procedure psb_nullify_sprec, psb_nullify_dprec,&
         & psb_nullify_cprec, psb_nullify_zprec
  end interface

  interface psb_check_def
    module procedure psb_icheck_def, psb_scheck_def, psb_dcheck_def
  end interface

  interface psb_precdescr
    module procedure psb_file_prec_descr, &
         &  psb_sfile_prec_descr, &
         &  psb_cfile_prec_descr, &
         &  psb_zfile_prec_descr
  end interface

  interface psb_sizeof
    module procedure psb_sprec_sizeof, &
         & psb_dprec_sizeof,&
         & psb_cprec_sizeof, &
         & psb_zprec_sizeof
  end interface



  interface psb_precaply
    subroutine psb_sprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only  : psb_desc_type, psb_spk_
      import psb_sprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_sprec_type), intent(in)  :: prec
      real(psb_spk_),intent(in)         :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_sprc_aply
    subroutine psb_sprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod, only  : psb_desc_type, psb_spk_
      import psb_sprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_sprec_type), intent(in)  :: prec
      real(psb_spk_),intent(inout)      :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_sprc_aply1
    subroutine psb_dprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      import psb_dprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(psb_dpk_),intent(in)         :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_dprc_aply
    subroutine psb_dprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      import psb_dprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_dprec_type), intent(in)  :: prec
      real(psb_dpk_),intent(inout)      :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_dprc_aply1
    subroutine psb_cprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only  : psb_desc_type, psb_spk_
      import psb_cprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_cprec_type), intent(in)  :: prec
      complex(psb_spk_),intent(in)      :: x(:)
      complex(psb_spk_),intent(inout)   :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      complex(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_cprc_aply
    subroutine psb_cprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod, only  : psb_desc_type, psb_spk_
      import psb_cprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_cprec_type), intent(in)  :: prec
      complex(psb_spk_),intent(inout)   :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_cprc_aply1
    subroutine psb_zprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      import psb_zprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(psb_dpk_),intent(in)      :: x(:)
      complex(psb_dpk_),intent(inout)   :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_zprc_aply
    subroutine psb_zprc_aply1(prec,x,desc_data,info,trans)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      import psb_zprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_zprec_type), intent(in)  :: prec
      complex(psb_dpk_),intent(inout)   :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_zprc_aply1
  end interface


contains

  


  subroutine psb_file_prec_descr(p,iout)
    use psb_base_mod
    type(psb_dprec_type), intent(in) :: p
    integer, intent(in), optional    :: iout
    integer :: iout_, info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)

  end subroutine psb_file_prec_descr

  subroutine psb_sfile_prec_descr(p,iout)
    use psb_base_mod
    type(psb_sprec_type), intent(in) :: p
    integer, intent(in), optional    :: iout
    integer :: iout_,info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)
    
  end subroutine psb_sfile_prec_descr


  subroutine psb_cfile_prec_descr(p,iout)
    use psb_base_mod
    type(psb_cprec_type), intent(in) :: p
    integer, intent(in), optional    :: iout
    integer :: iout_,info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)
    
  end subroutine psb_cfile_prec_descr


  subroutine psb_zfile_prec_descr(p,iout)
    use psb_base_mod
    type(psb_zprec_type), intent(in) :: p
    integer, intent(in), optional    :: iout
    integer :: iout_,info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)
    
  end subroutine psb_zfile_prec_descr


  function is_legal_prec(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_prec

    is_legal_prec = ((ip>=psb_noprec_).and.(ip<=psb_bjac_))
    return
  end function is_legal_prec
  function is_legal_ml_fact(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact

    is_legal_ml_fact = (ip==psb_f_ilu_n_)
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
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
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
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
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
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine psb_dcheck_def

  subroutine psb_s_precfree(p,info)
    use psb_base_mod
    type(psb_sprec_type), intent(inout) :: p 
    integer, intent(out) ::  info
    integer :: me, err_act,i 
    character(len=20) :: name
    if(psb_get_errstatus() /= 0) return
    info=0
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(p%prec)) then 
      call p%prec%precfree(info)
      if (info /= 0) goto 9999
      deallocate(p%prec,stat=info)
      if (info /= 0) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_s_precfree

  subroutine psb_nullify_sprec(p)
    type(psb_sprec_type), intent(inout) :: p

  end subroutine psb_nullify_sprec

  subroutine psb_d_precfree(p,info)
    use psb_base_mod
    type(psb_dprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer             :: me, err_act,i
    character(len=20)   :: name
    if(psb_get_errstatus() /= 0) return 
    info=0
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(p%prec)) then 
      call p%prec%precfree(info)
      if (info /= 0) goto 9999
      deallocate(p%prec,stat=info)
      if (info /= 0) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_d_precfree

  subroutine psb_nullify_dprec(p)
    type(psb_dprec_type), intent(inout) :: p


  end subroutine psb_nullify_dprec

  subroutine psb_c_precfree(p,info)
    use psb_base_mod
    type(psb_cprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer             :: err_act,i
    character(len=20)   :: name
    if(psb_get_errstatus() /= 0) return 
    info=0
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(p%prec)) then 
      call p%prec%precfree(info)
      if (info /= 0) goto 9999
      deallocate(p%prec,stat=info)
      if (info /= 0) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_precfree

  subroutine psb_nullify_cprec(p)
    type(psb_cprec_type), intent(inout) :: p


  end subroutine psb_nullify_cprec

  subroutine psb_z_precfree(p,info)
    use psb_base_mod
    type(psb_zprec_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer             :: err_act,i
    character(len=20)   :: name
    if(psb_get_errstatus() /= 0) return 
    info=0
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(p%prec)) then 
      call p%prec%precfree(info)
      if (info /= 0) goto 9999
      deallocate(p%prec,stat=info)
      if (info /= 0) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_z_precfree

  subroutine psb_nullify_zprec(p)
    type(psb_zprec_type), intent(inout) :: p

  end subroutine psb_nullify_zprec


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


  function psb_dprec_sizeof(prec) result(val)
    use psb_base_mod
    type(psb_dprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer :: i
    val = 0
    
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
  end function psb_dprec_sizeof

  function psb_sprec_sizeof(prec) result(val)
    use psb_base_mod
    type(psb_sprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
    
  end function psb_sprec_sizeof

  function psb_zprec_sizeof(prec) result(val)
    use psb_base_mod
    type(psb_zprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
    
  end function psb_zprec_sizeof

  function psb_cprec_sizeof(prec) result(val)
    use psb_base_mod
    type(psb_cprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
    
  end function psb_cprec_sizeof
    
 
  subroutine s_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_sprec_type), intent(in)  :: prec
    real(psb_spk_),intent(in)       :: x(:)
    real(psb_spk_),intent(inout)    :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    
    character     :: trans_ 
    real(psb_spk_), pointer :: work_(:)
    integer :: ictxt,np,me,err_act
    character(len=20)   :: name
    
    name='s_apply2v'
    info = 0
    call psb_erractionsave(err_act)
    
    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    if (present(trans)) then 
      trans_=trans
    else
      trans_='N'
    end if
    
    if (present(work)) then 
      work_ => work
    else
      allocate(work_(4*psb_cd_get_local_cols(desc_data)),stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999      
      end if
      
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    call prec%prec%apply(sone,x,szero,y,desc_data,info,trans_,work=work_)
    if (present(work)) then 
    else
      deallocate(work_,stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='DeAllocate')
        goto 9999      
      end if
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_apply2v

  subroutine s_apply1v(prec,x,desc_data,info,trans)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_sprec_type), intent(in)  :: prec
    real(psb_spk_),intent(inout)    :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans

    character     :: trans_
    integer :: ictxt,np,me, err_act
    real(psb_spk_), pointer :: WW(:), w1(:)
    character(len=20)   :: name
    name='s_apply1v'
    info = 0
    call psb_erractionsave(err_act)
    
    
    ictxt=psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    if (present(trans)) then 
      trans_=psb_toupper(trans)
    else
      trans_='N'
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    allocate(ww(size(x)),w1(size(x)),stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if
    call prec%prec%apply(sone,x,szero,ww,desc_data,info,trans_,work=w1)
    if(info /=0) goto 9999
    x(:) = ww(:)
    deallocate(ww,W1,stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine s_apply1v


  subroutine d_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_dprec_type), intent(in)  :: prec
    real(psb_dpk_),intent(in)       :: x(:)
    real(psb_dpk_),intent(inout)    :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    
    character     :: trans_ 
    real(psb_dpk_), pointer :: work_(:)
    integer :: ictxt,np,me,err_act
    character(len=20)   :: name
    
    name='d_apply2v'
    info = 0
    call psb_erractionsave(err_act)
    
    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    if (present(trans)) then 
      trans_=trans
    else
      trans_='N'
    end if
    
    if (present(work)) then 
      work_ => work
    else
      allocate(work_(4*psb_cd_get_local_cols(desc_data)),stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999      
      end if
      
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    call prec%prec%apply(done,x,dzero,y,desc_data,info,trans_,work=work_)
    if (present(work)) then 
    else
      deallocate(work_,stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='DeAllocate')
        goto 9999      
      end if
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_apply2v

  subroutine d_apply1v(prec,x,desc_data,info,trans)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_dprec_type), intent(in)  :: prec
    real(psb_dpk_),intent(inout)    :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans

    character     :: trans_
    integer :: ictxt,np,me, err_act
    real(psb_dpk_), pointer :: WW(:), w1(:)
    character(len=20)   :: name
    name='d_apply1v'
    info = 0
    call psb_erractionsave(err_act)
    
    
    ictxt=psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    if (present(trans)) then 
      trans_=psb_toupper(trans)
    else
      trans_='N'
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    allocate(ww(size(x)),w1(size(x)),stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if
    call prec%prec%apply(done,x,dzero,ww,desc_data,info,trans_,work=w1)
    if(info /=0) goto 9999
    x(:) = ww(:)
    deallocate(ww,W1,stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_apply1v

 
  subroutine c_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_cprec_type), intent(in) :: prec
    complex(psb_spk_),intent(in)      :: x(:)
    complex(psb_spk_),intent(inout)   :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    
    character     :: trans_ 
    complex(psb_spk_), pointer :: work_(:)
    integer :: ictxt,np,me,err_act
    character(len=20)   :: name
    
    name='c_apply2v'
    info = 0
    call psb_erractionsave(err_act)
    
    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    if (present(trans)) then 
      trans_=trans
    else
      trans_='N'
    end if
    
    if (present(work)) then 
      work_ => work
    else
      allocate(work_(4*psb_cd_get_local_cols(desc_data)),stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999      
      end if
      
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    call prec%prec%apply(cone,x,czero,y,desc_data,info,trans_,work=work_)
    if (present(work)) then 
    else
      deallocate(work_,stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='DeAllocate')
        goto 9999      
      end if
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine c_apply2v

  subroutine c_apply1v(prec,x,desc_data,info,trans)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_cprec_type), intent(in) :: prec
    complex(psb_spk_),intent(inout)   :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans

    character     :: trans_
    integer :: ictxt,np,me, err_act
    complex(psb_spk_), pointer :: WW(:), w1(:)
    character(len=20)   :: name
    name='c_apply1v'
    info = 0
    call psb_erractionsave(err_act)
    
    
    ictxt=psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    if (present(trans)) then 
      trans_=psb_toupper(trans)
    else
      trans_='N'
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    allocate(ww(size(x)),w1(size(x)),stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if
    call prec%prec%apply(cone,x,czero,ww,desc_data,info,trans_,work=w1)
    if(info /=0) goto 9999
    x(:) = ww(:)
    deallocate(ww,W1,stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_apply1v


  subroutine z_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_zprec_type), intent(in) :: prec
    complex(psb_dpk_),intent(in)      :: x(:)
    complex(psb_dpk_),intent(inout)   :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    
    character     :: trans_ 
    complex(psb_dpk_), pointer :: work_(:)
    integer :: ictxt,np,me,err_act
    character(len=20)   :: name
    
    name='z_apply2v'
    info = 0
    call psb_erractionsave(err_act)
    
    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    if (present(trans)) then 
      trans_=trans
    else
      trans_='N'
    end if
    
    if (present(work)) then 
      work_ => work
    else
      allocate(work_(4*psb_cd_get_local_cols(desc_data)),stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999      
      end if
      
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    call prec%prec%apply(zone,x,zzero,y,desc_data,info,trans_,work=work_)
    if (present(work)) then 
    else
      deallocate(work_,stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='DeAllocate')
        goto 9999      
      end if
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine z_apply2v

  subroutine z_apply1v(prec,x,desc_data,info,trans)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_zprec_type), intent(in) :: prec
    complex(psb_dpk_),intent(inout)   :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans

    character     :: trans_
    integer :: ictxt,np,me, err_act
    complex(psb_dpk_), pointer :: WW(:), w1(:)
    character(len=20)   :: name
    name='z_apply1v'
    info = 0
    call psb_erractionsave(err_act)
    
    
    ictxt=psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    if (present(trans)) then 
      trans_=psb_toupper(trans)
    else
      trans_='N'
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    allocate(ww(size(x)),w1(size(x)),stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if
    call prec%prec%apply(zone,x,zzero,ww,desc_data,info,trans_,work=w1)
    if(info /=0) goto 9999
    x(:) = ww(:)
    deallocate(ww,W1,stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine z_apply1v


  subroutine s_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_s_base_prec_type), intent(in)  :: prec
    real(psb_spk_),intent(in)         :: alpha, beta
    real(psb_spk_),intent(in)         :: x(:)
    real(psb_spk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_base_apply

  subroutine s_base_precinit(prec,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precinit'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_base_precinit

  subroutine s_base_precbld(a,desc_a,prec,info,upd)
    
    use psb_base_mod
    Implicit None
    
    type(psb_s_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precbld'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_base_precbld

  subroutine s_base_precseti(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precseti'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_base_precseti

  subroutine s_base_precsetr(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precsetr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_base_precsetr

  subroutine s_base_precsetc(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precsetc'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine s_base_precsetc

  subroutine s_base_precfree(prec,info)
    
    use psb_base_mod
    Implicit None

    class(psb_s_base_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precfree'
    
    call psb_erractionsave(err_act)
    
    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine s_base_precfree
  

  subroutine s_base_precdescr(prec,iout)
    
    use psb_base_mod
    Implicit None

    class(psb_s_base_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='s_base_precdescr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine s_base_precdescr
  

  function s_base_sizeof(prec) result(val)
    use psb_base_mod
    class(psb_s_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function s_base_sizeof


  subroutine d_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_d_base_prec_type), intent(in)  :: prec
    real(psb_dpk_),intent(in)         :: alpha, beta
    real(psb_dpk_),intent(in)         :: x(:)
    real(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine d_base_apply

  subroutine d_base_precinit(prec,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_d_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precinit'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_precinit

  subroutine d_base_precbld(a,desc_a,prec,info,upd)
    
    use psb_base_mod
    Implicit None
    
    type(psb_d_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_d_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precbld'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_precbld

  subroutine d_base_precseti(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_d_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precseti'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_precseti

  subroutine d_base_precsetr(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_d_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precsetr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_precsetr

  subroutine d_base_precsetc(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_d_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precsetc'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine d_base_precsetc

  subroutine d_base_precfree(prec,info)
    
    use psb_base_mod
    Implicit None

    class(psb_d_base_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='d_base_precfree'
    
    call psb_erractionsave(err_act)
    
    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_base_precfree
  

  subroutine d_base_precdescr(prec,iout)
    
    use psb_base_mod
    Implicit None

    class(psb_d_base_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='d_base_precdescr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine d_base_precdescr
  

  function d_base_sizeof(prec) result(val)
    use psb_base_mod
    class(psb_d_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function d_base_sizeof


  subroutine c_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)       :: desc_data
    class(psb_c_base_prec_type), intent(in)  :: prec
    complex(psb_spk_),intent(in)         :: alpha, beta
    complex(psb_spk_),intent(in)         :: x(:)
    complex(psb_spk_),intent(inout)      :: y(:)
    integer, intent(out)                 :: info
    character(len=1), optional           :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine c_base_apply

  subroutine c_base_precinit(prec,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precinit'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_precinit

  subroutine c_base_precbld(a,desc_a,prec,info,upd)
    
    use psb_base_mod
    Implicit None
    
    type(psb_c_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precbld'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_precbld

  subroutine c_base_precseti(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precseti'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_precseti

  subroutine c_base_precsetr(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precsetr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_precsetr

  subroutine c_base_precsetc(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precsetc'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine c_base_precsetc

  subroutine c_base_precfree(prec,info)
    
    use psb_base_mod
    Implicit None

    class(psb_c_base_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='c_base_precfree'
    
    call psb_erractionsave(err_act)
    
    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_base_precfree
  

  subroutine c_base_precdescr(prec,iout)
    
    use psb_base_mod
    Implicit None

    class(psb_c_base_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='c_base_precdescr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine c_base_precdescr
  

  function c_base_sizeof(prec) result(val)
    use psb_base_mod
    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function c_base_sizeof


  subroutine z_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_base_mod
    type(psb_desc_type),intent(in)       :: desc_data
    class(psb_z_base_prec_type), intent(in)  :: prec
    complex(psb_dpk_),intent(in)         :: alpha, beta
    complex(psb_dpk_),intent(in)         :: x(:)
    complex(psb_dpk_),intent(inout)      :: y(:)
    integer, intent(out)                 :: info
    character(len=1), optional           :: trans
    complex(psb_dpk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine z_base_apply

  subroutine z_base_precinit(prec,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_z_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precinit'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_base_precinit

  subroutine z_base_precbld(a,desc_a,prec,info,upd)
    
    use psb_base_mod
    Implicit None
    
    type(psb_z_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_z_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precbld'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_base_precbld

  subroutine z_base_precseti(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_z_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precseti'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_base_precseti

  subroutine z_base_precsetr(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_z_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precsetr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_base_precsetr

  subroutine z_base_precsetc(prec,what,val,info)
    
    use psb_base_mod
    Implicit None
    
    class(psb_z_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precsetc'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine z_base_precsetc

  subroutine z_base_precfree(prec,info)
    
    use psb_base_mod
    Implicit None

    class(psb_z_base_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='z_base_precfree'
    
    call psb_erractionsave(err_act)
    
    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine z_base_precfree
  

  subroutine z_base_precdescr(prec,iout)
    
    use psb_base_mod
    Implicit None

    class(psb_z_base_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='z_base_precdescr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine z_base_precdescr
  

  function z_base_sizeof(prec) result(val)
    use psb_base_mod
    class(psb_z_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function z_base_sizeof



end module psb_prec_type
