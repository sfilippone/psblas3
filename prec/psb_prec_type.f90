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
  use psb_const_mod
  use psb_base_mod, only : psb_sspmat_type, psb_cspmat_type,&
       & psb_dspmat_type, psb_zspmat_type,&
       & psb_desc_type, psb_sizeof

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


  type psb_sprec_type
    type(psb_sspmat_type), allocatable :: av(:) 
    real(psb_spk_), allocatable        :: d(:)  
    type(psb_desc_type)                :: desc_data 
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer                            :: prec, base_prec
  end type psb_sprec_type

  type psb_dprec_type
    type(psb_dspmat_type), allocatable :: av(:) 
    real(psb_dpk_), allocatable        :: d(:)  
    type(psb_desc_type)                :: desc_data 
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer                            :: prec, base_prec
  end type psb_dprec_type

  type psb_cprec_type
    type(psb_cspmat_type), allocatable :: av(:) 
    complex(psb_spk_), allocatable     :: d(:)  
    type(psb_desc_type)                :: desc_data 
    integer, allocatable               :: iprcparm(:) 
    real(psb_spk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer                            :: prec, base_prec
  end type psb_cprec_type


  type psb_zprec_type
    type(psb_zspmat_type), allocatable :: av(:) 
    complex(psb_dpk_), allocatable     :: d(:)  
    type(psb_desc_type)                :: desc_data 
    integer, allocatable               :: iprcparm(:) 
    real(psb_dpk_), allocatable        :: rprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer                            :: prec, base_prec
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

  interface psb_prec_descr
    module procedure psb_out_prec_descr, psb_file_prec_descr, &
         &  psb_sout_prec_descr, psb_sfile_prec_descr, &
         &  psb_cout_prec_descr, psb_cfile_prec_descr, &
         &  psb_zout_prec_descr, psb_zfile_prec_descr
  end interface

  interface psb_sizeof
    module procedure psb_sprec_sizeof,  psb_dprec_sizeof,&
         & psb_cprec_sizeof, psb_zprec_sizeof
  end interface

contains

  subroutine psb_sout_prec_descr(p)
    use psb_base_mod
    type(psb_sprec_type), intent(in) :: p
    call psb_sfile_prec_descr(6,p)
  end subroutine psb_sout_prec_descr

  subroutine psb_out_prec_descr(p)
    use psb_base_mod
    type(psb_dprec_type), intent(in) :: p
    call psb_file_prec_descr(6,p)
  end subroutine psb_out_prec_descr

  subroutine psb_cout_prec_descr(p)
    use psb_base_mod
    type(psb_cprec_type), intent(in) :: p
    call psb_cfile_prec_descr(6,p)
  end subroutine psb_cout_prec_descr

  subroutine psb_zout_prec_descr(p)
    use psb_base_mod
    type(psb_zprec_type), intent(in) :: p
    call psb_zfile_prec_descr(6,p)
  end subroutine psb_zout_prec_descr

  subroutine psb_file_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(psb_dprec_type), intent(in) :: p
    
    write(iout,*) 'Preconditioner description'
    select case(p%iprcparm(psb_p_type_))
    case(psb_noprec_)
      write(iout,*) 'No preconditioning'
    case(psb_diag_)
      write(iout,*) 'Diagonal scaling'
    case(psb_bjac_)
      write(iout,*) 'Block Jacobi with: ',&
           &  fact_names(p%iprcparm(psb_f_type_))
    end select
    
  end subroutine psb_file_prec_descr

  subroutine psb_sfile_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(psb_sprec_type), intent(in) :: p
    
    write(iout,*) 'Preconditioner description'
    select case(p%iprcparm(psb_p_type_))
    case(psb_noprec_)
      write(iout,*) 'No preconditioning'
    case(psb_diag_)
      write(iout,*) 'Diagonal scaling'
    case(psb_bjac_)
      write(iout,*) 'Block Jacobi with: ',&
           &  fact_names(p%iprcparm(psb_f_type_))
    end select
    
  end subroutine psb_sfile_prec_descr

  subroutine psb_cfile_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(psb_cprec_type), intent(in) :: p

    write(iout,*) 'Preconditioner description'
    select case(p%iprcparm(psb_p_type_))
    case(psb_noprec_)
      write(iout,*) 'No preconditioning'
    case(psb_diag_)
      write(iout,*) 'Diagonal scaling'
    case(psb_bjac_)
      write(iout,*) 'Block Jacobi with: ',&
           &  fact_names(p%iprcparm(psb_f_type_))
    end select
  end subroutine psb_cfile_prec_descr

  subroutine psb_zfile_prec_descr(iout,p)
    use psb_base_mod
    integer, intent(in)              :: iout
    type(psb_zprec_type), intent(in) :: p

    write(iout,*) 'Preconditioner description'
    select case(p%iprcparm(psb_p_type_))
    case(psb_noprec_)
      write(iout,*) 'No preconditioning'
    case(psb_diag_)
      write(iout,*) 'Diagonal scaling'
    case(psb_bjac_)
      write(iout,*) 'Block Jacobi with: ',&
           &  fact_names(p%iprcparm(psb_f_type_))
    end select
  end subroutine psb_zfile_prec_descr


  function is_legal_prec(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_prec

    is_legal_prec = ((ip>=psb_noprec_).and.(ip<=psb_bjac_))
    return
  end function is_legal_prec
  function is_legal_ml_fact(ip)
    use psb_base_mod
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact

    is_legal_ml_fact = (ip==psb_f_ilu_n_)
    return
  end function is_legal_ml_fact
  function is_legal_ml_eps(ip)
    use psb_base_mod
    real(psb_dpk_), intent(in) :: ip
    logical             :: is_legal_ml_eps

    is_legal_ml_eps = (ip>=0.0d0)
    return
  end function is_legal_ml_eps


  subroutine psb_icheck_def(ip,name,id,is_legal)
    use psb_base_mod
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
    use psb_base_mod
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
    use psb_base_mod
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

    ! Actually we migh just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    end if

    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)

    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_prec(p)

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
    use psb_base_mod
    type(psb_sprec_type), intent(inout) :: p

!!$    nullify(p%av,p%d,p%iprcparm,p%rprcparm,p%perm,p%invperm,p%mlia,&
!!$         & p%nlaggr,p%base_a,p%base_desc,p%dorig,p%desc_data, p%desc_ac)

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

    ! Actually we migh just deallocate the top level array, except 
    ! for the inner UMFPACK or SLU stuff

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)
    end if

    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)

    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_prec(p)

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
    use psb_base_mod
    type(psb_dprec_type), intent(inout) :: p

!!$    nullify(p%av,p%d,p%iprcparm,p%rprcparm,p%perm,p%invperm,p%mlia,&
!!$         & p%nlaggr,p%base_a,p%base_desc,p%dorig,p%desc_data, p%desc_ac)

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

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)

    end if
    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)

    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_prec(p)
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
    use psb_base_mod
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

    if (allocated(p%d)) then 
      deallocate(p%d,stat=info)
    end if

    if (allocated(p%av))  then 
      do i=1,size(p%av) 
        call psb_sp_free(p%av(i),info)
        if (info /= 0) then 
          ! Actually, we don't care here about this.
          ! Just let it go.
          ! return
        end if
      enddo
      deallocate(p%av,stat=info)

    end if
    if (allocated(p%desc_data%matrix_data)) &
         & call psb_cdfree(p%desc_data,info)

    if (allocated(p%rprcparm)) then 
      deallocate(p%rprcparm,stat=info)
    end if

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_prec(p)
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
    use psb_base_mod
    type(psb_zprec_type), intent(inout) :: p


  end subroutine psb_nullify_zprec


  function pr_to_str(iprec)
    use psb_base_mod

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


  function psb_dprec_sizeof(prec)
    use psb_base_mod
    type(psb_dprec_type), intent(in) :: prec
    integer             :: psb_dprec_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) val = val + psb_sizeof_int * size(prec%iprcparm)
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp  * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + psb_sizeof_dp  * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
    val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
    psb_dprec_sizeof = val 
    
  end function psb_dprec_sizeof

  function psb_sprec_sizeof(prec)
    use psb_base_mod
    type(psb_sprec_type), intent(in) :: prec
    integer             :: psb_sprec_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) val = val + psb_sizeof_int * size(prec%iprcparm)
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp  * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + psb_sizeof_sp  * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
    val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
    psb_sprec_sizeof = val 
    
  end function psb_sprec_sizeof

  function psb_zprec_sizeof(prec)
    use psb_base_mod
    type(psb_zprec_type), intent(in) :: prec
    integer             :: psb_zprec_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) val = val + psb_sizeof_int * size(prec%iprcparm)
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_dp  * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + 2 * psb_sizeof_dp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
    psb_zprec_sizeof = val 
    
  end function psb_zprec_sizeof

  function psb_cprec_sizeof(prec)
    use psb_base_mod
    type(psb_cprec_type), intent(in) :: prec
    integer             :: psb_cprec_sizeof
    integer             :: val,i
    
    val = 0
    if (allocated(prec%iprcparm)) val = val + psb_sizeof_int * size(prec%iprcparm)
    if (allocated(prec%rprcparm)) val = val + psb_sizeof_sp  * size(prec%rprcparm)
    if (allocated(prec%d))        val = val + 2 * psb_sizeof_sp * size(prec%d)
    if (allocated(prec%perm))     val = val + psb_sizeof_int * size(prec%perm)
    if (allocated(prec%invperm))  val = val + psb_sizeof_int * size(prec%invperm)
                                  val = val + psb_sizeof(prec%desc_data)
    if (allocated(prec%av))  then 
      do i=1,size(prec%av)
        val = val + psb_sizeof(prec%av(i))
      end do
    end if
    
    psb_cprec_sizeof = val 
    
  end function psb_cprec_sizeof
    

end module psb_prec_type
