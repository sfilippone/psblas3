!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
module psb_prec_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define PREC_DATA,           !!
!!      structure for preconditioning.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use psb_base_mod

  integer, parameter :: min_prec_=0, noprec_=0, diagsc_=1, bja_=2,&
       & asm_=3, ras_=5, ash_=4, rash_=6, ras2lv_=7, ras2lvm_=8,&
       & lv2mras_=9, lv2smth_=10, lv2lsm_=11, sl2sm_=12, superlu_=13,&
       & new_loc_smth_=14, new_glb_smth_=15, ag2lsm_=16,&
       & msy2l_=18, msy2g_=19, max_prec_=19

  ! Multilevel stuff.
  integer, parameter :: no_ml_=0, add_ml_prec_=1, mult_ml_prec_=2
  integer, parameter :: new_ml_prec_=3, max_ml_=new_ml_prec_
  integer, parameter :: pre_smooth_=1, post_smooth_=2, smooth_both_=3,&
       &    max_smooth_=smooth_both_
  integer, parameter :: loc_aggr_=0, glb_aggr_=1, new_loc_aggr_=2
  integer, parameter :: new_glb_aggr_=3, max_aggr_=new_glb_aggr_
  integer, parameter :: no_smth_=0, smth_omg_=1, smth_biz_=2
  integer, parameter :: lib_choice_=0, user_choice_=1
  integer, parameter :: mat_distr_=0, mat_repl_=1
  ! Entries in iprcparm: preconditioner type, factorization type,
  ! prolongation type, restriction type, renumbering algorithm,
  ! number of overlap layers, pointer to SuperLU factors, 
  ! levels of fill in for ILU(N), 
  integer, parameter :: p_type_=1, f_type_=2, restr_=3, prol_=4
  integer, parameter :: iren_=5, n_ovr_=6
  integer, parameter :: ilu_fill_in_=8, jac_sweeps_=9, ml_type_=10
  integer, parameter :: smth_pos_=11, aggr_alg_=12, smth_kind_=13
  integer, parameter :: om_choice_=14, glb_smth_=15, coarse_mat_=16
  !Renumbering. SEE BELOW
  integer, parameter :: renum_none_=0, renum_glb_=1, renum_gps_=2
  !! 2 ints for 64 bit versions
  integer, parameter :: slu_ptr_=17, umf_symptr_=17, umf_numptr_=19
  integer, parameter :: ifpsz=20
  ! Entries in dprcparm: ILU(E) epsilon, smoother omega
  integer, parameter :: fact_eps_=1, smooth_omega_=2
  integer, parameter :: dfpsz=4
  ! Factorization types: none, ILU(N), ILU(E), SuperLU, UMFPACK
  integer, parameter :: f_none_=0,f_ilu_n_=1,f_ilu_e_=2,f_slu_=3,f_umf_=4
  ! Fields for sparse matrices ensembles: 
  integer, parameter :: l_pr_=1, u_pr_=2, bp_ilu_avsz=2
  integer, parameter :: ap_nd_=3, ac_=4, sm_pr_t_=5, sm_pr_=6
  integer, parameter :: smth_avsz=6, max_avsz=smth_avsz 


  type psb_dbaseprc_type

    type(psb_dspmat_type), allocatable :: av(:) 
    real(kind(1.d0)), allocatable      :: d(:)  
    type(psb_desc_type)                :: desc_data , desc_ac
    integer, allocatable               :: iprcparm(:) 
    real(kind(1.d0)), allocatable      :: dprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_dspmat_type), pointer     :: base_a    => null() !
    type(psb_desc_type), pointer       :: base_desc=> null() ! 
    real(kind(1.d0)), allocatable      :: dorig(:) 

  end type psb_dbaseprc_type


  !
  !  Multilevel preconditioning
  !
  !  To each level I there corresponds a matrix A(I) and a preconditioner K(I)
  !
  !  A notational difference: in the DD reference above the preconditioner for 
  !  a given level K(I) is written out as a sum over the subdomains
  !
  !  SUM_k(R_k^T A_k R_k) 
  !
  !  whereas in this code the sum is implicit in the parallelization, 
  !  i.e. each process takes care of one subdomain, and for each level we have 
  !  as many subdomains as there are processes (except for the coarsest level where 
  !  we might have a replicated index space). Thus the sum apparently disappears 
  !  from our code, but only apparently, because it is implicit in the call 
  !  to psb_baseprc_aply. 
  !
  !  A bit of description of the baseprecv(:) data structure:
  !   1. Number of levels = NLEV = size(baseprecv(:))
  !   2. baseprecv(ilev)%av(:)    sparse matrices needed for the current level. 
  !      Includes:
  !   2.1.:  baseprecv(ilev)%av(l_pr_)    L factor of ILU preconditioners
  !   2.2.:  baseprecv(ilev)%av(u_pr_)    U factor of ILU preconditioners
  !   2.3.:  baseprecv(ilev)%av(ap_nd_)   Off-diagonal part of A for Jacobi sweeps
  !   2.4.:  baseprecv(ilev)%av(ac_)      Aggregated matrix of level ILEV 
  !   2.5.:  baseprecv(ilev)%av(sm_pr_t_) Smoother prolongator transpose; maps vectors  
  !                                          (ilev-1) --->  (ilev) 
  !   2.6.:  baseprecv(ilev)%av(sm_pr_)   Smoother prolongator; maps vectors  
  !                                          (ilev)   --->  (ilev-1) 
  !   Shouldn't we keep just one of them and handle transpose in the sparse BLAS? maybe 
  !
  !   3.    baseprecv(ilev)%desc_data     comm descriptor for level ILEV
  !   4.    baseprecv(ilev)%base_a        Pointer (really a pointer!) to the base matrix 
  !                                       of the current level, i.e.: if ILEV=1 then  A
  !                                       else the aggregated matrix av(ac_); so we have 
  !                                       a unified treatment of residuals. Need this to 
  !                                       avoid passing explicitly matrix A to the 
  !                                       outer prec. routine
  !   5.    baseprecv(ilev)%mlia          The aggregation map from (ilev-1)-->(ilev)
  !                                       if no smoother, it is used instead of sm_pr_
  !   6.    baseprecv(ilev)%nlaggr        Number of aggregates on the various procs. 
  !   
  type psb_dprec_type
    type(psb_dbaseprc_type), allocatable  :: baseprecv(:) 
    ! contain type of preconditioning to be performed
    integer                       :: prec, base_prec
  end type psb_dprec_type

  type psb_zbaseprc_type

    type(psb_zspmat_type), allocatable :: av(:) 
    complex(kind(1.d0)), allocatable   :: d(:)  
    type(psb_desc_type)                :: desc_data , desc_ac
    integer, allocatable               :: iprcparm(:) 
    real(kind(1.d0)), allocatable      :: dprcparm(:) 
    integer, allocatable               :: perm(:),  invperm(:) 
    integer, allocatable               :: mlia(:), nlaggr(:) 
    type(psb_zspmat_type), pointer     :: base_a    => null() !
    type(psb_desc_type), pointer       :: base_desc => null() !
    complex(kind(1.d0)), allocatable   :: dorig(:)

  end type psb_zbaseprc_type

  type psb_zprec_type
    type(psb_zbaseprc_type), allocatable  :: baseprecv(:) 
    ! contain type of preconditioning to be performed
    integer                       :: prec, base_prec
  end type psb_zprec_type


  character(len=15), parameter, private :: &
       &  smooth_names(1:3)=(/'Pre-smoothing ','Post-smoothing',&
       & 'Smooth both   '/)
  character(len=15), parameter, private :: &
       &  smooth_kinds(0:2)=(/'No  smoother  ','Omega smoother',&
       &           'Bizr. smoother'/)
  character(len=15), parameter, private :: &
       &  matrix_names(0:1)=(/'Distributed   ','Replicated    '/)
  character(len=18), parameter, private :: &
       &  aggr_names(0:3)=(/'Local aggregation ','Global aggregation',&
       &     'New local aggr.   ','New global aggr.  '/)
  character(len=6), parameter, private :: &
       &  restrict_names(0:4)=(/'None ','     ','     ','     ','Halo '/)
  character(len=12), parameter, private :: &
       &  prolong_names(0:3)=(/'None       ','Sum        ','Average    ','Square root'/)
  character(len=15), parameter, private :: &
       &  ml_names(0:3)=(/'None          ','Additive      ','Multiplicative',&
       & 'New ML        '/)
  character(len=15), parameter, private :: &
       &  fact_names(0:4)=(/'None          ','ILU(n)        ',&
       &  'ILU(eps)      ','Sparse SuperLU','UMFPACK Sp. LU'/)

  interface psb_base_precfree
    module procedure psb_dbase_precfree, psb_zbase_precfree
  end interface

  interface psb_nullify_baseprec
    module procedure psb_nullify_dbaseprec, psb_nullify_zbaseprec
  end interface

  interface psb_check_def
    module procedure psb_icheck_def, psb_dcheck_def
  end interface

  interface psb_prec_descr
    module procedure psb_out_prec_descr, psb_file_prec_descr, &
         &  psb_zout_prec_descr, psb_zfile_prec_descr
  end interface

  interface psb_prec_short_descr
    module procedure psb_prec_short_descr, psb_zprec_short_descr
  end interface

contains

  subroutine psb_out_prec_descr(p)
    type(psb_dprec_type), intent(in) :: p
    call psb_file_prec_descr(6,p)
  end subroutine psb_out_prec_descr

  subroutine psb_zout_prec_descr(p)
    type(psb_zprec_type), intent(in) :: p
    call psb_zfile_prec_descr(6,p)
  end subroutine psb_zout_prec_descr

  subroutine psb_file_prec_descr(iout,p)
    integer, intent(in)              :: iout
    type(psb_dprec_type), intent(in) :: p
    integer  :: ilev

    write(iout,*) 'Preconditioner description'
    if (allocated(p%baseprecv)) then 
      if (size(p%baseprecv)>=1) then 
        write(iout,*) 'Base preconditioner'
        select case(p%baseprecv(1)%iprcparm(p_type_))
        case(noprec_)
          write(iout,*) 'No preconditioning'
        case(diagsc_)
          write(iout,*) 'Diagonal scaling'
        case(bja_)
          write(iout,*) 'Block Jacobi with: ',&
               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
        case(asm_,ras_,ash_,rash_)
          write(iout,*) 'Additive Schwarz with: ',&
               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
          write(iout,*) 'Overlap:',&
               &  p%baseprecv(1)%iprcparm(n_ovr_)
          write(iout,*) 'Restriction: ',&
               &  restrict_names(p%baseprecv(1)%iprcparm(restr_))
          write(iout,*) 'Prolongation: ',&
               &  prolong_names(p%baseprecv(1)%iprcparm(prol_))
        end select
      end if
      if (size(p%baseprecv)>=2) then 
        do ilev = 2, size(p%baseprecv) 
          if (.not.allocated(p%baseprecv(ilev)%iprcparm)) then 
            write(iout,*) 'Inconsistent MLPREC part!'
            return
          endif
          
          write(iout,*) 'Multilevel: Level No', ilev
          write(iout,*) 'Multilevel type: ',&
               &   ml_names(p%baseprecv(ilev)%iprcparm(ml_type_))
          if (p%baseprecv(ilev)%iprcparm(ml_type_)>no_ml_) then 
            write(iout,*) 'Multilevel aggregation: ', &
                 &   aggr_names(p%baseprecv(ilev)%iprcparm(aggr_alg_))
            write(iout,*) 'Smoother:               ', &
                 &  smooth_kinds(p%baseprecv(ilev)%iprcparm(smth_kind_))
            if (p%baseprecv(ilev)%iprcparm(smth_kind_) /= no_smth_) then 
              write(iout,*) 'Smoothing omega: ', &
                   & p%baseprecv(ilev)%dprcparm(smooth_omega_)
              write(iout,*) 'Smoothing position: ',&
                   & smooth_names(p%baseprecv(ilev)%iprcparm(smth_pos_))
            end if
            write(iout,*) 'Coarse matrix: ',&
                 & matrix_names(p%baseprecv(ilev)%iprcparm(coarse_mat_))
            if (allocated(p%baseprecv(ilev)%nlaggr)) then 
              write(iout,*) 'Aggregation sizes: ', &
                   &  sum( p%baseprecv(ilev)%nlaggr(:)),' : ',p%baseprecv(ilev)%nlaggr(:)
            end if
            write(iout,*) 'Factorization type: ',&
                 & fact_names(p%baseprecv(ilev)%iprcparm(f_type_))
            select case(p%baseprecv(ilev)%iprcparm(f_type_))
            case(f_ilu_n_)      
              write(iout,*) 'Fill level :',p%baseprecv(ilev)%iprcparm(ilu_fill_in_)
            case(f_ilu_e_)         
              write(iout,*) 'Fill threshold :',p%baseprecv(ilev)%dprcparm(fact_eps_)
            case(f_slu_,f_umf_)         
            case default
              write(iout,*) 'Should never get here!'
            end select
            write(iout,*) 'Number of Jacobi sweeps: ', &
                 &   (p%baseprecv(ilev)%iprcparm(jac_sweeps_))
          end if
        end do
      end if

    else
      write(iout,*) 'No Base preconditioner available, something is wrong!'
      return
    endif

  end subroutine psb_file_prec_descr

  function  psb_prec_short_descr(p)
    type(psb_dprec_type), intent(in) :: p
    character(len=20) :: psb_prec_short_descr
    psb_prec_short_descr = ' '
!!$    write(iout,*) 'Preconditioner description'
!!$    if (associated(p%baseprecv)) then 
!!$      if (size(p%baseprecv)>=1) then 
!!$        write(iout,*) 'Base preconditioner'
!!$        select case(p%baseprecv(1)%iprcparm(p_type_))
!!$        case(noprec_)
!!$          write(iout,*) 'No preconditioning'
!!$        case(diagsc_)
!!$          write(iout,*) 'Diagonal scaling'
!!$        case(bja_)
!!$          write(iout,*) 'Block Jacobi with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
!!$        case(asm_,ras_,ash_,rash_)
!!$          write(iout,*) 'Additive Schwarz with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
!!$          write(iout,*) 'Overlap:',&
!!$               &  p%baseprecv(1)%iprcparm(n_ovr_)
!!$          write(iout,*) 'Restriction: ',&
!!$               &  restrict_names(p%baseprecv(1)%iprcparm(restr_))
!!$          write(iout,*) 'Prolongation: ',&
!!$               &  prolong_names(p%baseprecv(1)%iprcparm(prol_))
!!$        end select
!!$      end if
!!$      if (size(p%baseprecv)>=2) then 
!!$        if (.not.associated(p%baseprecv(2)%iprcparm)) then 
!!$          write(iout,*) 'Inconsistent MLPREC part!'
!!$          return
!!$        endif
!!$        write(iout,*) 'Multilevel: ',ml_names(p%baseprecv(2)%iprcparm(ml_type_))
!!$        if (p%baseprecv(2)%iprcparm(ml_type_)>no_ml_) then 
!!$          write(iout,*) 'Multilevel aggregation: ', &
!!$               &   aggr_names(p%baseprecv(2)%iprcparm(aggr_alg_))
!!$          write(iout,*) 'Smoother:               ', &
!!$               &  smooth_kinds(p%baseprecv(2)%iprcparm(smth_kind_))
!!$          write(iout,*) 'Smoothing omega: ', p%baseprecv(2)%dprcparm(smooth_omega_)
!!$          write(iout,*) 'Smoothing position: ',&
!!$               & smooth_names(p%baseprecv(2)%iprcparm(smth_pos_))
!!$          write(iout,*) 'Coarse matrix: ',&
!!$               & matrix_names(p%baseprecv(2)%iprcparm(coarse_mat_))
!!$          write(iout,*) 'Factorization type: ',&
!!$               & fact_names(p%baseprecv(2)%iprcparm(f_type_))
!!$          select case(p%baseprecv(2)%iprcparm(f_type_))
!!$          case(f_ilu_n_)      
!!$            write(iout,*) 'Fill level :',p%baseprecv(2)%iprcparm(ilu_fill_in_)
!!$          case(f_ilu_e_)         
!!$            write(iout,*) 'Fill threshold :',p%baseprecv(2)%dprcparm(fact_eps_)
!!$          case(f_slu_,f_umf_)         
!!$          case default
!!$            write(iout,*) 'Should never get here!'
!!$          end select
!!$          write(iout,*) 'Number of Jacobi sweeps: ', &
!!$               &   (p%baseprecv(2)%iprcparm(jac_sweeps_))
!!$
!!$        end if
!!$      end if
!!$
!!$    else
!!$      write(iout,*) 'No Base preconditioner available, something is wrong!'
!!$      return
!!$    endif

  end function psb_prec_short_descr


  subroutine psb_zfile_prec_descr(iout,p)
    integer, intent(in)              :: iout
    type(psb_zprec_type), intent(in) :: p

    write(iout,*) 'Preconditioner description'
    if (allocated(p%baseprecv)) then 
      if (size(p%baseprecv)>=1) then 
        write(iout,*) 'Base preconditioner'
        select case(p%baseprecv(1)%iprcparm(p_type_))
        case(noprec_)
          write(iout,*) 'No preconditioning'
        case(diagsc_)
          write(iout,*) 'Diagonal scaling'
        case(bja_)
          write(iout,*) 'Block Jacobi with: ',&
               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
        case(asm_,ras_,ash_,rash_)
          write(iout,*) 'Additive Schwarz with: ',&
               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
          write(iout,*) 'Overlap:',&
               &  p%baseprecv(1)%iprcparm(n_ovr_)
          write(iout,*) 'Restriction: ',&
               &  restrict_names(p%baseprecv(1)%iprcparm(restr_))
          write(iout,*) 'Prolongation: ',&
               &  prolong_names(p%baseprecv(1)%iprcparm(prol_))
        end select
      end if
      if (size(p%baseprecv)>=2) then 
        if (.not.allocated(p%baseprecv(2)%iprcparm)) then 
          write(iout,*) 'Inconsistent MLPREC part!'
          return
        endif
        write(iout,*) 'Multilevel: ',ml_names(p%baseprecv(2)%iprcparm(ml_type_))
        if (p%baseprecv(2)%iprcparm(ml_type_)>no_ml_) then 
          write(iout,*) 'Multilevel aggregation: ', &
               &   aggr_names(p%baseprecv(2)%iprcparm(aggr_alg_))
          write(iout,*) 'Smoother:               ', &
               &  smooth_kinds(p%baseprecv(2)%iprcparm(smth_kind_))
          if (p%baseprecv(2)%iprcparm(smth_kind_) /= no_smth_) then 
            write(iout,*) 'Smoothing omega: ', p%baseprecv(2)%dprcparm(smooth_omega_)
            write(iout,*) 'Smoothing position: ',&
                 & smooth_names(p%baseprecv(2)%iprcparm(smth_pos_))
          end if
          
          write(iout,*) 'Coarse matrix: ',&
               & matrix_names(p%baseprecv(2)%iprcparm(coarse_mat_))
          if (allocated(p%baseprecv(ilev)%nlaggr)) then 
            write(iout,*) 'Aggregation sizes: ', &
                 &  sum( p%baseprecv(2)%nlaggr(:)),' : ',p%baseprecv(2)%nlaggr(:)
          endif
          write(iout,*) 'Factorization type: ',&
               & fact_names(p%baseprecv(2)%iprcparm(f_type_))
          select case(p%baseprecv(2)%iprcparm(f_type_))
          case(f_ilu_n_)      
            write(iout,*) 'Fill level :',p%baseprecv(2)%iprcparm(ilu_fill_in_)
          case(f_ilu_e_)         
            write(iout,*) 'Fill threshold :',p%baseprecv(2)%dprcparm(fact_eps_)
          case(f_slu_,f_umf_)         
          case default
            write(iout,*) 'Should never get here!'
          end select
          write(iout,*) 'Number of Jacobi sweeps: ', &
               &   (p%baseprecv(2)%iprcparm(jac_sweeps_))

        end if
      end if

    else
      write(iout,*) 'No Base preconditioner available, something is wrong!'
      return
    endif

  end subroutine psb_zfile_prec_descr

  function  psb_zprec_short_descr(p)
    type(psb_zprec_type), intent(in) :: p
    character(len=20) :: psb_zprec_short_descr
    psb_zprec_short_descr = ' '
!!$    write(iout,*) 'Preconditioner description'
!!$    if (associated(p%baseprecv)) then 
!!$      if (size(p%baseprecv)>=1) then 
!!$        write(iout,*) 'Base preconditioner'
!!$        select case(p%baseprecv(1)%iprcparm(p_type_))
!!$        case(noprec_)
!!$          write(iout,*) 'No preconditioning'
!!$        case(diagsc_)
!!$          write(iout,*) 'Diagonal scaling'
!!$        case(bja_)
!!$          write(iout,*) 'Block Jacobi with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
!!$        case(asm_,ras_,ash_,rash_)
!!$          write(iout,*) 'Additive Schwarz with: ',&
!!$               &  fact_names(p%baseprecv(1)%iprcparm(f_type_))
!!$          write(iout,*) 'Overlap:',&
!!$               &  p%baseprecv(1)%iprcparm(n_ovr_)
!!$          write(iout,*) 'Restriction: ',&
!!$               &  restrict_names(p%baseprecv(1)%iprcparm(restr_))
!!$          write(iout,*) 'Prolongation: ',&
!!$               &  prolong_names(p%baseprecv(1)%iprcparm(prol_))
!!$        end select
!!$      end if
!!$      if (size(p%baseprecv)>=2) then 
!!$        if (.not.associated(p%baseprecv(2)%iprcparm)) then 
!!$          write(iout,*) 'Inconsistent MLPREC part!'
!!$          return
!!$        endif
!!$        write(iout,*) 'Multilevel: ',ml_names(p%baseprecv(2)%iprcparm(ml_type_))
!!$        if (p%baseprecv(2)%iprcparm(ml_type_)>no_ml_) then 
!!$          write(iout,*) 'Multilevel aggregation: ', &
!!$               &   aggr_names(p%baseprecv(2)%iprcparm(aggr_alg_))
!!$          write(iout,*) 'Smoother:               ', &
!!$               &  smooth_kinds(p%baseprecv(2)%iprcparm(smth_kind_))
!!$          write(iout,*) 'Smoothing omega: ', p%baseprecv(2)%dprcparm(smooth_omega_)
!!$          write(iout,*) 'Smoothing position: ',&
!!$               & smooth_names(p%baseprecv(2)%iprcparm(smth_pos_))
!!$          write(iout,*) 'Coarse matrix: ',&
!!$               & matrix_names(p%baseprecv(2)%iprcparm(coarse_mat_))
!!$          write(iout,*) 'Factorization type: ',&
!!$               & fact_names(p%baseprecv(2)%iprcparm(f_type_))
!!$          select case(p%baseprecv(2)%iprcparm(f_type_))
!!$          case(f_ilu_n_)      
!!$            write(iout,*) 'Fill level :',p%baseprecv(2)%iprcparm(ilu_fill_in_)
!!$          case(f_ilu_e_)         
!!$            write(iout,*) 'Fill threshold :',p%baseprecv(2)%dprcparm(fact_eps_)
!!$          case(f_slu_,f_umf_)         
!!$          case default
!!$            write(iout,*) 'Should never get here!'
!!$          end select
!!$          write(iout,*) 'Number of Jacobi sweeps: ', &
!!$               &   (p%baseprecv(2)%iprcparm(jac_sweeps_))
!!$
!!$        end if
!!$      end if
!!$
!!$    else
!!$      write(iout,*) 'No Base preconditioner available, something is wrong!'
!!$      return
!!$    endif

  end function psb_zprec_short_descr




  function is_legal_base_prec(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_base_prec

    is_legal_base_prec = ((ip>=noprec_).and.(ip<=rash_))
    return
  end function is_legal_base_prec
  function is_legal_n_ovr(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_n_ovr

    is_legal_n_ovr = (ip >=0) 
    return
  end function is_legal_n_ovr
  function is_legal_renum(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_renum
    ! For the time being we are disabling renumbering options. 
    is_legal_renum = (ip ==0) 
    return
  end function is_legal_renum
  function is_legal_jac_sweeps(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_jac_sweeps

    is_legal_jac_sweeps = (ip >= 1) 
    return
  end function is_legal_jac_sweeps
  function is_legal_prolong(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_prolong

    is_legal_prolong = ((ip>=psb_none_).and.(ip<=psb_square_root_))
    return
  end function is_legal_prolong
  function is_legal_restrict(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_restrict
    is_legal_restrict = ((ip==psb_nohalo_).or.(ip==psb_halo_))
    return
  end function is_legal_restrict
  function is_legal_ml_type(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_type

    is_legal_ml_type = ((ip>=no_ml_).and.(ip<=max_ml_))
    return
  end function is_legal_ml_type
  function is_legal_ml_aggr_kind(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_aggr_kind

    is_legal_ml_aggr_kind = ((ip>=loc_aggr_).and.(ip<=max_aggr_))
    return
  end function is_legal_ml_aggr_kind
  function is_legal_ml_smooth_pos(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_smooth_pos

    is_legal_ml_smooth_pos = ((ip>=pre_smooth_).and.(ip<=max_smooth_))
    return
  end function is_legal_ml_smooth_pos
  function is_legal_ml_smth_kind(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_smth_kind

    is_legal_ml_smth_kind = ((ip>=no_smth_).and.(ip<=smth_biz_))
    return
  end function is_legal_ml_smth_kind
  function is_legal_ml_coarse_mat(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_coarse_mat

    is_legal_ml_coarse_mat = ((ip>=mat_distr_).and.(ip<=mat_repl_))
    return
  end function is_legal_ml_coarse_mat
  function is_legal_ml_fact(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_fact

    is_legal_ml_fact = ((ip>=f_ilu_n_).and.(ip<=f_umf_))
    return
  end function is_legal_ml_fact
  function is_legal_ml_lev(ip)
    integer, intent(in) :: ip
    logical             :: is_legal_ml_lev

    is_legal_ml_lev = (ip>=0)
    return
  end function is_legal_ml_lev
  function is_legal_omega(ip)
    real(kind(1.d0)), intent(in) :: ip
    logical             :: is_legal_omega

    is_legal_omega = ((ip>=0.0d0).and.(ip<=2.0d0))
    return
  end function is_legal_omega
  function is_legal_ml_eps(ip)
    real(kind(1.d0)), intent(in) :: ip
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

  subroutine psb_dcheck_def(ip,name,id,is_legal)
    real(kind(1.d0)), intent(inout) :: ip
    real(kind(1.d0)), intent(in)    :: id
    character(len=*), intent(in) :: name
    interface 
      function is_legal(i)
        real(kind(1.d0)), intent(in) :: i
        logical             :: is_legal
      end function is_legal
    end interface

    if (.not.is_legal(ip)) then     
      write(0,*) 'Illegal value for ',name,' :',ip, '. defaulting to ',id
      ip = id
    end if
  end subroutine psb_dcheck_def

  subroutine psb_dbase_precfree(p,info)
    use psb_base_mod

    type(psb_dbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

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
    if (allocated(p%desc_ac%matrix_data)) &
         & call psb_cdfree(p%desc_ac,info)
    
    if (allocated(p%dprcparm)) then 
      deallocate(p%dprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%dorig)) then 
      deallocate(p%dorig,stat=info)
    endif

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(f_type_)==f_slu_) then 
        call psb_dslu_free(p%iprcparm(slu_ptr_),info)
      end if
      if (p%iprcparm(f_type_)==f_umf_) then 
        call psb_dumf_free(p%iprcparm(umf_symptr_),&
             & p%iprcparm(umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_baseprec(p)
  end subroutine psb_dbase_precfree

  subroutine psb_nullify_dbaseprec(p)
    use psb_base_mod

    type(psb_dbaseprc_type), intent(inout) :: p

    nullify(p%base_a) 
    nullify(p%base_desc) 
!!$    nullify(p%av,p%d,p%iprcparm,p%dprcparm,p%perm,p%invperm,p%mlia,&
!!$         & p%nlaggr,p%base_a,p%base_desc,p%dorig,p%desc_data, p%desc_ac)

  end subroutine psb_nullify_dbaseprec

  subroutine psb_zbase_precfree(p,info)
    use psb_base_mod
    type(psb_zbaseprc_type), intent(inout) :: p
    integer, intent(out)                :: info
    integer :: i

    info = 0

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
    ! call psb_cdfree(p%desc_data,info)
    ! call psb_cdfree(p%desc_ac,info)

    if (allocated(p%dprcparm)) then 
      deallocate(p%dprcparm,stat=info)
    end if
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_a) 
    ! This is a pointer to something else, must not free it here. 
    nullify(p%base_desc) 

    if (allocated(p%dorig)) then 
      deallocate(p%dorig,stat=info)
    endif

    if (allocated(p%mlia)) then 
      deallocate(p%mlia,stat=info)
    endif

    if (allocated(p%nlaggr)) then 
      deallocate(p%nlaggr,stat=info)
    endif

    if (allocated(p%perm)) then 
      deallocate(p%perm,stat=info)
    endif

    if (allocated(p%invperm)) then 
      deallocate(p%invperm,stat=info)
    endif

    if (allocated(p%iprcparm)) then 
      if (p%iprcparm(f_type_)==f_slu_) then 
        call psb_zslu_free(p%iprcparm(slu_ptr_),info)
      end if
      if (p%iprcparm(f_type_)==f_umf_) then 
        call psb_zumf_free(p%iprcparm(umf_symptr_),&
             & p%iprcparm(umf_numptr_),info)
      end if
      deallocate(p%iprcparm,stat=info)
    end if
    call psb_nullify_baseprec(p)
  end subroutine psb_zbase_precfree

  subroutine psb_nullify_zbaseprec(p)
    use psb_base_mod

    type(psb_zbaseprc_type), intent(inout) :: p


    nullify(p%base_a) 
    nullify(p%base_desc) 

  end subroutine psb_nullify_zbaseprec


  function pr_to_str(iprec)

    integer, intent(in)  :: iprec
    character(len=10)     :: pr_to_str

    select case(iprec)
    case(noprec_)
      pr_to_str='NOPREC'
    case(diagsc_)         
      pr_to_str='DIAGSC'
    case(bja_)         
      pr_to_str='BJA'
    case(asm_)      
      pr_to_str='ASM'
    case(ash_)      
      pr_to_str='ASM'
    case(ras_)      
      pr_to_str='ASM'
    case(rash_)      
      pr_to_str='ASM'
    end select

  end function pr_to_str

end module psb_prec_type
