!!$ 
!!$ 
!!$              MPcube: Multilevel Parallel Preconditioners Package 
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela Di Serafino    II University of Naples
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
!!$    3. The name of the MPCUBE group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MPCUBE GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  

subroutine psb_dprecset(p,ptype,iv,rs,rv,info)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  implicit none

  type(psb_dprec_type), intent(inout)    :: p
  character(len=*), intent(in)           :: ptype
  integer, optional, intent(in)          :: iv(:)
  real(kind(1.d0)), optional, intent(in) :: rs
  real(kind(1.d0)), optional, intent(in) :: rv(:)
  integer, optional, intent(out)         :: info

  type(psb_dbase_prec), pointer          :: bpv(:)=>null()
  character(len=len(ptype))              :: typeup
  integer                                :: isz, err

  if (present(info)) info = 0

  if (.not.associated(p%baseprecv)) then 
    allocate(p%baseprecv(1),stat=err)
    call psb_nullify_baseprec(p%baseprecv(1))
  endif

  if (.not.associated(p%baseprecv(1)%iprcparm)) then 
    allocate(p%baseprecv(1)%iprcparm(ifpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure',err
    endif
  end if

  select case(toupper(ptype(1:len_trim(ptype))))
  case ('NONE','NOPREC') 
    p%baseprecv(1)%iprcparm(p_type_)     = noprec_
    p%baseprecv(1)%iprcparm(f_type_)     = f_none_
    p%baseprecv(1)%iprcparm(restr_)      = psb_none_
    p%baseprecv(1)%iprcparm(prol_)       = psb_none_
    p%baseprecv(1)%iprcparm(iren_)       = 0
    p%baseprecv(1)%iprcparm(n_ovr_)      = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_) = 1

  case ('DIAG','DIAGSC')
    p%baseprecv(1)%iprcparm(p_type_)     = diagsc_
    p%baseprecv(1)%iprcparm(f_type_)     = f_none_
    p%baseprecv(1)%iprcparm(restr_)      = psb_none_
    p%baseprecv(1)%iprcparm(prol_)       = psb_none_
    p%baseprecv(1)%iprcparm(iren_)       = 0 
    p%baseprecv(1)%iprcparm(n_ovr_)      = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_) = 1

  case ('BJA','ILU') 
    p%baseprecv(1)%iprcparm(p_type_)      = bja_
    p%baseprecv(1)%iprcparm(f_type_)      = f_ilu_n_
    p%baseprecv(1)%iprcparm(restr_)       = psb_none_
    p%baseprecv(1)%iprcparm(prol_)        = psb_none_
    p%baseprecv(1)%iprcparm(iren_)        = 0
    p%baseprecv(1)%iprcparm(n_ovr_)       = 0
    p%baseprecv(1)%iprcparm(ilu_fill_in_) = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_)  = 1

  case ('ASM','AS')
    ! Defaults first 
    p%baseprecv(1)%iprcparm(p_type_)      = asm_
    p%baseprecv(1)%iprcparm(f_type_)      = f_ilu_n_
    p%baseprecv(1)%iprcparm(restr_)       = psb_halo_
    p%baseprecv(1)%iprcparm(prol_)        = psb_none_
    p%baseprecv(1)%iprcparm(iren_)        = 0
    p%baseprecv(1)%iprcparm(n_ovr_)       = 1
    p%baseprecv(1)%iprcparm(ilu_fill_in_) = 0
    p%baseprecv(1)%iprcparm(jac_sweeps_)  = 1
    if (present(iv)) then 
      isz = size(iv) 
      if (isz >= 1) p%baseprecv(1)%iprcparm(n_ovr_)  = iv(1)
      if (isz >= 2) p%baseprecv(1)%iprcparm(restr_)  = iv(2)
      if (isz >= 3) p%baseprecv(1)%iprcparm(prol_)   = iv(3)
      if (isz >= 4) p%baseprecv(1)%iprcparm(f_type_) = iv(4) 
      ! Do not consider renum for the time being. 
!!$      if (isz >= 5) p%baseprecv(1)%iprcparm(iren_) = iv(5)
    end if


  case ('ML', '2LEV')

    select case (size(p%baseprecv)) 
    case(1)
      ! Reallocate
      allocate(bpv(2),stat=err)
      if (err/=0) then 
        write(0,*)'Precset Memory Failure 2l:1',err
      endif
      bpv(1) = p%baseprecv(1)
      call psb_nullify_baseprec(bpv(2))
      deallocate(p%baseprecv)
      p%baseprecv => bpv
      nullify(bpv)

    case(2)
      ! Do nothing

    case default
      ! Error

    end select

    allocate(p%baseprecv(2)%iprcparm(ifpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure 2l:2',err
    endif
    allocate(p%baseprecv(2)%dprcparm(dfpsz),stat=err)
    if (err/=0) then 
      write(0,*)'Precset Memory Failure 2l:3',err
    endif



    p%baseprecv(2)%iprcparm(p_type_)       = bja_
    p%baseprecv(2)%iprcparm(ml_type_)      = mult_ml_prec_
    p%baseprecv(2)%iprcparm(aggr_alg_)     = loc_aggr_
    p%baseprecv(2)%iprcparm(smth_kind_)    = smth_omg_
    p%baseprecv(2)%iprcparm(coarse_mat_)   = mat_distr_
    p%baseprecv(2)%iprcparm(smth_pos_)     = post_smooth_
    p%baseprecv(2)%iprcparm(glb_smth_)     = 1
    p%baseprecv(2)%iprcparm(om_choice_)    = lib_choice_
    p%baseprecv(2)%iprcparm(f_type_)       = f_ilu_n_
    p%baseprecv(2)%iprcparm(ilu_fill_in_)  = 0
    p%baseprecv(2)%dprcparm(smooth_omega_) = 4.d0/3.d0         
    p%baseprecv(2)%iprcparm(jac_sweeps_)   = 1


    if (present(iv)) then 
      isz = size(iv)
      if (isz >= 1) p%baseprecv(2)%iprcparm(ml_type_)      = iv(1)
      if (isz >= 2) p%baseprecv(2)%iprcparm(aggr_alg_)     = iv(2) 
      if (isz >= 3) p%baseprecv(2)%iprcparm(smth_kind_)    = iv(3) 
      if (isz >= 4) p%baseprecv(2)%iprcparm(coarse_mat_)   = iv(4) 
      if (isz >= 5) p%baseprecv(2)%iprcparm(smth_pos_)     = iv(5)
      if (isz >= 6) p%baseprecv(2)%iprcparm(glb_smth_)     = iv(6)
      if (isz >= 7) p%baseprecv(2)%iprcparm(f_type_)       = iv(7)
      if (isz >= 8) p%baseprecv(2)%iprcparm(jac_sweeps_)   = iv(8)

    end if

    if (present(rs)) then 
      p%baseprecv(2)%iprcparm(om_choice_)    = user_choice_
      p%baseprecv(2)%dprcparm(smooth_omega_) = rs      
    end if


  case default
    write(0,*) 'Unknown preconditioner type request "',ptype,'"'
    err = 2

  end select

  if (present(info)) info = err

end subroutine psb_dprecset
