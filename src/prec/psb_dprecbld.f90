!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
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
!!$  
subroutine psb_dprecbld(a,desc_a,p,info,upd)

  use psb_serial_mod
  Use psb_spmat_type
  use psb_descriptor_type
  use psb_prec_type
  use psb_tools_mod
  use psb_comm_mod
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  Implicit None

  type(psb_dspmat_type), target              :: a
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_dprec_type),intent(inout)         :: p
  integer, intent(out)                       :: info
  character, intent(in), optional            :: upd

  interface psb_baseprc_bld
    subroutine psb_dbaseprc_bld(a,desc_a,p,info,upd)
      Use psb_spmat_type
      use psb_descriptor_type
      use psb_prec_type
      type(psb_dspmat_type), target              :: a
      type(psb_desc_type), intent(in)            :: desc_a
      type(psb_dbaseprc_type),intent(inout)      :: p
      integer, intent(out)                       :: info
      character, intent(in), optional            :: upd
    end subroutine psb_dbaseprc_bld
  end interface

  interface psb_mlprc_bld
    subroutine psb_dmlprc_bld(a,desc_a,p,info)
      use psb_serial_mod
      use psb_descriptor_type
      use psb_prec_type
      use psb_const_mod
      implicit none 

      type(psb_dspmat_type), intent(in)      :: a
      type(psb_desc_type), intent(in)        :: desc_a
      type(psb_dbaseprc_type), intent(inout) :: p
      integer, intent(out)                   :: info
    end subroutine psb_dmlprc_bld
  end interface

  ! Local scalars
  Integer      :: err, nnzero, I,j,k,icontxt,&
       & me,mycol,nprow,npcol,lw, mtype, nrg, nzg, err_act
  real(kind(1.d0))         :: temp, real_err(5)
  integer      :: int_err(5)
  character    :: iupd

  logical, parameter :: debug=.false. 
  integer,parameter  :: iroot=0,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_precbld'

  if (debug) write(0,*) 'Entering precbld',P%prec,desc_a%matrix_data(:)
  info = 0
  int_err(1) = 0
  icontxt = desc_a%matrix_data(psb_ctxt_)

  if (debug) write(0,*) 'Preconditioner Blacs_gridinfo'
  call blacs_gridinfo(icontxt, nprow, npcol, me, mycol)

  if (present(upd)) then 
    if (debug) write(0,*) 'UPD ', upd
    if ((upd.eq.'F').or.(upd.eq.'T')) then
      iupd=upd
    else
      iupd='F'
    endif
  else
    iupd='F'
  endif

  if (.not.associated(p%baseprecv)) then 
    !! Error 1: should call precset
      info=4010
      ch_err='unassociated bpv'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
  end if
  !
  ! Should add check to ensure all procs have the same... 
  !
  ! ALso should define symbolic names for the preconditioners. 
  !

  if (size(p%baseprecv) >= 1) then 
    call init_baseprc_av(p%baseprecv(1),info)
    if (info /= 0) then 
      info=4010
      ch_err='allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif
    
    call psb_baseprc_bld(a,desc_a,p%baseprecv(1),info,iupd)

  else
      info=4010
      ch_err='size bpv'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999

  endif

  if (size(p%baseprecv) > 1) then
    call init_baseprc_av(p%baseprecv(2),info)
    if (info /= 0) then 
      info=4010
      ch_err='allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif


    call psb_mlprc_bld(a,desc_a,p%baseprecv(2),info)

    if(info /= 0) then
      info=4010
      ch_err='psb_mlprc_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    !
    ! Note: this here has not been tried out. We probably need 
    ! a different baseprc field %desc_ac, in case we try RAS on lev. 2 of 
    ! a 3-level  prec. 
    !
    do i=3, size(p%baseprecv)
      call init_baseprc_av(p%baseprecv(i),info)
      if (info /= 0) then 
        info=4010
        ch_err='allocate'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif

      call psb_mlprc_bld(p%baseprecv(i-1)%av(ac_),p%baseprecv(i-1)%desc_data,&
           & p%baseprecv(i),info)
    end do
    
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

contains

  subroutine init_baseprc_av(p,info)
    type(psb_dbaseprc_type), intent(inout) :: p
    integer                                :: info
    if (associated(p%av)) then 
      ! Have not decided what to do yet
    end if
    allocate(p%av(max_avsz),stat=info)
    if (info /= 0) return
    do k=1,size(p%av)
      call psb_nullify_sp(p%av(k))
    end do
  end subroutine init_baseprc_av

end subroutine psb_dprecbld

