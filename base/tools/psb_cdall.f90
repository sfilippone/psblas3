!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
! File: psb_cdall.f90
!
! Subroutine: psb_cdall
!    Allocate descriptor Outer routine
!
subroutine psb_cdall(ctxt, desc, info,mg,ng,parts,&
     & vg,vl,flag,nl,repl,globalcheck,lidx,usehash)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_cd_tools_mod, psb_protect_name => psb_cdall
  use psi_mod
  implicit None
  procedure(psb_parts)             :: parts
  integer(psb_lpk_), intent(in)    :: mg,ng, vl(:)
  type(psb_ctxt_type), intent(in)  :: ctxt
  integer(psb_ipk_), intent(in)    :: vg(:), lidx(:),nl
  integer(psb_ipk_), intent(in)    :: flag
  logical, intent(in)              :: repl, globalcheck,usehash
  integer(psb_ipk_), intent(out)   :: info
  type(psb_desc_type), intent(out) :: desc

  optional :: mg,ng,parts,vg,vl,flag,nl,repl, globalcheck,lidx, usehash

  interface 
    subroutine psb_cdals(m, n, parts, ctxt, desc, info)
      use psb_desc_mod
      procedure(psb_parts)               :: parts
      integer(psb_lpk_), intent(in)      :: m,n
      type(psb_ctxt_type), intent(in) :: ctxt
      Type(psb_desc_type), intent(out)   :: desc
      integer(psb_ipk_), intent(out)     :: info
    end subroutine psb_cdals
    subroutine psb_cdalv(v, ctxt, desc, info, flag)
      use psb_desc_mod
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_ipk_), intent(in)               :: v(:)
      integer(psb_ipk_), intent(in), optional     :: flag
      integer(psb_ipk_), intent(out)              :: info
      Type(psb_desc_type), intent(out)  :: desc
    end subroutine psb_cdalv
    subroutine psb_cd_inloc(v, ctxt, desc, info, globalcheck,idx, usehash)
      use psb_desc_mod
      implicit None
      type(psb_ctxt_type), intent(in) :: ctxt
      integer(psb_lpk_), intent(in)               :: v(:)
      integer(psb_ipk_), intent(out)              :: info
      type(psb_desc_type), intent(out)  :: desc
      logical, intent(in), optional     :: globalcheck, usehash
      integer(psb_ipk_), intent(in), optional     :: idx(:)
    end subroutine psb_cd_inloc
    subroutine psb_cdrep(m, ctxt, desc,info)
      use psb_desc_mod
      integer(psb_lpk_), intent(in)     :: m
      type(psb_ctxt_type), intent(in) :: ctxt


      Type(psb_desc_type), intent(out)  :: desc
      integer(psb_ipk_), intent(out)    :: info
    end subroutine psb_cdrep
  end interface
  character(len=20)   :: name
  integer(psb_ipk_) :: err_act, flag_, i, me, np, nnv, lr
  integer(psb_lpk_) :: n_, nlp
  logical            :: usehash_
  integer(psb_ipk_), allocatable :: itmpv(:)
  integer(psb_lpk_), allocatable :: lvl(:)
  logical, parameter :: timings=.false.
  real(psb_dpk_)  :: t0, t1


  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  name = 'psb_cdall'
  call psb_erractionsave(err_act)
  if (timings) t0 = psb_wtime()
  call psb_info(ctxt, me, np)
  if (count((/ present(vg),present(vl),&
       &  present(parts),present(nl), present(repl) /)) < 1) then 
    !
    ! Allocate a NULL descriptor, its only role is to
    ! store a CTXT for future reference
    !
    allocate(psb_indx_map      :: desc%indxmap, stat=info)
    call desc%indxmap%init_null(ctxt,info)
    if (info /= 0) goto 9999
    goto 9998
  endif

  desc%base_desc => null() 
  if (allocated(desc%indxmap)) then 
    write(0,*) 'Allocated on an intent(OUT) var?'
  end if

  if (present(parts)) then 

    if (.not.present(mg)) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 9999 
    end if
    if (present(ng)) then 
      n_ = ng
    else
      n_ = mg 
    endif
    call  psb_cdals(mg, n_, parts, ctxt, desc, info)

  else if (present(repl)) then 

    if (.not.present(mg)) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 9999 
    end if
    if (.not.repl) then 
      info=psb_err_no_optional_arg_
      call psb_errpush(info,name)
      goto 9999 
    end if

    call  psb_cdrep(mg, ctxt, desc, info)


  else if (present(vg)) then 

    if (present(flag)) then 
      flag_=flag
    else
      flag_=0
    endif
    if (present(mg)) then 
      nnv = min(mg,size(vg))
    else
      nnv = size(vg)
    end if

    call psb_cdalv(vg(1:nnv), ctxt, desc, info, flag=flag_)

  else if (present(vl)) then 

    if (present(nl)) then 
      nnv = min(nl,size(vl))
    else
      nnv = size(vl)
    end if

    call psb_cd_inloc(vl(1:nnv),ctxt,desc,info, globalcheck=globalcheck,idx=lidx)

  else if (present(nl)) then 

    if (present(usehash)) then
      usehash_ = usehash
    else
      usehash_ = .false.
    end if

    if (usehash_) then
      nlp = nl
      call psb_exscan_sum(ctxt,nlp)
      lvl = [ (i,i=1,nl) ] + nlp
      call psb_cd_inloc(lvl(1:nl),ctxt,desc,info, globalcheck=.false.)

    else
      if (np == 1) then 
        allocate(psb_repl_map      :: desc%indxmap, stat=info)
      else
        allocate(psb_gen_block_map :: desc%indxmap, stat=info)
      end if
      if (info == psb_success_) then 
        select type(aa => desc%indxmap) 
        type is (psb_repl_map)
          n_ = nl
          call aa%repl_map_init(ctxt,n_,info)
        type is (psb_gen_block_map) 
          call aa%gen_block_map_init(ctxt,nl,info)
          class default 
            ! This cannot happen 
          info = psb_err_internal_error_
          goto 9999
        end select
      end if
    end if

    call psb_realloc(1,itmpv, info)
    if (info /= 0) then 
      write(0,*) 'Error reallocating itmspz'
      goto 9999
    end if
    itmpv(:) = -1
    call psi_bld_tmpovrl(itmpv,desc,info)

  endif

  if (info /= psb_success_) goto 9999
  if (timings) then
    t1 = psb_wtime()
    write(0,*) name,' 1st phase:',t1-t0
    t0 = psb_wtime()
  end if
  ! Finish off 
  lr = desc%indxmap%get_lr()
  call psb_realloc(max(1,lr/2),desc%halo_index, info)
  if (info == psb_success_) call psb_realloc(1,desc%ext_index, info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_realloc')
    Goto 9999
  end if
  desc%halo_index(:)           = -1
  desc%ext_index(:)            = -1
  call psb_cd_set_bld(desc,info)
  if (info /= psb_success_) goto 9999
  if (timings) then
    t1 = psb_wtime()
    write(0,*) name,' 2nd phase:',t1-t0
    t0 = psb_wtime()
  end if

9998 continue
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_cdall
