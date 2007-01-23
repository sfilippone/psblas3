!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
! File:  psb_dcssm.f90 
! Subroutine: 
! Parameters:

subroutine psb_dcssm(alpha,t,b,beta,c,info,trans,unitd,d)
  use psb_spmat_type
  use psb_error_mod
  implicit none

  type(psb_dspmat_type) :: t
  real(kind(1.d0))      :: alpha, beta, b(:,:), c(:,:)
  integer               :: info
  character, optional   :: trans, unitd
  real(kind(1.d0)), optional, target :: d(:)
  
  real(kind(1.d0)), allocatable :: work(:)
  real(kind(1.d0)), pointer :: ddl(:)
  character :: lt, lu
  integer   :: iwsz,m,n,lb,lc,err_act
  character(len=20) :: name

  name='psb_dcssm'
  info  = 0
  call psb_erractionsave(err_act)

  
  if (present(trans)) then 
    lt = trans
  else
    lt = 'N'
  endif
  if (present(unitd)) then 
    lu = unitd
  else
    lu = 'U'
  endif
  if (present(d)) then 
    ddl => d
  else
    allocate(ddl(1))
  endif

  m = t%m
  n = min(size(b,2),size(c,2))
  lb = size(b,1)
  lc = size(c,1)
  iwsz = 2*m*n
  allocate(work(iwsz))
  
  call dcssm(lt,m,n,alpha,lu,ddl,&
       & t%pl,t%fida,t%descra,t%aspk,t%ia1,t%ia2,t%infoa,t%pr,&
       & b,lb,beta,c,lc,work,iwsz,info)
  
  if (.not.present(d)) then 
    deallocate(ddl)
  endif
  deallocate(work)
  call psb_erractionrestore(err_act)

  if(info.ne.0) then
     if (err_act.eq.psb_act_abort_) then
        call psb_error()
        return
     end if
  end if
     
  return

end subroutine psb_dcssm
