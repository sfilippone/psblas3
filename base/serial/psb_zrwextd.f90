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
! File:  psb_zrwextd.f90 
! Subroutine: 
! Parameters:
!
! We have a problem here: 1. How to handle well all the formats? 
!                         2. What should we do with rowscale? Does it only 
!                            apply when a%fida='COO' ?????? 
!
!
subroutine psb_zrwextd(nr,a,info,b,rowscale)
  use psb_spmat_type
  use psb_error_mod
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer, intent(in)                            :: nr
  type(psb_zspmat_type), intent(inout)           :: a
  integer,intent(out)                            :: info
  type(psb_zspmat_type), intent(in), optional    :: b
  logical,intent(in), optional                   :: rowscale

  integer :: i,j,ja,jb,err_act,nza,nzb
  character(len=20)                 :: name, ch_err
  logical  rowscale_ 

  name='psb_zrwextd'
  info  = 0
  call psb_erractionsave(err_act)

  if (present(rowscale)) then 
    rowscale_ = rowscale
  else
    rowscale_ = .true.
  end if

  if (nr > a%m) then 
    if (a%fida == 'CSR') then 
      call psb_realloc(nr+1,a%ia2,info)
      if (present(b)) then 
        nzb = psb_sp_get_nnzeros(b)
        call psb_realloc(size(a%ia1)+nzb,a%ia1,info)
        call psb_realloc(size(a%aspk)+nzb,a%aspk,info)
        if (b%fida=='CSR') then 

          do i=1, min(nr-a%m,b%m)
            a%ia2(a%m+i+1) =  a%ia2(a%m+i) + b%ia2(i+1) - b%ia2(i)
            ja = a%ia2(a%m+i)
            jb = b%ia2(i)
            do 
              if (jb >=  b%ia2(i+1)) exit
              a%aspk(ja) = b%aspk(jb)
              a%ia1(ja) = b%ia1(jb)
              ja = ja + 1
              jb = jb + 1
            end do
          end do
          do j=i,nr-a%m
            a%ia2(a%m+i+1) = a%ia2(a%m+i)
          end do
        else 
          write(0,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'
        endif
      else
        do i=a%m+2,nr+1
          a%ia2(i) = a%ia2(i-1)
        end do
      end if
      a%m = nr
    else if (a%fida == 'COO') then 
      if (present(b)) then 
        nza = psb_sp_get_nnzeros(a)
        nzb = psb_sp_get_nnzeros(b)
        call psb_sp_reall(a,nza+nzb,info)
        if (b%fida=='COO') then 
          if (rowscale_) then 
            do j=1,nzb
              if ((a%m + b%ia1(j)) <= nr) then 
                a%ia1(nza+j)  = a%m + b%ia1(j)
                a%ia2(nza+j)  = b%ia2(j)
                a%aspk(nza+j) = b%aspk(j)
              end if
            enddo
          else
            do j=1,nzb
              if ((b%ia1(j)) <= nr) then 
                a%ia1(nza+j)  = b%ia1(j)
                a%ia2(nza+j)  = b%ia2(j)
                a%aspk(nza+j) = b%aspk(j)
              endif
            enddo
          endif
          a%infoa(psb_nnz_) = nza+nzb
        else if(b%fida=='CSR') then 
          do i=1, min(nr-a%m,b%m)
            do 
              jb = b%ia2(i)
              if (jb >=  b%ia2(i+1)) exit
              nza = nza + 1 
              a%aspk(nza) = b%aspk(jb)
              a%ia1(nza)  = a%m + i
              a%ia2(nza)  = b%ia1(jb)
              jb = jb + 1
            end do
          end do
          a%infoa(psb_nnz_) = nza
        else
          write(0,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'
        endif
      endif
      a%m = nr
    else if (a%fida == 'JAD') then 
       info=135
       ch_err=a%fida(1:3)
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    else
       info=136
       ch_err=a%fida(1:3)
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end subroutine psb_zrwextd
