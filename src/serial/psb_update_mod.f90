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
module psb_update_mod


  use psb_spmat_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_string_mod

  interface psb_srch_upd
    module procedure psb_d_srch_upd, psb_z_srch_upd
  end interface

  interface coo_srch_upd
    module procedure d_coo_srch_upd, z_coo_srch_upd
  end interface

  interface csr_srch_upd
    module procedure d_csr_srch_upd, z_csr_srch_upd
  end interface

contains


  subroutine psb_d_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    
    info  = 0

    select case(toupper(a%fida))
    case ('CSR') 
!!$      write(0,*) 'Calling csr_srch_upd'
      call  csr_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
           & imin,imax,jmin,jmax,nzl,info)
!!$      write(0,*) 'From csr_srch_upd:',info
    case ('COO') 
      call  coo_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
           & imin,imax,jmin,jmax,nzl,info)

    case default

      info = -9

    end select

  end subroutine psb_d_srch_upd

  subroutine psb_z_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info

    info  = 0

    select case(toupper(a%fida))
    case ('CSR') 
!!$      write(0,*) 'Calling csr_srch_upd'
      call  csr_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
           & imin,imax,jmin,jmax,nzl,info)
!!$      write(0,*) 'From csr_srch_upd:',info
    case ('COO') 
      call  coo_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
           & imin,imax,jmin,jmax,nzl,info)

    case default

      info = -9

    end select

  end subroutine psb_z_srch_upd


  subroutine d_csr_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.

      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc=i2-i1


          if (.true.) then 
            call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
            if (ip>0) then 
              a%aspk(i1+ip-1) = val(i)
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = i
              return
            end if

          else
!!$ 
            ip = -1
            lb = i1
            ub = i2-1
            do
              if (lb > ub) exit
              m = (lb+ub)/2
!!$              write(0,*) 'Debug: ',lb,m,ub
              if (ic == a%ia1(m))  then
                ip   = m 
                lb   = ub + 1
              else if (ic < a%ia1(m))  then
                ub = m-1
              else 
                lb = m + 1
              end if
            enddo

            if (ip>0) then 
              a%aspk(ip) = val(i)                
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = i
              return
            end if

          end if
        end if
      end do

    case(psb_dupl_add_)
      ! Add
      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc = i2-i1
          call issrch(ip,ic,nc,a%ia1(i1:i2-1))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = i
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select
  end subroutine d_csr_srch_upd

  subroutine d_coo_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl,isrt

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (psb_sp_getifld(psb_srtd_,a,info) /= psb_isrtdcoo_) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = psb_sp_getifld(psb_nnz_,a,info)


    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit
              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = val(i)
          else
            info = i 
            return
          end if
        end if
      end do
    case(psb_dupl_add_)
      ! Add
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit
              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = i 
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select

  end subroutine d_coo_srch_upd



  subroutine z_csr_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.

      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc=i2-i1


          if (.true.) then 
            call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
            if (ip>0) then 
              a%aspk(i1+ip-1) = val(i)
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = i
              return
            end if

          else
!!$ 
            ip = -1
            lb = i1
            ub = i2-1
            do
              if (lb > ub) exit
              m = (lb+ub)/2
!!$              write(0,*) 'Debug: ',lb,m,ub
              if (ic == a%ia1(m))  then
                ip   = m 
                lb   = ub + 1
              else if (ic < a%ia1(m))  then
                ub = m-1
              else 
                lb = m + 1
              end if
            enddo

            if (ip>0) then 
              a%aspk(ip) = val(i)                
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = i
              return
            end if

          end if
        end if
      end do

    case(psb_dupl_add_)
      ! Add
      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc = i2-i1
          call issrch(ip,ic,nc,a%ia1(i1:i2-1))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = i
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select
  end subroutine z_csr_srch_upd

  subroutine z_coo_srch_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl,isrt

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (psb_sp_getifld(psb_srtd_,a,info) /= psb_isrtdcoo_) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = psb_sp_getifld(psb_nnz_,a,info)


    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit

              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = val(i)
          else
            info = i
            return
          end if
        end if
      end do
    case(psb_dupl_add_)
      ! Add
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit
              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = i
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select

  end subroutine z_coo_srch_upd


end module psb_update_mod

