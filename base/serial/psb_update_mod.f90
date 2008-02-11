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
module psb_update_mod

  interface psb_srch_upd
    module procedure psb_d_srch_upd, psb_z_srch_upd
  end interface

  interface coo_srch_upd
    module procedure d_coo_srch_upd, z_coo_srch_upd
  end interface

  interface csr_srch_upd
    module procedure d_csr_srch_upd, z_csr_srch_upd
  end interface

  interface jad_srch_upd
    module procedure d_jad_srch_upd, z_jad_srch_upd
  end interface

contains


  subroutine psb_d_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    info  = 0

    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif
    end if
    select case(tolower(a%fida))
    case ('csr') 
      call  csr_srch_upd(nz,ia,ja,val,nza,a,&
           & imin,imax,jmin,jmax,nzl,info,gtl,ng)
    case ('coo') 
      call  coo_srch_upd(nz,ia,ja,val,nza,a,&
           & imin,imax,jmin,jmax,nzl,info,gtl,ng)
    case ('jad') 
      call  jad_srch_upd(nz,ia,ja,val,nza,a,&
           & imin,imax,jmin,jmax,nzl,info,gtl,ng)
      
    case default
      
      info = -9
      
    end select

  end subroutine psb_d_srch_upd

  subroutine psb_z_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    info  = 0

    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif

      select case(toupper(a%fida))
      case ('CSR') 
!!$      write(0,*) 'Calling csr_srch_upd'
        call  csr_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info,gtl,ng)
!!$      write(0,*) 'From csr_srch_upd:',info
      case ('COO') 
        call  coo_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info,gtl,ng)

      case default

        info = -9

      end select
    else
      select case(toupper(a%fida))
      case ('CSR') 
!!$      write(0,*) 'Calling csr_srch_upd'
        call  csr_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info)
!!$      write(0,*) 'From csr_srch_upd:',info
      case ('COO') 
        call  coo_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info)

      case default

        info = -9

      end select
    end if

  end subroutine psb_z_srch_upd


  subroutine d_csr_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    use psb_error_mod
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    integer  :: debug_level, debug_unit
    character(len=20)    :: name='d_csr_srch_upd'
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,dupl

    info = 0
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif

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
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%ia2(ir)
              i2 = a%ia2(ir+1)
              nc=i2-i1


              if (.true.) then 
                call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
                if (ip>0) then 
                  a%aspk(i1+ip-1) = val(i)
                else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
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
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                  info = i
                  return
                end if

              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
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
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%ia2(ir)
              i2 = a%ia2(ir+1)
              nc = i2-i1
              call issrch(ip,ic,nc,a%ia1(i1:i2-1))
              if (ip>0) then 
                a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
              else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
            end if
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    else

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 

          if ((ir > 0).and.(ir <= a%m)) then 

            i1 = a%ia2(ir)
            i2 = a%ia2(ir+1)
            nc=i2-i1


            if (.true.) then 
              call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
              if (ip>0) then 
                a%aspk(i1+ip-1) = val(i)
              else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if

            else
              ip = -1
              lb = i1
              ub = i2-1
              do
                if (lb > ub) exit
                m = (lb+ub)/2
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
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if
            end if
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if

        end do

      case(psb_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 
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
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select
    end if

  end subroutine d_csr_srch_upd

  subroutine d_coo_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,nnz,dupl
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='d_coo_srch_upd'

    info = 0
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (psb_sp_getifld(psb_srtd_,a,info) /= psb_isrtdcoo_) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = psb_sp_getifld(psb_nnz_,a,info)


    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            if ((ir > 0).and.(ir <= a%m)) then 
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
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              call issrch(ip,ic,nc,a%ia2(i1:i2))
              if (ip>0) then 
                a%aspk(i1+ip-1) = val(i)
              else
                info = i 
                return
              end if
            else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
            endif
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
            if ((ir > 0).and.(ir <= a%m)) then 

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
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              call issrch(ip,ic,nc,a%ia2(i1:i2))
              if (ip>0) then 
                a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
              else
                info = i 
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'              
            end if
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    else

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 

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
            else
              i1 = 1
              i2 = 1
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
          if ((ir > 0).and.(ir <= a%m)) then 

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
            else
              i1 = 1
              i2 = 1
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
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine d_coo_srch_upd

  subroutine d_jad_srch_upd(nz,ia,ja,val,nza,a,imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    
    implicit none

    type(psb_dspmat_type), intent(inout), target :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    integer, pointer     :: ia1(:), ia2(:), ia3(:),&
         & ja_(:), ka_(:)
    integer, allocatable :: indices(:), blks(:), rows(:) 
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ngr,&
         & i,j,nr,dupl, ii, ir, ic

    info = 0
    dupl = psb_sp_getifld(psb_dupl_,a,info)

    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)

    ngr  =  a%ia2(png)              ! the number of blocks
    ja_  => a%ia2(pja:)            ! the array containing the pointers to ka and aspk
    ka_  => a%ia1(:)               ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index 
    ! of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. 
    ! in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. 
    ! in ja to the first csr column


    if (a%pl(1) /= 0) then

      if (present(gtl)) then 
        if (.not.present(ng)) then 
          info = -1
          return
        endif


        allocate(rows(nz),stat=info) 
        if (info /= 0) then 
          info = -4010
          return
        endif
        j=0
        do i=1,nz
          ir = ia(i)
          if ((ir >=1).and.(ir<=ng)) then 
            ir = gtl(ir)
            j = j + 1
            rows(j) = ir
          endif
        enddo
        call psb_msort_unique(rows(1:j),nr)
        allocate(indices(nr),blks(nr),stat=info)
        if (info /= 0) then 
          info = -4010
          return
        endif

        do i=1,nr
          indices(i)=a%pl(rows(i))
          j=0
          blkfnd_gtl: do
            j=j+1
            if(ia1(j) == indices(i)) then
              blks(i)=j
              ipx = ia1(j)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j)+rb
              exit blkfnd_gtl
            else if(ia1(j) > indices(i)) then
              blks(i)=j-1
              ipx = ia1(j-1)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j-1)+rb
              exit blkfnd_gtl
            end if
          end do blkfnd_gtl
        end do


        ! cycle over elements
        update_gtl: do ii=1,nz
          ir = ia(ii)
          ic = ja(ii) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            if ((ir > 0).and.(ir <= a%m)) then 
              ic = gtl(ic) 
              call ibsrch(i,ir,nr,rows)

              ! find which block the row belongs to
              blk = blks(i)

              ! extract first part of the row from the jad block
              ipx = ia1(blk)             ! the first row index of the block
              k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
              rb  = indices(i)-ipx       ! the row offset within the block
              npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

              do  col = ia2(blk), ia3(blk)-1 
                if (ic ==  ka_(ja_(col)+rb)) then 
                  select case(dupl)
                  case(psb_dupl_ovwrt_,psb_dupl_err_)
                    a%aspk(ja_(col)+rb) = val(ii)
                  case(psb_dupl_add_)
                    a%aspk(ja_(col)+rb) = a%aspk(ja_(col)+rb) + val(ii)
                  end select
                  cycle update_gtl
                endif
              end do

              ! extract second part of the row from the csr tail just in case
              row=ia3(blk)+rb
              do j=ja_(row), ja_(row+1)-1
                if (ic ==  ka_(j)) then 
                  select case(dupl)
                  case(psb_dupl_ovwrt_,psb_dupl_err_)
                    a%aspk(j) = val(ii)
                  case(psb_dupl_add_)
                    a%aspk(j) = a%aspk(j) + val(ii)
                  end select
                  cycle update_gtl
                endif
              end do
            end if
          end if
        end do update_gtl

      else

        allocate(rows(nz),stat=info) 
        if (info /= 0) then 
          info = -4010
          return
        endif
        j=0
        do i=1,nz
          ir = ia(i)
          if ((ir >=1).and.(ir<=a%m)) then 
            j = j + 1
            rows(j) = ir
          endif
        enddo
        call psb_msort_unique(rows(1:j),nr)
        allocate(indices(nr),blks(nr),stat=info)
        if (info /= 0) then 
          info = -4010
          return
        endif

        do i=1,nr
          indices(i)=a%pl(rows(i))
          j=0
          blkfnd: do
            j=j+1
            if(ia1(j) == indices(i)) then
              blks(i)=j
              ipx = ia1(j)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j)+rb
              exit blkfnd
            else if(ia1(j) > indices(i)) then
              blks(i)=j-1
              ipx = ia1(j-1)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j-1)+rb
              exit blkfnd
            end if
          end do blkfnd
        end do


        ! cycle over elements
        update: do ii=1,nz
          ir = ia(ii)
          ic = ja(ii) 
          if ((ir >=1).and.(ir<=a%m).and.(ic>=1).and.(ic<=a%k)) then 
            call ibsrch(i,ir,nr,rows)
            ! find which block the row belongs to
            blk = blks(i)

            ! extract first part of the row from the jad block
            ipx = ia1(blk)             ! the first row index of the block
            k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
            rb  = indices(i)-ipx       ! the row offset within the block
            npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

            do  col = ia2(blk), ia3(blk)-1 
              if (ic ==  ka_(ja_(col)+rb)) then 
                select case(dupl)
                case(psb_dupl_ovwrt_,psb_dupl_err_)
                  a%aspk(ja_(col)+rb) = val(ii)
                case(psb_dupl_add_)
                  a%aspk(ja_(col)+rb) = a%aspk(ja_(col)+rb) + val(ii)
                end select
                cycle update
              endif
            end do

            ! extract second part of the row from the csr tail just in case
            row=ia3(blk)+rb
            do j=ja_(row), ja_(row+1)-1
              if (ic ==  ka_(j)) then 
                select case(dupl)
                case(psb_dupl_ovwrt_,psb_dupl_err_)
                  a%aspk(j) = val(ii)
                case(psb_dupl_add_)
                  a%aspk(j) = a%aspk(j) + val(ii)
                end select
                cycle update
              endif
            end do
          end if
        end do update

      end if
    else
      ! There might be some problems
      info=134
    end if

  end subroutine d_jad_srch_upd


  subroutine z_csr_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,dupl

    integer  :: debug_level, debug_unit
    character(len=20)    :: name='z_csr_srch_upd'

    info = 0
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif

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
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%ia2(ir)
              i2 = a%ia2(ir+1)
              nc=i2-i1


              if (.true.) then 
                call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
                if (ip>0) then 
                  a%aspk(i1+ip-1) = val(i)
                else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
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
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                  info = i
                  return
                end if

              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
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
            if ((ir > 0).and.(ir <= a%m)) then 
              i1 = a%ia2(ir)
              i2 = a%ia2(ir+1)
              nc = i2-i1
              call issrch(ip,ic,nc,a%ia1(i1:i2-1))
              if (ip>0) then 
                a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
              else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'
            end if
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    else

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.

        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 

          if ((ir > 0).and.(ir <= a%m)) then 

            i1 = a%ia2(ir)
            i2 = a%ia2(ir+1)
            nc=i2-i1


            if (.true.) then 
              call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
              if (ip>0) then 
                a%aspk(i1+ip-1) = val(i)
              else
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if

            else
              ip = -1
              lb = i1
              ub = i2-1
              do
                if (lb > ub) exit
                m = (lb+ub)/2
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
                  if (debug_level >= psb_debug_serial_) &
                       & write(debug_unit,*) trim(name),&
                       & ': Was searching ',ic,' in: ',i1,i2,&
                       & ' : ',a%ia1(i1:i2-1)
                info = i
                return
              end if
            end if
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if

        end do

      case(psb_dupl_add_)
        ! Add
        ilr = -1 
        ilc = -1 
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 
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
          else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select
    end if

  end subroutine z_csr_srch_upd

  subroutine z_coo_srch_upd(nz,ia,ja,val,nza,a,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    implicit none 

    type(psb_zspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)
    integer  :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nc,nnz,dupl
    integer              :: debug_level, debug_unit
    character(len=20)    :: name='z_coo_srch_upd'


    info = 0
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()
    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (psb_sp_getifld(psb_srtd_,a,info) /= psb_isrtdcoo_) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = psb_sp_getifld(psb_nnz_,a,info)


    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            if ((ir > 0).and.(ir <= a%m)) then 
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
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              call issrch(ip,ic,nc,a%ia2(i1:i2))
              if (ip>0) then 
                a%aspk(i1+ip-1) = val(i)
              else
                info = i 
                return
              end if
            else
            if (debug_level >= psb_debug_serial_) &
                 & write(debug_unit,*) trim(name),&
                 & ': Discarding row that does not belong to us.'
            endif
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
            if ((ir > 0).and.(ir <= a%m)) then 

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
              else
                i1 = 1
                i2 = 1
              end if
              nc = i2-i1+1
              call issrch(ip,ic,nc,a%ia2(i1:i2))
              if (ip>0) then 
                a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
              else
                info = i 
                return
              end if
            else
              if (debug_level >= psb_debug_serial_) &
                   & write(debug_unit,*) trim(name),&
                   & ': Discarding row that does not belong to us.'              
            end if
          end if
        end do

      case default
        info = -3
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    else

      select case(dupl)
      case(psb_dupl_ovwrt_,psb_dupl_err_)
        ! Overwrite.
        ! Cannot test for error, should have been caught earlier.
        do i=1, nz
          ir = ia(i)
          ic = ja(i) 
          if ((ir > 0).and.(ir <= a%m)) then 

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
            else
              i1 = 1
              i2 = 1
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
          if ((ir > 0).and.(ir <= a%m)) then 

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
            else
              i1 = 1
              i2 = 1
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
        if (debug_level >= psb_debug_serial_) &
             & write(debug_unit,*) trim(name),&
             & ': Duplicate handling: ',dupl
      end select

    end if

  end subroutine z_coo_srch_upd



  subroutine z_jad_srch_upd(nz,ia,ja,val,nza,a,imin,imax,jmin,jmax,nzl,info,gtl,ng)

    use psb_spmat_type
    use psb_const_mod
    use psb_realloc_mod
    use psb_string_mod
    use psb_serial_mod
    
    implicit none

    type(psb_zspmat_type), intent(inout), target :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    complex(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)

    integer, pointer     :: ia1(:), ia2(:), ia3(:),&
         & ja_(:), ka_(:)
    integer, allocatable :: indices(:), blks(:), rows(:) 
    integer  :: png, pia, pja, ipx, blk, rb, row, k_pt, npg, col, ngr,&
         & i,j,nr,dupl, ii, ir, ic

    info = 0
    dupl = psb_sp_getifld(psb_dupl_,a,info)

    png = a%ia2(1) ! points to the number of blocks
    pia = a%ia2(2) ! points to the beginning of ia(3,png)
    pja = a%ia2(3) ! points to the beginning of ja(:)

    ngr  =  a%ia2(png)              ! the number of blocks
    ja_  => a%ia2(pja:)            ! the array containing the pointers to ka and aspk
    ka_  => a%ia1(:)               ! the array containing the column indices
    ia1 => a%ia2(pia:pja-1:3)      ! the array containing the first row index 
    ! of each block
    ia2 => a%ia2(pia+1:pja-1:3)    ! the array containing a pointer to the pos. 
    ! in ja to the first jad column
    ia3 => a%ia2(pia+2:pja-1:3)    ! the array containing a pointer to the pos. 
    ! in ja to the first csr column


    if (a%pl(1) /= 0) then

      if (present(gtl)) then 
        if (.not.present(ng)) then 
          info = -1
          return
        endif


        allocate(rows(nz),stat=info) 
        if (info /= 0) then 
          info = -4010
          return
        endif
        j=0
        do i=1,nz
          ir = ia(i)
          if ((ir >=1).and.(ir<=ng)) then 
            ir = gtl(ir)
            j = j + 1
            rows(j) = ir
          endif
        enddo
        call psb_msort_unique(rows(1:j),nr)
        allocate(indices(nr),blks(nr),stat=info)
        if (info /= 0) then 
          info = -4010
          return
        endif

        do i=1,nr
          indices(i)=a%pl(rows(i))
          j=0
          blkfnd_gtl: do
            j=j+1
            if(ia1(j) == indices(i)) then
              blks(i)=j
              ipx = ia1(j)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j)+rb
              exit blkfnd_gtl
            else if(ia1(j) > indices(i)) then
              blks(i)=j-1
              ipx = ia1(j-1)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j-1)+rb
              exit blkfnd_gtl
            end if
          end do blkfnd_gtl
        end do


        ! cycle over elements
        update_gtl: do ii=1,nz
          ir = ia(ii)
          ic = ja(ii) 
          if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
            ir = gtl(ir)
            if ((ir > 0).and.(ir <= a%m)) then 
              ic = gtl(ic) 
              call ibsrch(i,ir,nr,rows)

              ! find which block the row belongs to
              blk = blks(i)

              ! extract first part of the row from the jad block
              ipx = ia1(blk)             ! the first row index of the block
              k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
              rb  = indices(i)-ipx       ! the row offset within the block
              npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

              do  col = ia2(blk), ia3(blk)-1 
                if (ic ==  ka_(ja_(col)+rb)) then 
                  select case(dupl)
                  case(psb_dupl_ovwrt_,psb_dupl_err_)
                    a%aspk(ja_(col)+rb) = val(ii)
                  case(psb_dupl_add_)
                    a%aspk(ja_(col)+rb) = a%aspk(ja_(col)+rb) + val(ii)
                  end select
                  cycle update_gtl
                endif
              end do

              ! extract second part of the row from the csr tail just in case
              row=ia3(blk)+rb
              do j=ja_(row), ja_(row+1)-1
                if (ic ==  ka_(j)) then 
                  select case(dupl)
                  case(psb_dupl_ovwrt_,psb_dupl_err_)
                    a%aspk(j) = val(ii)
                  case(psb_dupl_add_)
                    a%aspk(j) = a%aspk(j) + val(ii)
                  end select
                  cycle update_gtl
                endif
              end do
            end if
          end if
        end do update_gtl

      else

        allocate(rows(nz),stat=info) 
        if (info /= 0) then 
          info = -4010
          return
        endif
        j=0
        do i=1,nz
          ir = ia(i)
          if ((ir >=1).and.(ir<=a%m)) then 
            j = j + 1
            rows(j) = ir
          endif
        enddo
        call psb_msort_unique(rows(1:j),nr)
        allocate(indices(nr),blks(nr),stat=info)
        if (info /= 0) then 
          info = -4010
          return
        endif

        do i=1,nr
          indices(i)=a%pl(rows(i))
          j=0
          blkfnd: do
            j=j+1
            if(ia1(j) == indices(i)) then
              blks(i)=j
              ipx = ia1(j)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j)+rb
              exit blkfnd
            else if(ia1(j) > indices(i)) then
              blks(i)=j-1
              ipx = ia1(j-1)         ! the first row index of the block
              rb  = indices(i)-ipx   ! the row offset within the block
              row = ia3(j-1)+rb
              exit blkfnd
            end if
          end do blkfnd
        end do


        ! cycle over elements
        update: do ii=1,nz
          ir = ia(ii)
          ic = ja(ii) 
          if ((ir >=1).and.(ir<=a%m).and.(ic>=1).and.(ic<=a%k)) then 
            call ibsrch(i,ir,nr,rows)
            ! find which block the row belongs to
            blk = blks(i)

            ! extract first part of the row from the jad block
            ipx = ia1(blk)             ! the first row index of the block
            k_pt= ia2(blk)             ! the pointer to the beginning of a column in ja
            rb  = indices(i)-ipx       ! the row offset within the block
            npg = ja_(k_pt+1)-ja_(k_pt)  ! the number of rows in the block

            do  col = ia2(blk), ia3(blk)-1 
              if (ic ==  ka_(ja_(col)+rb)) then 
                select case(dupl)
                case(psb_dupl_ovwrt_,psb_dupl_err_)
                  a%aspk(ja_(col)+rb) = val(ii)
                case(psb_dupl_add_)
                  a%aspk(ja_(col)+rb) = a%aspk(ja_(col)+rb) + val(ii)
                end select
                cycle update
              endif
            end do

            ! extract second part of the row from the csr tail just in case
            row=ia3(blk)+rb
            do j=ja_(row), ja_(row+1)-1
              if (ic ==  ka_(j)) then 
                select case(dupl)
                case(psb_dupl_ovwrt_,psb_dupl_err_)
                  a%aspk(j) = val(ii)
                case(psb_dupl_add_)
                  a%aspk(j) = a%aspk(j) + val(ii)
                end select
                cycle update
              endif
            end do
          end if
        end do update

      end if
    else
      ! There might be some problems
      info=134
    end if

  end subroutine z_jad_srch_upd


end module psb_update_mod

