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

subroutine psi_i_inner_cnvs(x,hashmask,hashv,glb_lc)
  use psi_mod, psi_protect_name => psi_i_inner_cnvs

  integer(psb_ipk_), intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
  integer(psb_ipk_), intent(inout) :: x

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices 
  ! are kept in a (hashed) list of ordered lists. 
  ! Thus we first hash the index, then we do a binary search on the 
  ! ordered sublist. The hashing is based on the low-order bits 
  ! for a width of psb_hash_bits 
  !

  key = x
  ih  = iand(key,hashmask)
  idx = hashv(ih)
  nh  = hashv(ih+1) - hashv(ih) 
  if (nh > 0) then 
    tmp = -1 
    lb = idx
    ub = idx+nh-1
    do 
      if (lb>ub) exit
      lm = (lb+ub)/2
      if (key == glb_lc(lm,1)) then 
        tmp = lm
        exit
      else if (key<glb_lc(lm,1)) then 
        ub = lm - 1
      else
        lb = lm + 1
      end if
    end do
  else 
    tmp = -1
  end if
  if (tmp > 0) then 
    x = glb_lc(tmp,2)
  else         
    x = tmp 
  end if
end subroutine psi_i_inner_cnvs

subroutine psi_i_inner_cnvs2(x,y,hashmask,hashv,glb_lc)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnvs2
  integer(psb_ipk_), intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
  integer(psb_ipk_), intent(in)  :: x
  integer(psb_ipk_), intent(out) :: y

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices 
  ! are kept in a (hashed) list of ordered lists. 
  ! Thus we first hash the index, then we do a binary search on the 
  ! ordered sublist. The hashing is based on the low-order bits 
  ! for a width of psb_hash_bits 
  !

  key = x
  ih  = iand(key,hashmask)
  idx = hashv(ih)
  nh  = hashv(ih+1) - hashv(ih) 
  if (nh > 0) then 
    tmp = -1 
    lb = idx
    ub = idx+nh-1
    do 
      if (lb>ub) exit
      lm = (lb+ub)/2
      if (key == glb_lc(lm,1)) then 
        tmp = lm
        exit
      else if (key<glb_lc(lm,1)) then 
        ub = lm - 1
      else
        lb = lm + 1
      end if
    end do
  else 
    tmp = -1
  end if
  if (tmp > 0) then 
    y = glb_lc(tmp,2)
  else         
    y = tmp 
  end if
end subroutine psi_i_inner_cnvs2


subroutine psi_i_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnv1
  integer(psb_ipk_), intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
  logical, intent(in), optional    :: mask(:)
  integer(psb_ipk_), intent(inout) :: x(:)

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices 
  ! are kept in a (hashed) list of ordered lists. 
  ! Thus we first hash the index, then we do a binary search on the 
  ! ordered sublist. The hashing is based on the low-order bits 
  ! for a width of psb_hash_bits 
  !
  if (present(mask)) then 
    do i=1, n
      if (mask(i)) then 
        key = x(i) 
        ih  = iand(key,hashmask)
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih) 
        if (nh > 0) then 
          tmp = -1 
          lb = idx
          ub = idx+nh-1
          do 
            if (lb>ub) exit
            lm = (lb+ub)/2
            if (key == glb_lc(lm,1)) then 
              tmp = lm
              exit
            else if (key<glb_lc(lm,1)) then 
              ub = lm - 1
            else
              lb = lm + 1
            end if
          end do
        else 
          tmp = -1
        end if
        if (tmp > 0) then 
          x(i) = glb_lc(tmp,2)
        else         
          x(i) = tmp 
        end if
      end if
    end do
  else
    do i=1, n
      key = x(i) 
      ih  = iand(key,hashmask)
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih) 
      if (nh > 0) then 
        tmp = -1 
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm,1)) then 
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else 
        tmp = -1
      end if
      if (tmp > 0) then 
        x(i) = glb_lc(tmp,2)
      else         
        x(i) = tmp 
      end if
    end do
  end if
end subroutine psi_i_inner_cnv1

subroutine psi_i_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask)
  use psi_mod, psi_protect_name =>  psi_i_inner_cnv2
  integer(psb_ipk_), intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
  logical, intent(in),optional  :: mask(:)
  integer(psb_ipk_), intent(in)  :: x(:)
  integer(psb_ipk_), intent(out) :: y(:)

  integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
  !
  ! When a large descriptor is assembled the indices 
  ! are kept in a (hashed) list of ordered lists. 
  ! Thus we first hash the index, then we do a binary search on the 
  ! ordered sublist. The hashing is based on the low-order bits 
  ! for a width of psb_hash_bits 
  !
  if (present(mask)) then 
    do i=1, n
      if (mask(i)) then 
        key = x(i) 
        ih  = iand(key,hashmask)
        if (ih > ubound(hashv,1) ) then 
          write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
        end if
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih) 
        if (nh > 0) then 
          tmp = -1 
          lb = idx
          ub = idx+nh-1
          do 
            if (lb>ub) exit
            lm = (lb+ub)/2
            if (key == glb_lc(lm,1)) then 
              tmp = lm
              exit
            else if (key<glb_lc(lm,1)) then 
              ub = lm - 1
            else
              lb = lm + 1
            end if
          end do
        else 
          tmp = -1
        end if
        if (tmp > 0) then 
          y(i) = glb_lc(tmp,2)
        else         
          y(i) = tmp 
        end if
      end if
    end do
  else
    do i=1, n
      key = x(i) 
      ih  = iand(key,hashmask)
      if (ih > ubound(hashv,1) ) then 
        write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
      end if
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih) 
      if (nh > 0) then 
        tmp = -1 
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm,1)) then 
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else 
        tmp = -1
      end if
      if (tmp > 0) then 
        y(i) = glb_lc(tmp,2)
      else         
        y(i) = tmp 
      end if
    end do
  end if
end subroutine psi_i_inner_cnv2
