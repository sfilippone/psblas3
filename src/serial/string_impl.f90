function  tolowerc(string)
  character(len=*), intent(in)  :: string
  character(len=len(string))    :: tolowerc
  character(len=*), parameter   :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter   :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer  :: i,k

  do i=1,len(string)
    k = index(ucase,string(i:i))
    if (k /=0 ) then 
      tolowerc(i:i) = lcase(k:k)
    else          
      tolowerc(i:i) = string(i:i)
    end if
  enddo
end function tolowerc
function  toupperc(string)
  character(len=*), intent(in)  :: string
  character(len=len(string))    :: toupperc
  character(len=*), parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer  :: i,k

  do i=1,len(string)
    k = index(lcase,string(i:i))
    if (k /=0 ) then 
      toupperc(i:i) = ucase(k:k)
    else          
      toupperc(i:i) = string(i:i)
    end if
  enddo
end function toupperc

subroutine   sub_toupperc(string,strout)
  character(len=*), intent(in)  :: string
  character(len=*), intent(out)  :: strout
  character(len=*), parameter :: lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer  :: i,k

  do i=1,len(string)
    k = index(lcase,string(i:i))
    if (k /=0 ) then 
      strout(i:i) = ucase(k:k)
    else          
      strout(i:i) = string(i:i)
    end if
  enddo
end subroutine sub_toupperc
