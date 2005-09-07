      subroutine lowerc(string,pos,len)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c Convert uppercase letters to lowercase letters in string with
c starting postion pos and length len.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      integer pos, len
      character*(*) string

      character*26 lcase, ucase
      save lcase,ucase
      data lcase/'abcdefghijklmnopqrstuvwxyz'/
      data ucase/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      do i=pos,len
        k = index(ucase,string(i:i))
        if (k.ne.0) string(i:i) = lcase(k:k)
      enddo
      return
      end
