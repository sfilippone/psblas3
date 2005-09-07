#include <ctype.h>
#define    FTRUE     1
#define    FFALSE    0

#ifdef Add_
#define lsame  lsame_
#endif 

#ifdef UpCase
#define lsame  LSAME
#endif

#ifdef NoChange
#define lsame  lsame
#endif


int lsame(a,b,la,lb)
char *a, *b;
int la,lb;
{
  if ((tolower(*a))==(tolower(*b))) {
    return(FTRUE);
  } else {
    return(FFALSE);
  }
}
