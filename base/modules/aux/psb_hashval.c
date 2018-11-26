#include <stdint.h>
/*
   This is based on the djb2 hashing algorithm
   see e.g. http://www.cse.yorku.ca/~oz/hash.html
*/  

#define IVAL     5381
#define H32MASK  0x7FFFFFFF
#define H64MASK  0x7FFFFFFFFFFFFFFF
#define BMASK    0xFF

int32_t psb_chashval_32(int32_t inkey)
{
  uint32_t key, val, i;
  key = inkey;
  val = IVAL;
  for (i=0; i<4; i++) {
    val = ((val<<5)+val)+(key & BMASK);
    key >>= 8;
  }
  val &= H32MASK;
  return(val);
}

int64_t psb_chashval_64(int64_t inkey)
{
  uint64_t key, val, i;
  key = inkey;
  val = IVAL;
  for (i=0; i<8; i++) {
    val = ((val<<5)+val)+(key & BMASK);
    key >>= 8;
  }
  val &= H64MASK;
  return(val);
}

int32_t psb_chashval_64_32(int64_t inkey)
{
  uint32_t key, val, i;
  key = inkey;
  val = IVAL;
  for (i=0; i<8; i++) {
    val = ((val<<5)+val)+(key & BMASK);
    key >>= 8;
  }
  val &= H32MASK;
  return(val);
}
