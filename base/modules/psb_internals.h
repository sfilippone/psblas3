#ifndef PSB_INTERNALS_H
#define PSB_INTERNALS_H

/* providing a default mangling scheme */
#ifndef LowerUnderscore
#ifndef LowerDoubleUnderscore
#ifndef LowerCase
#ifndef UpperUnderscore
#ifndef UpperDoubleUnderscore 
#ifndef UpperCase
#define LowerUnderscore 1  /* 20110404 the default */
/* #error "should specify a default mangling scheme" */
#endif
#endif
#endif
#endif
#endif
#endif

#endif
