#ifndef PSB_UTIL_CBIND_
#define PSB_UTIL_CBIND_

#include "psb_c_sutil.h"
#include "psb_c_dutil.h"
#include "psb_c_cutil.h"
#include "psb_c_zutil.h"

psb_i_t psb_c_idx2ijk(psb_i_t *i, psb_i_t *j, psb_i_t idx, psb_i_t nx, psb_i_t ny, psb_i_t base );
psb_i_t psb_c_lidx2ijk(psb_i_t *i, psb_i_t *j, psb_l_t idx, psb_i_t nx, psb_i_t ny, psb_i_t base );


#endif
