#ifndef PSB_UTIL_CBIND_
#define PSB_UTIL_CBIND_

#include "psb_c_sutil.h"
#include "psb_c_dutil.h"
#include "psb_c_cutil.h"
#include "psb_c_zutil.h"

psb_i_t psb_c_i_idx2ijk(psb_i_t ijk[],psb_i_t idx,psb_i_t sizes[],psb_i_t modes,psb_i_t base);
psb_i_t psb_c_l_idx2ijk(psb_i_t ijk[],psb_l_t idx,psb_i_t sizes[],psb_i_t modes,psb_i_t base);
psb_i_t psb_c_i_ijk2idx(psb_i_t ijk[],psb_i_t sizes[],psb_i_t modes,psb_i_t base);
psb_l_t psb_c_l_ijk2idx(psb_i_t ijk[],psb_i_t sizes[],psb_i_t modes,psb_i_t base);
psb_i_t psb_c_dist1didx(psb_i_t n, psb_i_t np, psb_i_t base, psb_i_t indexsize, psb_i_t v[]);

#endif
