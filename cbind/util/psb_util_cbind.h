#ifndef PSB_UTIL_CBIND_
#define PSB_UTIL_CBIND_

#include "psb_c_sutil.h"
#include "psb_c_dutil.h"
#include "psb_c_cutil.h"
#include "psb_c_zutil.h"

psb_c_i_t psb_c_i_idx2ijk(psb_c_i_t ijk[],psb_c_i_t idx,psb_c_i_t sizes[],psb_c_i_t modes,psb_c_i_t base);
psb_c_i_t psb_c_l_idx2ijk(psb_c_i_t ijk[],psb_c_l_t idx,psb_c_i_t sizes[],psb_c_i_t modes,psb_c_i_t base);
psb_c_i_t psb_c_i_ijk2idx(psb_c_i_t ijk[],psb_c_i_t sizes[],psb_c_i_t modes,psb_c_i_t base);
psb_c_l_t psb_c_l_ijk2idx(psb_c_i_t ijk[],psb_c_i_t sizes[],psb_c_i_t modes,psb_c_i_t base);
psb_c_i_t psb_c_dist1didx(psb_c_i_t n, psb_c_i_t np, psb_c_i_t base, psb_c_i_t indexsize, psb_c_i_t v[]);

#endif
