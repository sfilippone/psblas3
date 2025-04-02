#ifndef PSB_TYPES_H
#define PSB_TYPES_H


#include <stdint.h>

#ifdef  __cplusplus
#include <complex>
#else
#include <math.h>
#include <complex.h>
#endif
  typedef int32_t psb_m_t;

#if defined(PSB_IPK4) &&  defined(PSB_LPK4)
  typedef int32_t psb_i_t;
  typedef int32_t psb_l_t;
#elif defined(PSB_IPK4) &&  defined(PSB_LPK8)
  typedef int32_t psb_i_t;
  typedef int64_t psb_l_t;
#elif defined(PSB_IPK8) &&  defined(PSB_LPK8)
  typedef int64_t psb_i_t;
  typedef int64_t psb_l_t;
#else
#endif
  typedef int64_t psb_e_t;

  typedef float  psb_s_t;
  typedef double psb_d_t;

#ifdef  __cplusplus
   typedef std::complex<float> psb_c_t;
   typedef std::complex<double> psb_z_t;
#else
   typedef float complex psb_c_t;
   typedef double complex psb_z_t;
#endif
#endif
