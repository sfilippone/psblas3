#include "elldev.h"

// sparse Ell matrix-vector product
int spmvEllDeviceFloat(void *deviceMat, float* alpha, void* deviceX, float* beta, void* deviceY);
int spmvEllDeviceDouble(void *deviceMat, double* alpha, void* deviceX, double* beta, void* deviceY);


int writeVecDeviceFloat(void* deviceVec, float* hostVec);
int writeVecDeviceDouble(void* deviceVec, double* hostVec);

int readVecDeviceFloat(void* deviceVec, float* hostVec);
int readVecDeviceDouble(void* deviceVec, double* hostVec);

int dotVecDeviceFloat(float* y_res, void* devVecA, void* devVecB);
int dotVecDeviceDouble(double* y_res, void* devVecA, void* devVecB);

int axpbyVecDeviceFloat(float* alpha, void* devVecX, float* beta, void* devVecY);
int axpbyVecDeviceDouble(double* alpha, void* devVecX, double* beta, void* devVecY);




int spmvEllDeviceFloat(void *deviceMat, float* alpha, void* deviceX, float* beta, void* deviceY)
{
#ifdef HAVE_ELL_GPU
  return spmvEllDevice(deviceMat, (void *) alpha, deviceX, (void *) beta, deviceY);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

int spmvEllDeviceDouble(void *deviceMat, double* alpha, void* deviceX, double* beta, void* deviceY)
{
#ifdef HAVE_ELL_GPU
  return spmvEllDevice(deviceMat, (void *) alpha, deviceX, (void *) beta, deviceY);
#else
  return CINTRF_UNSUPPORTED;
#endif
}


int writeVecDeviceFloat(void* deviceVec, float* hostVec)
{
#ifdef HAVE_ELL_GPU
  return writeVecDevice(deviceVec, (void *) hostVec);
#else
    return CINTRF_UNSUPPORTED;
#endif
}

int writeVecDeviceDouble(void* deviceVec, double* hostVec)
{
#ifdef HAVE_ELL_GPU
  return writeVecDevice(deviceVec, (void *) hostVec);
#else
    return CINTRF_UNSUPPORTED;
#endif
}

int readVecDeviceFloat(void* deviceVec, float* hostVec)
{
#ifdef HAVE_ELL_GPU
  return readVecDevice(deviceVec, (void *) hostVec);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

int readVecDeviceDouble(void* deviceVec, double* hostVec)
{
#ifdef HAVE_ELL_GPU
  return readVecDevice(deviceVec, (void *) hostVec);
#else
  return CINTRF_UNSUPPORTED;
#endif

}

int dotVecDeviceFloat(float* y_res, void* devVecA, void* devVecB)
{
#ifdef HAVE_ELL_GPU
  return dotVecDevice((void *) y_res, devVecA, devVecB);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

int dotVecDeviceDouble(double* y_res, void* devVecA, void* devVecB)
{
#ifdef HAVE_ELL_GPU
  return dotVecDevice((void *) y_res, devVecA, devVecB);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

int axpbyVecDeviceFloat(float* alpha, void* devVecX, float* beta, void* devVecY)
{
#ifdef HAVE_ELL_GPU
  return axpbyVecDevice((void *) alpha, devVecX, (void *) beta, devVecY);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

int axpbyVecDeviceDouble(double* alpha, void* devVecX, double* beta, void* devVecY)
{
#ifdef HAVE_ELL_GPU
  return axpbyVecDevice((void *) alpha, devVecX, (void *) beta, devVecY);
#else
  return CINTRF_UNSUPPORTED;
#endif
}

