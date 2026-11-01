#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#include "cudadebug.h"
#include "cudalang.h"
#include <cuda_runtime.h>
#include "core.h"
#include "cublas_v2.h"
#include "mathbase.cuh"

#undef GEN_SPGPU_FUNC_NAME1
#define GEN_SPGPU_FUNC_NAME1(x) CONCAT(CONCAT(spgpu,x), axpby)

#undef GEN_SPGPU_FUNC_NAME2
#define GEN_SPGPU_FUNC_NAME2(x) CONCAT(CONCAT(spgpu,x), axpbycz)

#undef GEN_SPGPU_MFUNC_NAME
#define GEN_SPGPU_MFUNC_NAME(x) CONCAT(CONCAT(spgpu,x), maxpby)

#undef GEN_SPGPU_KERN_NAME1
#define GEN_SPGPU_KERN_NAME1(x) CONCAT(CONCAT(CONCAT(spgpu, x), axpby), _krn)

#undef GEN_SPGPU_KERN_NAME2
#define GEN_SPGPU_KERN_NAME2(x) CONCAT(CONCAT(CONCAT(spgpu, x), axpbycz), _krn)

#define BLOCK_SIZE 512

#if 1
__global__ void GEN_SPGPU_KERN_NAME1(TYPE_SYMBOL)(VALUE_TYPE* z, int n, VALUE_TYPE beta, VALUE_TYPE* y, VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE * blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	if(CONCAT(VALUE_TYPE, _isZero)(beta))
		for(; id < n; id += gridSize)
	    	z[id] = CONCAT(VALUE_TYPE, _mul)(alpha, x[id]);
	else
	  	for(; id < n; id += gridSize)
			z[id] = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(alpha, x[id]), CONCAT(VALUE_TYPE, _mul)(beta, y[id]));
}

void GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL)(spgpuHandle_t handle,
	__device VALUE_TYPE* z,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE* y,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x)
{
	int num_mp         = getGPUMultiProcessors();
	int max_threads_mp = getGPUMaxThreadsPerMP();
	int num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	int num_blocks     = num_blocks_mp * num_mp;
	dim3 block(BLOCK_SIZE);
	dim3 grid(num_blocks);

	GEN_SPGPU_KERN_NAME1(TYPE_SYMBOL)<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

__global__ void GEN_SPGPU_KERN_NAME2(TYPE_SYMBOL)(VALUE_TYPE* w, int n, VALUE_TYPE gamma, VALUE_TYPE* z, VALUE_TYPE beta, VALUE_TYPE* y, VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE * blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	//TODO: alpha, beta, gamma == 0/1
	for(; id < n; id += gridSize)
		w[id] = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(alpha, x[id]), CONCAT(VALUE_TYPE, _mul)(beta, y[id])), CONCAT(VALUE_TYPE, _mul)(gamma, z[id]));
}

void GEN_SPGPU_FUNC_NAME2(TYPE_SYMBOL)(spgpuHandle_t handle,
	__device VALUE_TYPE* w,
	int n,
	VALUE_TYPE gamma,
	__device VALUE_TYPE* z,
	VALUE_TYPE beta,
	__device VALUE_TYPE* y,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x)
{
	int num_mp         = getGPUMultiProcessors();
	int max_threads_mp = getGPUMaxThreadsPerMP();
	int num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	int num_blocks     = num_blocks_mp * num_mp;
	dim3 block(BLOCK_SIZE);
	dim3 grid(num_blocks);

	GEN_SPGPU_KERN_NAME2(TYPE_SYMBOL)<<<grid, block, 0, handle->currentStream>>>(w, n, gamma, z, beta, y, alpha, x);
}


#else

__global__ void GEN_SPGPU_KERN_NAME1(TYPE_SYMBOL)(VALUE_TYPE *z, int n, VALUE_TYPE beta, VALUE_TYPE *y, VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if(id < n)
	{
		// Since z, x and y are accessed with the same offset by the same thread,
		// and the write to z follows the x and y read, x, y and z can share the same base address (in-place computing).

		if(CONCAT(VALUE_TYPE, _isZero)(beta))
			z[id] = CONCAT(VALUE_TYPE, _mul)(alpha,x[id]);
		else
			z[id] = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(alpha, x[id]), CONCAT(VALUE_TYPE, _mul)(beta,y[id]));
	}
}


void CONCAT(GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL), _)(spgpuHandle_t handle,
	__device VALUE_TYPE *z,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE *y,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);

	GEN_SPGPU_KERN_NAME1(TYPE_SYMBOL)<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

void GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL)(spgpuHandle_t handle,
	__device VALUE_TYPE *z,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE *y,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	while (n > maxNForACall) //managing large vectors
	{
		CONCAT(GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL), _)(handle, z, maxNForACall, beta, y, alpha, x);
		
		x = x + maxNForACall;
		y = y + maxNForACall;
		z = z + maxNForACall;
		n -= maxNForACall;
	}
	
	CONCAT(GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL), _)(handle, z, n, beta, y, alpha, x);

	cudaCheckError("CUDA error on daxpby");
}

#endif

void GEN_SPGPU_MFUNC_NAME(TYPE_SYMBOL)(spgpuHandle_t handle,
	__device VALUE_TYPE* z,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE* y,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x, 
	int count, int pitch)
{
  for(int i = 0; i < count; ++i)
    GEN_SPGPU_FUNC_NAME1(TYPE_SYMBOL)(handle, z + (pitch * i), n, beta, y + (pitch * i), alpha, x + (pitch * i));  
}
