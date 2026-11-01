#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#include "cudadebug.h"
#include "cudalang.h"
#include <cuda_runtime.h>
#include "core.h"
#include "cublas_v2.h"
#include "mathbase.cuh"
extern "C"
{
	#include "vector.h"
}

#undef GEN_SPGPU_FUNC_NAME
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x), mlt)

#undef GEN_SPGPU_KERN_NAME
#define GEN_SPGPU_KERN_NAME(x, type) CONCAT(CONCAT(CONCAT(CONCAT(spgpu, x), mlt), type), _krn)

#define BLOCK_SIZE 512
#define USE_CUBLAS 1

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, pwise)(int n, VALUE_TYPE* z, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y, VALUE_TYPE beta, int count, int pitchX, int pitchY, int pitchZ, bool conjgX, bool conjgY)
{
	const size_t gridSize = blockDim.x * gridDim.x;
	size_t id = threadIdx.x + BLOCK_SIZE * blockIdx.x;

	//TODO: check beta = 0
	for(; id < n; id += gridSize)
	{
		size_t baseX = 0, baseY = 0, baseZ = 0;
		for(int idx = 0; idx < count; ++idx)
		{
			//z[baseZ + id] = alpha * x[baseX + id] * y[baseY + id] + beta * z[baseZ + id]
			z[baseZ + id] = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(alpha, CONCAT(VALUE_TYPE, _mul)(x[baseX + id], y[baseY + id])), CONCAT(VALUE_TYPE, _mul)(beta, z[baseZ + id]));
			baseX += pitchX;
			baseY += pitchY;
			baseZ += pitchZ;
		}
	}
}

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, dscal)(int n, VALUE_TYPE* z, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y, VALUE_TYPE beta, int count, int pitchY, int pitchZ, bool conjgX, bool conjgY)
{
	const size_t gridSize = blockDim.x * gridDim.x;
	size_t id = threadIdx.x + BLOCK_SIZE * blockIdx.x;

	//TODO: check beta = 0
	for(; id < n; id += gridSize)
	{
		VALUE_TYPE ax = CONCAT(VALUE_TYPE, _mul)(alpha, x[id]);
		size_t baseY = 0, baseZ = 0;
		for(int idx = 0; idx < count; ++idx)
		{
			//z[baseZ + id] = alpha * x[baseX + id] * y[baseY + id] + beta * z[baseZ + id]
			z[baseZ + id] = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(ax, y[baseY + id]), CONCAT(VALUE_TYPE, _mul)(beta, z[baseZ + id]));
			baseY += pitchY;
			baseZ += pitchZ;
		}
	}
}

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, bcast)(int n, VALUE_TYPE* z, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y, VALUE_TYPE beta, int count, int pitchZ, bool conjgX, bool conjgY)
{
	const size_t gridSize = blockDim.x * gridDim.x;
	size_t id = threadIdx.x + BLOCK_SIZE * blockIdx.x;

	//TODO: check beta = 0
	for(; id < n; id += gridSize)
	{
		VALUE_TYPE axy = CONCAT(VALUE_TYPE, _mul)(alpha, CONCAT(VALUE_TYPE, _mul)(x[id], y[id]));
		size_t baseZ = 0;
		for(int idx = 0; idx < count; ++idx)
		{
			//z[baseZ + id] = alpha * x[baseX + id] * y[baseY + id] + beta * z[baseZ + id]
			z[baseZ + id] = CONCAT(VALUE_TYPE, _add)(axy, CONCAT(VALUE_TYPE, _mul)(beta, z[baseZ + id]));
			baseZ += pitchZ;
		}
	}
}

void GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(spgpuHandle_t handle,  __device VALUE_TYPE* z, int countZ, int pitchZ, int n, VALUE_TYPE alpha, __device VALUE_TYPE* x, int countX, int pitchX, 
										__device VALUE_TYPE* y, int countY, int pitchY, VALUE_TYPE beta, bool conjgX, bool conjgY)
{
	dim3 block(BLOCK_SIZE);
	int num_mp         = getGPUMultiProcessors();
	int max_threads_mp = getGPUMaxThreadsPerMP();
	int num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	int num_blocks     = num_blocks_mp * num_mp;
	dim3 grid(num_blocks);

	//TODO: if alpha = 0 call dscal

	if((countZ == countY) && (countY == countX))
		GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, pwise)<<<grid, block, 0, handle->currentStream>>>(n, z, alpha, x, y, beta, countZ, pitchX, pitchY, pitchZ, conjgX, conjgY);
	else if((countZ == countY) && (countY != countX))
		GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, dscal)<<<grid, block, 0, handle->currentStream>>>(n, z, alpha, x, y, beta, countZ, pitchY, pitchZ, conjgX, conjgY);
	else if((countZ != countY) && (countY == countX))
		GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, bcast)<<<grid, block, 0, handle->currentStream>>>(n, z, alpha, x, y, beta, countZ, pitchZ, conjgX, conjgY);
	else
		return;
}