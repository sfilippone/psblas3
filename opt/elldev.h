#ifndef SPGPU_INTERFACE_H
#define SPGPU_INTERFACE_H

////////////// 
// legenda:
// cM : compressed Matrix
// rP : row pointers
// rS : row size
/////////////

// element types
#define TYPE_FLOAT 0
#define TYPE_DOUBLE 1
// TYPE_COMPLEX?
// TYPE_INT?
// TYPE_BOOLEAN?

// return codes
#define CINTRF_SUCCESS 0
#define CINTRF_NOMEMORY -1
#define CINTRF_UNSUPPORTED -2

typedef struct EllDeviceParams
{			
	// The resulting allocation for cM and rP will be pitch*maxRowSize*elementSize
	unsigned int elementType;
	
	// Pitch (in number of elements)
	unsigned int pitch;

	// Number of rows.
	// Used to allocate rS array
	unsigned int rows; 
		
	// Largest row size
	unsigned int maxRowSize;
	
	// First index (e.g 0 or 1)
	unsigned int firstIndex;
} EllDeviceParams;

#ifdef HAVE_ELL_GPU
// Generate ELLPACK format matrix parameters
EllDeviceParams getEllDeviceParams(unsigned int rows, unsigned int maxRowSize, unsigned int elementType, unsigned int firstIndex);

// Allocate/Free matrix on device
// remote pointer returned in *remoteMatrix
// (device struct type is hidden, use device pointer as id)

// can return CINTRF_SUCCESS or CINTRF_NOMEMORY or CINTRF_UNSUPPORTED
// return the matrix pointer (use the pointer just as an id) in *deviceMat
int allocEllDevice(void** deviceMat, EllDeviceParams* params);
void freeEllDevice(void* deviceMat);

// Update device copy with host copy
int writeEllDevice(void *deviceMat, void* cM, int* rP, int* rS);

// Update host copy with device copy
int readEllDevice(void *deviceMat, void* cM, int* rP, int* rS);

// sparse Ell matrix-vector product
int spmvEllDevice(void *deviceMat, void* alpha, void* deviceX, void* beta, void* deviceY);

//// Vector management
// can return CINTRF_SUCCESS or CINTRF_NOMEMORY or CINTRF_UNSUPPORTED
// return the vector pointer (use the pointer just as an id) in *deviceVec
int allocVecDevice(void** deviceVec, int size, unsigned int elementType);
void freeVecDevice(void* deviceVec);

// Host Vector -> Device Vector
// if deviceVec is a float device vector, hostVec will be a float array
// if deviceVec is a double device vector, hostVec will be a double array
int writeVecDevice(void* deviceVec, void* hostVec);
// Device Vector -> Host Vector
int readVecDevice(void* deviceVec, void* hostVec);

// dot product (y_res = a * b)
// y_res: pointer to result (e.g. float*/double*)
// devVecA, devVecB: device vectors
// if devVecA and devVecB are float prec. device vectors, y_res should be a pointer to a float value
// if devVecA and devVecB are double prec. device vectors, y_res should be a pointer to a double value
int dotVecDevice(void* y_res, void* devVecA, void* devVecB);

// if devVecX and devVecY are float prec. device vectors, alpha and beta should be pointers to a float value 
// if devVecX and devVecY are double prec. device vectors, alpha and beta should be pointers to a double value 
int axpbyVecDevice(void* alpha, void* devVecX, void* beta, void* devVecY);
#endif
#endif
