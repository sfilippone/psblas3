#include <stdlib.h>
#include <stdio.h>

#include "cintrf.h"
#include "vectordev.h"


#define VALUE_TYPE cuDoubleComplex
#define TYPE_SYMBOL z
#include "psi_cuda_CopyCooToElg.cuh"
