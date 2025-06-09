# Computational Routines Test
This is a directory containing all the tests done in order to analyze the correctness of the computational routines present in PSBLAS. 

## Test Environment
These tests are developed using a linux environment, in particular Rocky Linux 9.5 (Blue Onyx). 

The compiler used is:
- gnu 12.2.1

The necessary dependnces are:
- mpich 4.2.2
- PSBLAS 3.9
- CUDA 12.5

## Routines
In this test suite were considered only computational routines implemented by PSBLAS, according to the version 3.9 of the documentation. In the following table are reported all the kernels, their implementation and wheter or not they were tested yet.
|**Kernel**| **PSBLAS Subroutine**|**Description**|**Test**|
| ------------------------------- | :--------------------------: | ---------------------------------------------------------------------- | :---------------: |
|**General Dense Matrix Sum**| `psb_geaxpby`| This subroutine is an interface to the computational kernel for dense matrix sum: $$Y \leftarrow \alpha  X + \beta  Y $$|Work in progress  :hammer_and_wrench:|
| **Dot product**|`psb_gedot`|This function computes dot product between two vectors x and y. $$ dot \leftarrow x^T y $$ If x and y are real vectors it computes dot-product as: $$ dot \leftarrow x^H y $$|No ❌|
| **Generalized Dot Product** |`psb_gedots`|This subroutine computes a series of dot products among the columns of two dense matrices x and y: $$ res(i) \leftarrow x(:,i)^T  y(:,i) $$ If the matrices are complex, then the usual convention applies, i.e. the conjugate transpose of x is used. If x and y are of rank one, then res is a scalar, else it is a rank one array.|No ❌|
|**Infinity-Norm of Vector**|`psb_normi`/`psb_geamax`|This function computes the infinity-norm of a vector x. If x is a real vector it computes infinity norm as: $$ amax \leftarrow max \mid x_i \mid $$ else if x is a complex vector then it computes the infinity-norm as: $$ amax \leftarrow max(\mid re(x_i) \mid + \mid im(x_i) \mid)$$|No ❌|
|**Generalized Infinity Norm**|`psb_geamaxs`|This subroutine computes a series of infinity norms on the columns of a dense matrix x: $$ res(i) \leftarrow max_k \mid x(k,i) \mid $$|        No ❌        |
| **1-Norm of Vector**| `psb_norm1` / `psb_geasums`|This function computes the 1-norm of a vector x. If x is a real vector it computes 1-norm as: $$ asum \leftarrow \mid \mid x_i \mid \mid $$ else if x is a complex vector then it computes 1-norm as: $$ asum \leftarrow \mid \mid re(x) \mid \mid_1 + \mid \mid im(x) \mid \mid_1 $$ |No ❌|
|**Generalized 1-Norm of Vector**|`psb_geasums`|This subroutine computes a series of 1-norms on the columns of a dense matrix x: $$ res(i) \leftarrow max_k \mid x(k,i) \mid $$ This function computes the 1-norm of a vector x. If x is a real vector it computes 1-norm as: $$ res(i) \leftarrow \mid \mid x_i \mid \mid $$ else if x is a complex vector then it computes 1-norm as: $$ res(i) \leftarrow \mid \mid re(x) \mid \mid_\ + \mid \mid im(x) \mid \mid_1 $$ |No ❌|
| **2-Norm of Vector**|`psb_norm2` / `psb_genrm2`| This function computes the 2-norm of a vector x. If x is a real vector it computes 2-norm as: $$ nrm2 \leftarrow \sqrt{x^T  x}$$ else if x is a complex vector then it computes 2-norm as: $$ nrm2 \leftarrow \sqrt{x^H  x}$$|No ❌|
|**Generalized 2-Norm of Vector**|`psb_genrm2s` / `psb_spnrm1` |This subroutine computes a series of 2-norms on the columns of a dense matrix x: $$ res(i) \leftarrow \mid \mid x(:,i) \mid \mid_2 $$|No ❌|
|**1-Norm of Sparse Matrix**|`psb_norm1`|This function computes the 1-norm of a matrix A: $$ nrm1 \leftarrow \mid \mid A \mid \mid_1 $$ where A represents the global matrix A|No ❌|
|**Infinity Norm of Sparse Matrix**|`psb_normi` / `psb_spnrmi`|This function computes the infinity-norm of a matrix A: $$ nrmi \leftarrow \mid \mid A \mid \mid_{\infty}$$ where: A represents the global matrix A|No ❌|
|**Sparse Matrix by Dense Matrix Product**| `psb_spmm`|This subroutine computes the Sparse Matrix by Dense Matrix Product: $$ y \leftarrow \alpha  A  x + \beta  y$$ $$ y \leftarrow \alpha  A^T  x + \beta  y$$ $$ y \leftarrow \alpha  A^H  x + \beta  y$$ where: <br> x is the global dense matrix x_{:,:} <br> y is the global dense matrix y_{:,:} <br> A is the global sparse matrix A|Work in progress ✅|
|**Triangular System Solve**|`psb_spsm`|This subroutine computes the Triangular System Solve: $$ y \leftarrow \alpha T^{-1}  x + \beta  y $$ $$ y \leftarrow \alpha D^{-1}  x + \beta  y $$ $$ y \leftarrow \alpha T^{-1}  D  x + \beta  y $$ $$ y \leftarrow \alpha T^{-T}  x + \beta  y $$ $$ y \leftarrow \alpha D  T^{-T}  x + \beta  y $$ $$ y \leftarrow \alpha T^{-T}  D  x + \beta  y $$ $$ y \leftarrow \alpha T^{-H}  x + \beta  y $$ $$ y \leftarrow \alpha D  T^{-H}  x + \beta  y $$ $$ y \leftarrow \alpha T^{-H}  D  x + \beta  y $$ where: <br> x is the global dense matrix x_{:,:} <br> y is the global dense matrix y_{:,:} <br> T is the global sparse block triangular submatrix T <br> D is the scaling diagonal matrix|No ❌|
|**Entrywise Product**|`psb_gemlt`|This function computes the entrywise product between two vectors x and y $$ dot \leftarrow x(i)y(i)$$|No ❌|
|**Entrywise Division**|`psb_gediv`|This function computes the entrywise division between two vectors x and y $$ div \leftarrow \frac{x(i)}{y(i)}$$|No ❌|
|**Entrywise Inversion**|`psb_geinv`|This function computes the entrywise inverse of a vector x and puts it into y $$ inv \leftarrow \frac{1}{x(i)}$$|No ❌|
