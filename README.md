# PSBLAS library, version 3.9


The PSBLAS library, developed with the aim to facilitate the parallelization of computationally intensive scientific applications, is designed to address parallel implementation of iterative solvers for sparse linear systems through the distributed memory paradigm. It includes routines for multiplying sparse matrices by dense matrices, solving block diagonal systems with triangular diagonal entries, preprocessing sparse matrices, and contains additional routines for dense matrix operations. The current implementation of PSBLAS addresses a distributed memory execution model operating with message passing.

The PSBLAS library version 3 is implemented in the Fortran 2008 programming language, with reuse and/or adaptation of existing Fortran 77 and Fortran 95 software, plus a handful of C routines. 

## References


The architecture, philosophy and implementation details of the library are contained in the following papers:

- The architecture of the Fortran 2003 sparse BLAS is described in:
  >S. Filippone, A. Buttari. Object-Oriented Techniques for Sparse Matrix
  >Computations in Fortran 2003, ACM Trans. on Math. Software, vol. 38, No.
  4, 2012.

- The ideas are explored further with the paper:
  >V. Cardellini, S. Filippone and D. Rouson. Design Patterns for
  >sparse-matrix computations on hybrid CPU/GPU platforms, Scientific 
  >Programming, 22(2014), pp.1-19.

- Version 1.0 of the library is described in:
  >S. Filippone, M. Colajanni. PSBLAS: A library for parallel linear
  >algebra computation on sparse matrices, ACM Trans. on Math. Software,
  >26(4), Dec. 2000, pp. 527-550.
- The software infrastructure changes required to accommodate the implementation of the
  Additive-Schwarz preconditioners available in [AMG4PSBLAS](https://github.com/sfilippone/amg4psblas/) are detailed in:
  > A. Buttari, P. D'Ambra, D. di Serafino, S. Filippone, Extending PSBLAS to build parallel Schwarz preconditioners, Applied Parallel Computing. State of the Art in Scientific Computing: 7th International Workshop, PARA 2004, LNCS 3732, 2006, pp. 593-602.
  
  > A. Buttari,  P. D'Ambra, D. Di Serafino, S. Filippone, 2LEV-D2P4: A package of high-performance preconditioners for scientific and engineering applications, Applicable Algebra in Engineering, Communications and Computing, 2007, 18(3), pp. 223-239.
  
  > P. D'Ambra, D. Di Serafino, S. Filippone, MLD2P4: A package of parallel algebraic multilevel domain decomposition preconditioners in Fortran 95 ACM Transactions on Mathematical Software, 2010, 37(3), 30

PSBLAS is the backbone of the Parallel Sparse Computation Toolkit ([PSCToolkit](https://psctoolkit.github.io/)) suite of libraries. See the paper:
> D’Ambra, P., Durastante, F., & Filippone, S. (2023). Parallel Sparse Computation Toolkit. Software Impacts, 15, 100463.

### Other Software credits 

We originally included a modified implementation of some of the Sparker
(serial sparse BLAS)  material; this has been completely rewritten, way
beyond the intention(s) and responsibilities of the original developers.
The main reference for the serial sparse BLAS is:
>Duff, I., Marrone, M., Radicati, G., and Vittoli, C. Level 3 basic 
>linear algebra subprograms for sparse matrices: a user level interface,
>ACM Trans. Math. Softw., 23(3), 379-401, 1997.

## Installing

To compile and run our software you will need the following
prerequisites (see also SERIAL below):

1. A working version of MPI

2. A version of the BLAS; if you don't have a specific version for your
   platform you may try ATLAS available from
   http://math-atlas.sourceforge.net/ 

3. We have had good results with  the METIS library, from 
   http://www-users.cs.umn.edu/~karypis/metis/metis/main.html.
   This is optional; it is  used in the util and test/fileread
   directories but only if you specify `--with-metis`.

4. If you have the AMD package of Davis, Duff and Amestoy, you can
   specify `--with-amd` (see `./configure --help` for more details).
   We use the C interface to AMD.

5. If you have CUDA available, use
   --enable-cuda           to compile CUDA-enabled methods
   --with-cudadir=<path>   to specify the CUDA toolkit location
   --with-cudacc=XX,YY,ZZ  to specify a list of target CCs (compute
   			   capabilities) to compile the CUDA code for.

The configure script will generate a Make.inc file suitable for building
the library. The script is capable of recognizing the needed libraries
with their default names; if they are in unusual places consider adding
the paths with `--with-libs`, or explicitly specifying the names in
`--with-blas`, etc. 

>[!CAUTION]
> Please note that a common way for the configure script
> to fail is to specify inconsistent MPI vs. plain compilers, either
> directly or indirectly via environment variables; e.g. specifying the
> Intel compiler with `FC=ifort` while at the same time having an 
> `MPIFC=mpif90` which points to GNU Fortran. 

>[!TIP]
> The best way to avoid this
> situation is (in our opinion) to use the environment modules package
> (see [http://modules.sourceforge.net/](http://modules.sourceforge.net/)), and load the relevant
> variables with (e.g.) 
> ```
> module load gcc/13.2.0 openmpi/4.1.6
> ```
> This will delegate to the modules setup to make sure that the version of
> openmpi in use is the one compiled with the gnu46 compilers. After the
> configure script has completed you can always tweak the Make.inc file
> yourself.

After you have Make.inc fixed,  run 
```
make
``` 
to  compile the library; go to the test directory and its subdirectories
to get test programs done. If you specify `--prefix=/path` you can do make
install and the libraries will be installed under `/path/lib`, while the
module files will be installed under `/path/modules`. The regular and
experimental C interface header files are under `/path/include`.

### CUDA and GPU support

This version of PSBLAS incorporates into a single package three
entities that were previouslty separated:
| Library |                    |
|---------|--------------------|
| PSBLAS  | the base library   |
| PSBLAS-EXT | a library providing additional storage formats for matrices and vectors |
| SPGPU      | a package of kernels for NVIDIA GPUs originally written by Davide Barbieri and Salvatore Filippone; see the license file [cuda/License-spgpu.md](cuda/License-spgpu.md) |

### OpenACC
There is a highly experimental version of an OpenACC interface,
you can access it by speficifying
```bash
--enable-openacc  --with-extraopenacc="-foffload=nvptx-none=-march=sm_70"
```
where the argument to the extraopenacc option depends on the compiler
you are using (the example shown here is relevant for the GNU
compiler). 

### Serial

Configuring with `--enable-serial` will provide a fake MPI stub library
that enables running in pure serial mode; no MPI installation is needed
in this case (but note that the fake MPI stubs are only guaranteed to
cover what we use internally, it's not a complete replacement). 

### Integers

We have two kind of integers: IPK for local indices, and LPK for
global indices. They can be specified independently at configure time,
e.g.
```bash
--with-ipk=4 --with-lpk=8
```
which is asking for 4-bytes local indices, and 8-bytes global indices
(this is the default). 


## Documentation

Further information on installation and configuration can be found in the documentation.
See [docs/psblas-3.9.pdf](docs/psblas-3.9.pdf); an HTML version of the same document is
available in docs/html. Please consult the sample programs, especially
- [test/pargen/psb_s_pde2d.F90](test/pargen/psb_s_pde2d.F90) [test/pargen/psb_d_pde2d.F90](test/pargen/psb_d_pde2d.F90)
- [test/pargen/psb_s_pde2d.F90](test/pargen/psb_s_pde3d.F90) [test/pargen/psb_d_pde2d.F90](test/pargen/psb_d_pde3d.F90)

which contain examples for the solution of linear systems obtained by the discretization of a generic second-order differential equation in two:
```math
- a_1 \frac{\partial^2 u}{\partial x^2} 
- a_2 \frac{\partial^2 u}{\partial y^2} 
+ b_1 \frac{\partial u}{\partial x} 
+ b_2 \frac{\partial u}{\partial y} 
+ c u = f
```
or three
```math
- a_1 \frac{\partial^2 u}{\partial x^2} 
- a_2 \frac{\partial^2 u}{\partial y^2} 
- a_3 \frac{\partial^2 u}{\partial z^2} 
+ b_1 \frac{\partial u}{\partial x} 
+ b_2 \frac{\partial u}{\partial y} 
+ b_3 \frac{\partial u}{\partial z} 
+ c u = f
```
dimensions on the unit square/cube with Dirichlet boundary conditions.

### Utilities

The [test/util](test/util) directory contains some utilities to convert to/from
Harwell-Boeing and MatrixMarket file formats.

## TODO and bugs

- [ ] Improving OpenACC support
- [ ] Improving OpenMP support
- [X] Fix all reamining bugs. Bugs? We dont' have any ! 🤓

> [!NOTE]
> To report bugs 🐛 or issues ❓ please use the [GitHub issue system](https://github.com/sfilippone/psblas3/issues).



## The PSBLAS team. 
**Project lead:**
Salvatore Filippone

**Contributors** (_roughly reverse cronological order_):

- Theophane  Loloum
- Fabio      Durastante
- Dimitri    Walther
- Andea      Di Iorio
- Stefano    Petrilli
- Soren 	   Rasmussen
- Zaak       Beekman
- Ambra	   Abdullahi Hassan
- Pasqua	   D'Ambra
- Alfredo    Buttari
- Daniela    di Serafino
- Michele    Martone
- Michele    Colajanni
- Fabio      Cerioni
- Stefano    Maiolatesi
- Dario      Pascucci


## RELATED SOFTWARE
If you are looking for more sophisticated preconditioners, you may be
interested in the package AMG4PSBLAS from
<http://github.com/sfilippone/amg4psblas> and the whole [PSCTooolkit suite](https://psctoolkit.github.io/).


Contact: <https://github.com/sfilippone/psblas3>
