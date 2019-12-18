PSBLAS library, version 3.6
===========================

The architecture of the Fortran 2003 sparse BLAS is described in:
>S. Filippone, A. Buttari. Object-Oriented Techniques for Sparse Matrix
>Computations in Fortran 2003, ACM Trans. on Math. Software, vol. 38, No.
4, 2012.

The ideas are explored further with the paper:
>V. Cardellini, S. Filippone and D. Rouson. Design Patterns for
>sparse-matrix computations on hybrid CPU/GPU platforms, Scientific 
>Programming, 22(2014), pp.1-19.

Version 1.0 of the library is described in:
>S. Filippone, M. Colajanni. PSBLAS: A library for parallel linear
>algebra computation on sparse matrices, ACM Trans. on Math. Software,
>26(4), Dec. 2000, pp. 527-550.


UTILITIES
---------
The `test/util` directory contains some utilities to convert to/from
Harwell-Boeing and MatrixMarket file formats.


DOCUMENTATION
-------------
See docs/psblas-3.5.pdf; an HTML version of the same document is
available in docs/html. Please consult the sample programs, especially
test/pargen/psb_[sd]_pde[23]d.f90


OTHER SOFTWARE CREDITS 
----------------------
We originally included a modified implementation of some of the Sparker
(serial sparse BLAS)  material; this has been completely rewritten, way
beyond the intention(s) and responsibilities of the original developers.
The main reference for the serial sparse BLAS is:
>Duff, I., Marrone, M., Radicati, G., and Vittoli, C. Level 3 basic 
>linear algebra subprograms for sparse matrices: a user level interface,
>ACM Trans. Math. Softw., 23(3), 379-401, 1997.


INSTALLING
----------
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

The configure script will generate a Make.inc file suitable for building
the library. The script is capable of recognizing the needed libraries
with their default names; if they are in unusual places consider adding
the paths with `--with-libs`, or explicitly specifying the names in
`--with-blas`, etc. Please note that a common way for the configure script
to fail is to specify inconsistent MPI vs. plain compilers, either
directly or indirectly via environment variables; e.g. specifying the
Intel compiler with `FC=ifort` while at the same time having an 
`MPIFC=mpif90` which points to GNU Fortran. The best way to avoid this
situation is (in our opinion) to use the environment modules package
(see http://modules.sourceforge.net/), and load the relevant
variables with (e.g.) 
```
module load gnu46 openmpi
```
This will delegate to the modules setup to make sure that the version of
openmpi in use is the one compiled with the gnu46 compilers. After the
configure script has completed you can always tweak the Make.inc file
yourself. 

After you have Make.inc fixed,  run 
```
make
``` 
to  compile the library; go to the test directory and its subdirectories
to get test programs done. If you specify `--prefix=/path` you can do make
install and the libraries will be installed under `/path/lib`, while the
module files will be installed under `/path/modules`. The regular and
experimental C interface header files are under `/path/include`.

SERIAL
------
Configuring with `--enable-serial` will provide a fake MPI stub library
that enables running in pure serial mode; no MPI installation is needed
in this case (but note that the fake MPI stubs are only guaranteed to
cover what we use internally, it's not a complete replacement). 

LONG INTEGERS
-------------
We have an experimental flag `--enable-long-integers` that will enable 
having 8-byte integer data, allowing an index space larger than 2G; some
small cases have been tested but we do not offer full guarantee (yet).


TODO
----
Fix all reamining bugs. Bugs? We dont' have any ! ;-)


The PSBLAS team. 
---------------
Project lead:
Salvatore Filippone

Contributors (roughly reverse cronological order):

Soren 	   Rasmussen
Zaak       Beekman
Ambra	   Abdullahi Hassan
Pasqua	   D'Ambra
Alfredo    Buttari
Daniela    di Serafino
Michele    Martone
Michele    Colajanni
Fabio      Cerioni
Stefano    Maiolatesi
Dario      Pascucci



RELATED SOFTWARE
----------------
If you are looking for more sophisticated preconditioners, you may be
interested in the package MLD2P4 from
<http://github.com/sfilippone/mld2p4-2>


Contact: <https://github.com/sfilippone/psblas3>
