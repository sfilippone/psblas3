# Introduction
This is a directory developed by Luca Pepè Sciarria and Simone Staccone froma Tor Vergata University to start to create some unit tests for PSBLAS 3.9, in particular for psb_spmm routine.

## Environment
These tests are developed using a linux environment, in particular Rocky Linux 9.5 (Blue Onyx). 

The compiler used is:
- gnu 12.2.1

The necessary dependnces are:
- mpich 4.2.2
- psblas 3.9

In order to have the exact same environment used for testing compile PSBALS library using cuda 12.5.


## Getting started
Steps to reproduce the tests:
- make
- insert the matrix files inside the matrix/ directory (or create one if it doesn't exists; psblas3/test/spmm/matrix/)
- run ./runs/psb_spmm_test
- ... 

## Test goal
Check the correctness of the matrix-vector multiplication $y = Ax$ using the **psb_spmm** routine, checking for all the test suite cases.

## Test Suite
**A** matrixes are located in the matrix/ directory
|Matrix|File Name|Type|
|:-:|:-:|:-:|
|$A_1$|1138_bus.mtx|Sparse| 

**x** vectors are located in the vectors/ directory. They are generated randomly using the same seed and then saved on different files based on their characteristics. The size of the vector is choosen accordingly to the size of the matrix column space considered for the single test instance.
|Vector|File Name|Coefficients|Coefficients Description|
|:-:|:-:|:-:|:-:|
|$x_1$|x1.txt|$x_i> 0, \forall i$|Positive coefficients| 
|$x_2$|x2.txt|$x_i < 0, \forall i$|Negative coefficients
|$x_3$|x3.txt|$x_i \ne 0, \forall i$|Random coefficients
|$x_4$|x4.txt|$x_i = 0, \forall i$|Null coefficients

## Output
The results of the computation will be saved on different files based on the instance of the test considered. In particular the naming conventiona format the output file as sol_m#_x#_y#.mtx, where each # is a number choosen w.r.t. the test instance. (Ex. sol_m1_x1_y1.mtx is the solution computed using the first matrix file, the first x vector file and the first y vector file).

## Notes
For now only integer multiplication is tested and on a single matrix