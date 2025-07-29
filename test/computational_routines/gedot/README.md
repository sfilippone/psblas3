# Introduction
This is a directory developed by Luca Pepè Sciarria and Simone Staccone from Tor Vergata University to start to create some unit tests for PSBLAS 3.9, in particular for ```psb_gedot``` routine.


## Getting started
Steps to reproduce the tests:

- Compile the code using ``` make ``` (Optional)
- Launch the script ```./autotest.sh``` or with ```source ./autotest.sh``` if you want to add modules to the .bashrc file permanently.
- Check the output log file ```psb_gedot_test.log``` to collect results and check for errors. In any case a summarization of tests passed should be shown in stdout.

NOTE: If the code is changed and a new compilation is needed to show the changes, the ```autotest.sh``` script isn't aware of this scenario, therefore it is necessary to manually recompile the code. Moreover, if you are running manually the script and not launching the main ```test.sh``` script, be careful to use the last compiled version of the utils.

## Test Suite
The ```psb_gedot``` is a function that performs the scalar product between two vectors giving a value as result. The signature of the function is:

```fortran
psb_gedot(x, y, desc_a, info [,global])
```

The strategy to validate the correctness of the computation is to compare single precision result and double precision result in the test cases in which the test should not give an error. In this way it is possible to have a correctness check of the computation comparing the two results considering a number of significant digits which is tuned on the single precision computation.

### Parameters Values
**x** vectors are located in the vectors/ directory. They are generated randomly using the same seed and then saved on different files based on their characteristics. The size of the vector is choosen accordingly to the size of the matrix column space considered for the single test instance.
|Variable|File Name|Coefficients|Coefficients Description|
|:-:|:-:|:-:|:-:|
|$x_1$|x1.txt|$x_i> 0, \forall i$|Positive coefficients| 
|$x_2$|x2.txt|$x_i < 0, \forall i$|Negative coefficients
|$x_3$|x3.txt|$x_i \ne 0, \forall i$|Random coefficients
|$x_4$|x4.txt|$x_i = 0, \forall i$|Null coefficients

**y** vectors are located in the vectors/ directory. They are generated randomly using the same seed and then saved on different files based on their characteristics. The size of the vector is choosen accordingly to the size of the matrix rows space considered for the single test instance.
|Variable|File Name|Coefficients|Coefficients Description|
|:-:|:-:|:-:|:-:|
|$y_1$|y1.txt|$y_i> 0, \forall i$|Positive coefficients| 
|$y_2$|y2.txt|$y_i < 0, \forall i$|Negative coefficients
|$y_3$|y3.txt|$y_i \ne 0, \forall i$|Random coefficients
|$y_4$|y4.txt|$y_i = 0, \forall i$|Null coefficients

**global** logical value indicating if the computation is global among all processes or it is only local, so that a global reduction is needed afterwards.
|Variable|Value|Value Description|
|:-:|:-:|:-:|
|$global_1$|True|Global computation is required implicitly| 
|$global_2$|False|Global computation is required explicitly, so only a local computation is done|


## TODO
List of things still to add in the test:
- Add computation with broken descriptor and catch the errore result
- Test using complex data ($dot \leftarrow x^H \cdot y$)
- Test also GPU excecution
- Try multiple distributions