# Tests for GOSC 2022 GeomScale: Sampling correlation matrices

## Description
All the implementations for the tests are given in the `C++` file `test_sampling/sampler.cpp`.

Some functions are taken and simplified from `volesti` (removing template, struct modified to normal function, etc) to make the test more
self-contained. The main calls to `volesti` are for `Spectrahedron` class and its methods `positiveLinearIntersection` 
and `compute_reflection`.

Test 3 is stilled in testing and benchmarking.
## To launch

* Copy the folder to `volesti/examples/`
* Build with cmake: `cmake . -DLP_SOLVE=/path/to/lpsolve/ && make`
* Launch `./sampler`