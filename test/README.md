# Tests for GOSC 2022 GeomScale: Sampling correlation matrices

### These are tests required to qualify the phase of project proposals for GeomScale and are not the work done during the coding phase of Google Summer of Code 2022!

## Description
All the implementations for the tests are given in the `C++` file `test_sampling/sampler.cpp`.

Some functions are taken and simplified from `volesti` (removing template, struct modified to normal function, etc) to make the test more
self-contained. The main calls to `volesti` are for `Spectrahedron` class and its methods `positiveLinearIntersection`
and `compute_reflection`.

## To launch

* Copy the folder to `volesti/examples/`
* Build with cmake: `cmake . -DLP_SOLVE=/path/to/lpsolve/ && make`
* Launch `./sampler`
