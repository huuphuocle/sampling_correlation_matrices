# Sampling correlation matrices - Google Summer of Code 2022

### **Contributor:** &emsp; Huu Phuoc Le &emsp; (Sorbonne Université)
### **Mentors:** &emsp; &emsp; &nbsp; Apostolos Chalkis &emsp; (GeomScale)
<br><br>

**Work repositories**
*  Fork branches: [https://github.com/huuphuocle/volesti/tree/feature/sampling_correlation_matrices](https://github.com/huuphuocle/volesti/tree/feature/sampling_correlation_matrices)
* Pull requests: [https://github.com/GeomScale/volesti/pull/233](https://github.com/GeomScale/volesti/pull/233)
* Project proposal: 

**Mathematical notations**
* $A$: a correlation matrix
* $n$: the size of $A$
* $K$: the spectrahedron associated to the set of $n\times n$ correlation matrices
* $A_{i,j}$: a matrix in the linear matrix inequality of $K$ where $a_{i,j} = a_{j,i} = 1$ and $0$ elsewhere

**Data Type:**
* `VT`: Vector Type `Eigen::Matrix<NT, Eigen::Dynamic, 1>`
* `MT`: Matrix Type `Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>`

## I - Implementation summary & Usage

In this section, I summarize the new implementation to be integrated into `volesti` for sampling correlation matrices.

### **Two new main spectrahedron classes & ReHMC Walks**
I implemented two new classes for representing the spectrahedron $K$ that inherit `volesti`'s `Spectrahedron` class in `include/convex_bodies/correlation_matrices`:  

* `CorreSpectra` (`corre_spectra.hpp`)
* `CorreSpectra_MT` (`corre_spectra_MT.hpp`)

and ReHMC Walk for sampling from log-concave distributions:

*
*

I added several functions in `include/sampling/sample_correlation_matrices.hpp` to facilitate the use of sampling correlation matrices in `volesti`. These functions are template function with template 

`template<typename WalkType, typename PointType, typename RNGType, typename PointList>`

where

* `WalkType`: type of the random walk in use, e.g., BallWalk, RDHRWalk, BilliardWalk
* `PointType`: type of points in use
* `RNGType`: type of random number generators
* `PointList`: type of lists of sample points, e.g., `std::vector<PointType>` or `std::list<PointType>`

These user functions are listed below, classified according to the sampling distribution:

**Uniform sampling** 
* `direct_uniform_sampling<...>(...)`: uniform sampling with the general `Spectrahedron` class (already in `volesti`) (used for comparison purposes)
* `uniform_correlation_sampling<...>(...)`: uniform sampling with the `CorreSpectra` class
* `uniform_correlation_sampling_MT<...>(...)`: uniform sampling with the `CorreSpectra_MT` class

**Exponential sampling**

**Gaussian sampling**

### **Examples**

Uniform sampling of `num_points` `n * n` correlation matrices 
`uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, 0)`

## II - More implementation details

The main operations for sampling correlation matrices are the oracles:

* `membership` tests whether the current point lies in the 
* `intersection` computes the intersection of $\bm{x} + t \cdot \bm{v}$ and the boundary of $K$
* `reflection` computes the direction of reflection

### CorreSpectra

This class inherits: The main instance variables of this class are

* `MT` : A correlation matrix $A = (a_{i,j})_{1 \leq i,j \leq n}$
* `VT` : A vector of independent coefficients of a correlation matrix
* `unsigned int d` :
* `unsigned int n` :

The membership oracle is changed to 

### CorreSpectra_MT

This class serves as the same but we want to experiment the idea of storing the whole matrix . We need to say that 

For this class, I needed to implement also a new method GetDirection (found in):

Based on the use of PointType template variables in random walk implementations of `volesti`, to handle MT PointType, I implemented a `CorreMatrix` class that operates as and its basic operators.

The old GetDirection method return a Vector Type (VT) . My new method returns a Matrix Type (MT). The way this method being implemented is .

### Geometric oracles

I added some functions in `include/matrix_operations/EigenvaluesProblems.h` to accelerate random walk algorithms:

* **Membership:** Testing $A \in K$ means testing whether $A$ is positive semidefinite. The old function in `volesti` computes the largest eigenvalue of $A$ with `Eigen::SelfAdjointEigenSolver` or `Spectra::SymEigsSolver`.<br><br> 
I changed this to `Eigen::LDLT` that computes `A`'s Cholesky decomposition when `A` is positive semidefinite or signals an error otherwise (see `isPositiveSemidefinite` function).

* **Intersection:** The intersection oracle needs to solve a generalized eigenvalue problem for symmetric matrices. This principle is not changed but by experiments, I found that `Eigen::GeneralizedSelfAdjointEigenSolver` is faster than the currently in-use `Spectra::SymGEigsSolver`.<br><br>
Hence, I added a function `minPosLinearEigenvalue_EigenSymSolver` which is used by the intersection oracles in `CorreSpectra` and `CorreSpectra_MT` classes.  

* **Reflection:** Each matrix $A_{i,j}$ in the linear matrix inequality of $K$ has only two non-zero entries $a_{i,j} = a_{j,i} = 1$. Substituting these matrices in the formula for the reflection direction we obtain: $$ a.$$
This new implementation is much faster than the naive one and helps improve the performance of random walks using reflection. 

### ReHMC Walk

## Testing

Testing in `volesti` is done in the `test/` folder which uses `doctest` (v. ) and `CMake` to generate tests.

I wrote ` test/sampling_correlation_matrices_test.cpp` including test functions and several `doctest`-format tests and added the aforementioned new tests in `test/CMakeLists.txt` for building. 

My tests assert whether
* the generated matrices are correlation matrices 
* the PRSF of the whole sample set is smaller than $1.1$.

I tested Ball Walk, RDHR Walk and Billiard Walk on `Spectrahedron`, `CorreSpectra` and `CorreSpectra_MT` for various matrix size $3 \leq n \leq 100$ and obtained:
* :heavy_check_mark: The matrices generated by random walk algorithms all satisfy to be correlation matrices
* The PRSF value varies following Walk Type, $n$ and the number of samples.


## Experiments with random walk algorithms

Note that 

First, I will explain in details the implementation changes in Ball Walk, RDHR Walk and Billiard Walk. These random walks already exist in `volesti` but the new operations in `CorreSpectra` and `CorreSpectra_MT` bring better performance.

### Ball Walk algorithm
Ball Walk is already implemented in `volesti`. To use this implementation to run with ampling correlation matrices, we needed to change . Ball Walk algorithm is done by 
*
*

Hence, the membership oracle `is_in` is used intensively. The old implementation of uses . My new `is_in` method uses . The table below illustrates the speed-up using this new `is_in` method. 

<center>

Timings (in milliseconds) for Ball Walk to sample $100000$ $n\times n$ correlation matrices

| $n$ | Old implementation | CorreSpectra | CorreSpectra_MT |
|-----|--------------------|--------------|-----------------|
| 5   | 170                | 42           | 51              |
| 10  | 492                | 104          | 145             |
| 15  | 1069               | 206          | 291             |
| 20  | 2129               | 345          | 490             |
| 25  | 3981               | 525          | 743             |
| 30  | 6702               | 740          | 1005            |
| 35  | 11506              | 1003         | 1405            |
| 40  | 28597              | 1327         | 1852            |
| 45  | 51578              | 1675         | 2349            |
| 50  | 73294              | 2116         | 2958            |
| 55  | 113227             | 2554         | 3593            |
| 60  | 143268             | 3094         | 4300            |
| 65  | 186087             | 3658         | 5073            |
| 70  | 245637             | 4306         | 5975            |
| 75  | 317582             | 5017         | 7109            |
| 80  | 411608             | 5949         | 8195            |
| 85  | 510693             | 6806         | 9317            |

</center>

### RDHR Walk algorithm

<center>

Timings (in milliseconds) for RDHR Walk to sample $10000$ $n\times n$ correlation matrices

| $n$ | Old implemetation | CorreSpectra | CorreSpectra_MT |
|-----|-------------------|--------------|-----------------|
| 5   | 66                | 29           | 31              |
| 10  | 180               | 195          | 203             |
| 15  | 427               | 262          | 282             |
| 20  | 650               | 312          | 335             |
| 25  | 1104              | 386          | 415             |
| 30  | 1787              | 448          | 492             |
| 35  | 3321              | 561          | 634             |
| 40  | 8301              | 624          | 706             |
| 45  | 16750             | 761          | 860             |
| 50  | 23656             | 829          | 957             |
| 55  | 39411             | 984          | 1143            |
| 60  | 46494             | 1063         | 1242            |
| 65  | 61584             | 1262         | 1476            |
| 70  | 85654             | 1348         | 1596            |
| 75  | 106953            | 1550         | 1853            |
| 80  | 140273            | 1712         | 1990            |
| 85  | 179069            | 1918         | 2278            |

</center>

### Billiard Walk algorithm

Billiard Walk mainly uses the `intersection` and `reflection` oracles. The new implementations of these oracles in `CorreSpectra` and `CorreSpectra_MT` accelerate Billiard Walk in `volesti` for sampling correlation matrices.

The main acceleration come from the new `reflection` which is specialized for correlation matrices. Since Billiard Walk is a random walk that uses intensively `reflection`, this results in an impressive performance illustrated in the table below.

<center>

Timings (in milliseconds) for Billiard Walk to sample $1000$ $n\times n$ correlation matrices

| $n$ | Old implementation | CorreSpectra | CorreSpectra_MT |
|-----|--------------------|--------------|-----------------|
| 5   | 40                 | 0            | 0               |
| 10  | 630                | 70           | 70              |
| 15  | 4110               | 360          | 350             |
| 20  | 12650              | 920          | 900             |
| 25  | 37330              | 2150         | 2280            |
| 30  | 121150             | 4780         | 4640            |
| 35  | 350640             | 8660         | 8790            |
| 40  | 868060             | 14200        | 14690           |

</center>

### Doctest bug in `volesti`

While building `volesti`'s test on an Ubuntu 22.04 system, I encountered an error: . This comes from `doctest.h` v1.2.9 being used in `test/`. A proposed fix is to use `doctest.h` v2.4.8 found at [https://github.com/doctest/doctest/blob/v2.4.8/doctest/doctest.h](https://github.com/doctest/doctest/blob/v2.4.8/doctest/doctest.h).

## IV - Remaining tasks