# Sampling correlation matrices - Google Summer of Code 2022

### **Contributor:** &emsp; Huu Phuoc Le &emsp; (Sorbonne Université)
### **Mentors:** &emsp; &emsp; &nbsp; Apostolos Chalkis, Elias Tsigaridas, Vissarion Fisikopoulos (GeomScale)
<br>

**<span style="color:red">This repository serves only as a summary for my work. All the codes produced during the coding phase are given in the links right below!</span>**

**Project proposal:** [https://www-polsys.lip6.fr/~phuoc/includes/gsoc22.pdf](https://www-polsys.lip6.fr/~phuoc/includes/gsoc22.pdf)
* Pull request 1: [https://github.com/GeomScale/volesti/pull/233](https://github.com/GeomScale/volesti/pull/233) <br>
Branch: [https://github.com/huuphuocle/volesti/tree/feature/sampling_correlation_matrices](https://github.com/huuphuocle/volesti/tree/feature/sampling_correlation_matrices)

* Pull request 2: [https://github.com/GeomScale/volesti/pull/240](https://github.com/GeomScale/volesti/pull/240)<br>
Branch: [https://github.com/huuphuocle/volesti/tree/feature/ReHMCWalk](https://github.com/huuphuocle/volesti/tree/feature/ReHMCWalk)

**Notations**
* $A$: a correlation matrix with size $n \times n$
* $A_{i,j}$: a matrix in the linear matrix inequality of $K$ where $a_{i,j} = a_{j,i} = 1$ and $0$ elsewhere
* $K$: the spectrahedron associated to the set of $n\times n$ correlation matrices
* `VT`: Vector Type `Eigen::Matrix<NT, Eigen::Dynamic, 1>`
* `MT`: Matrix Type `Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic>`

## I - Implementation summary & Usage

In this section, I summarize the new implementation to be integrated into `volesti` for sampling correlation matrices.

### **Two new main spectrahedron classes & ReHMC Walks**
I implemented two new classes for representing the spectrahedron $K$ that inherit `volesti`'s `Spectrahedron` class in `include/convex_bodies/correlation_matrices`:

* `CorrelationSpectrahedron` (`correlation_spectrahedron.hpp`)
* `CorrelationSpectrahedron_MT` (`correlation_spectrahedron_MT.hpp`) <br>
with `CorreMatrix` (`corre_matrix.hpp`) for encoding matrices as `Point`.

Some modifications are also made in existing classes of `volesti` to enable `BallWalk`, `RDHRWalk` and `BilliardWalk` on these new classes.

I implemented ReHMC Walk for sampling from log-concave distributions in `include/random_walks`:

* `GaussianReHMCWalk` (`gaussian_ReHMC_correlation.hpp`)
* `ExponentialReHMCWalk` (`exponential_ReHMC_correlation.hpp`)

I added several functions in `include/sampling/sample_correlation_matrices.hpp` to facilitate the use of sampling correlation matrices in `volesti`. These functions are listed below according to the sampling distribution and templated by

* `WalkType`: type of the random walk in use, e.g., BallWalk, RDHRWalk, BilliardWalk
* `PointType`: type of points in use (normal `Point` or `CorreMatrix`)
* `RNGType`: type of random number generators
* `PointList`: type of lists of sample points, e.g., `std::vector<PointType>` or `std::list<PointType>`.

**Common parameters**
* `n`: size of matrices
* `randPoints`: a list of `Point` to store output correlation matrices
* `walkL`: walk length of the random walks
* `num_points`: the requested number of ouput matrices
* `nburns`: the number of burnt samples

**Uniform sampling**
* `uniform_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, nburns)`: uniform sampling with the `CorrelationSpectrahedron` class
* `uniform_correlation_sampling_MT<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, nburns)`: uniform sampling with the `CorrelationSpectrahedron_MT` class



**Gaussian sampling:**
For a distribution of density proportional with
$$ \exp(-a \langle x, x \rangle)$$
where $a \in \mathbb{R}$

* `gaussian_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, a, nburns)`
* `gaussian_correlation_sampling_MT<WalkType, PointMT, RNGType>(n, randPoints, walkL, num_points, a, nburns)`

**Exponential sampling**

For a distribution of density proportional with
$$\exp(- \langle c, x \rangle/T)$$
where $c \in \mathbb{R}^d$ and $T \in \mathbb{R}$

* `exponential_correlation_sampling<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, c, T, nburns)`
* `exponential_correlation_sampling_MT<WalkType, Point, RNGType>(n, randPoints, walkL, num_points, c, T, nburns)`
### **Examples**

Please look at `examples/correlation_matrices/sampler.cpp` for a full example.

```
// Define some types

typedef BoostRandomNumberGenerator<boost::mt19937, double, 3>   RNGType;
typedef point<double>                                           Point;

// Declare the output list of points

std::vector<Point> randPoints;

// Set walk length

int walkL = 1;

// Uniform sampling of 1000 3x3 correlation matrices with `AcceleratedBilliardWalk`:

uniform_correlation_sampling<AcceleratedBilliardWalk, point<double>, RNGType>(3, randPoints, 1, num_points, 0);

// Gaussian sampling of 1000 3x3 correlation matrices with `GaussianReHMCWalk` using `CorreMatrix<double>` Point:

gaussian_correlation_sampling_MT<GaussianReHMCWalk, CorreMatrix<double>, RNGType>(3, randPoints, walkL, 1000, a, 0);
```

## II - More implementation details

The main operations for sampling correlation matrices are the oracles:

* **membership** tests whether the current point lies in the
* **intersection** computes the intersection of $x + t \cdot v$ and the boundary of $K$
* **reflection** computes the direction of reflection.

## CorrelationSpectrahedron

This class inherits `Spectrahedron` class. The main attributes of this class are

* `MT` : A correlation matrix $A = (a_{i,j})_{1 \leq i,j \leq n}$
* `VT` : A vector of independent coefficients of a correlation matrix
* `unsigned int n` : the size of matrices
* `unsigned int d` : the dimension of vectors of coefficients, $d = n(n-1)/2$.

In this class, the correlation matrices are still manipulated through vectors of coefficients with the help of `point` class and its basic operators. Therefore, integrating this new class `CorrelationSpectrahedron` into random walk algorithms existing in `volesti` is straightforward.

However, the important innovation here is that some member functions are modified to use new geometric oracles (for efficiency purposes). For instance, **membership** oracle `is_in` is changed to accelerate `BallWalk` algorithms, **intersection** (`line_positive_intersection`, `line_intersect`) and **reflection** (`compute_reflection`) oracles are changed to accelerate `BilliardWalk` algorithms.
These new aspects are detailed below in **Geometric oracles** section.

## CorrelationSpectrahedron_MT

This class serves for the same purposes of `CorrelationSpectrahedron` class. The main difference is that it uses a new PointType class `CorreMatrix` which stores a matrix in its attribute `mat`. This helps avoid switching between matrices and vectors of coefficients during computation. The implementation of this class includes the following:
* `CorreMatrix` class with basic operators (`+`,`-`,`*`,`/`,`+=`,`-=`, $\ldots$).
* `CorrelationSpectrahedron_MT` class with member functions modified to use `CorreMatrix` as PointType
* New `GetDirection` functions in `include/sampling/sphere.hpp` for `CorreMatrix` that returns a direction stored in matrix form.

Since we operate only on symmetric matrices and only the lower triangular part of the symmetric matrices are referenced by `Eigen::SelfAdjointEigenSolver`, c.f.,

[https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html](https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html),

the `CorreMatrix` class stores only the lower part of matrices and uses `triangularView<Eigen::Lower>` to operate on them.

## Geometric oracles

I added some functions in `include/matrix_operations/EigenvaluesProblems.h` to accelerate random walk algorithms:

* **Membership:** Testing $A \in K$ means testing whether $A$ is positive semidefinite. The old function in `volesti` computes the largest eigenvalue of $A$ with `Eigen::SelfAdjointEigenSolver` or `Spectra::SymEigsSolver`.<br><br>
I changed this to `Eigen::LDLT` that computes `A`'s Cholesky decomposition when `A` is positive semidefinite or signals an error otherwise (see `isPositiveSemidefinite` function).

* **Intersection:** The intersection oracle needs to solve a generalized eigenvalue problem for symmetric matrices. This principle is not changed but by experiments, I found that `Eigen::GeneralizedSelfAdjointEigenSolver` is faster than the currently in-use `Spectra::SymGEigsSolver`.<br><br>
Hence, I added a function `minPosLinearEigenvalue_EigenSymSolver` which is used by the intersection oracles in `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT` classes.

* **Reflection:** In the project proposal, we explained a formula for the normal vector at the reflection point that helps compute the reflection direction. However, each matrix $A_{i,j}$ in the linear matrix inequality of $K$ has only two non-zero entries $a_{i,j} = a_{j,i} = 1$. Substituting these matrices in the normal vector formula gives a much simpler formula:
$$n_{normal} \propto (e_2e_1,e_3e_2,e_3e_1,\ldots)$$
where $e = (e_1,\ldots,e_n)$ is an eigenvector obtained from **intersection** oracle. <br>
Using this formula is much faster than the naive one and helps improve the performance of random walks which intensively use reflections, e.g., `BilliardWalk` or `ReHMCWalk`.

## ReHMC Walk

As explained in the project proposal, Reflective Hamiltonian Monte-Carlo Walk allows to sample correlation matrices from a distribution over $\mathbb{R}^d$. Here, we suppose the distribution over $\mathbb{R}^d$ is log-concave and is defined by a density function
$$Z \cdot \exp(-f(x)).$$

In principle, the only difference between ReHMC walks comes from $f(x)$ which appears in the Hamiltonian for Metropolis filter
$$H(x,v) = f(x) + \frac{1}{2}\langle v, v \rangle$$
and how the direction $v$ is updated
$$v = v - \frac{\eta}{2}\nabla f(x).$$

I implemented ReHMC Walk for two log-concave distributions in `include/random_walks`:
* `GaussianReHMCWalk`: Gaussian distribution with log-density $- a \langle x, x \rangle$
* `ExponentialReHMCWalk`: Exponential distribution with log-density $-\langle c, x \rangle/T$.

Some small optimizations are done to skip repetitive computations (but the main speed-up is still coming from faster **reflection** oracles).

## III - Testing

Tests for my new implementations in `volesti` is written in `test/sampling_correlation_matrices_test.cpp` which consists of new tests using `doctest`. I also modified `test/CMakeLists.txt` to generate tests by `CMake`.

My tests assert whether
* the generated matrices are correlation matrices
* the Potential Scale Reduction Factor (PSRF) of the whole sample set is smaller than $1.1$.

I tested Ball Walk, RDHR Walk and Billiard Walk on `Spectrahedron`, `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT` for various matrix size $3 \leq n \leq 100$ and obtained:
* :heavy_check_mark: The matrices generated by random walk algorithms all satisfy to be correlation matrices
* The PSRF varies according to Walk Type, $n$ and the number of samples. The experiment section below will detail the difference of random walks in this aspect.


## IV - Experiments with random walk algorithms

In `examples/correlation_matrices/sampler.cpp`, I implemented the function:
* `old_uniform_sampling<WalkType, Point>(n, randPoints, walkL, num_points, nburns)`

that uses the general `Spectrahedron` class (and the LMI formed by $A_{i,j}$) to sample correlation matrices. This function is used in comparison with the new classes `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT`.

First, I will explain in details the speed-up for Ball Walk, RDHR Walk and Billiard Walk. These random walks already exist in `volesti` but the new operations in `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT` bring better performance.

### **Ball Walk algorithm**
Ball Walk is already implemented in `volesti`. To use this implementation with `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT`, we needed to change . Ball Walk algorithm depends mainly on two operations:
* Generating a random direction (see `GetDirection` functions in `include/sampling/sphere.hpp`).<br>
New `GetDirection` functions are added to return a direction in matrix form.
* Checking the positive semi-definiteness of a generated matrix.<br>
Here, the membership oracle `is_in` of `Spectrahedron` class and its subclasses is used intensively. As explained above, the old implementation of `is_in` is replaced by `Eigen::LDLT`. The table below illustrates the speed-up using this new `is_in` method.

<center>

Timings (in milliseconds) for Ball Walk to sample $100000$ $n\times n$ correlation matrices

| $n$ | Old implementation | CorrelationSpectrahedron | CorrelationSpectrahedron_MT |
|-----|--------------------|--------------------------|-----------------------------|
| 5   | 170                | 42                       | 51                          |
| 10  | 492                | 104                      | 145                         |
| 15  | 1069               | 206                      | 291                         |
| 20  | 2129               | 345                      | 490                         |
| 25  | 3981               | 525                      | 743                         |
| 30  | 6702               | 740                      | 1005                        |
| 35  | 11506              | 1003                     | 1405                        |
| 40  | 28597              | 1327                     | 1852                        |
| 45  | 51578              | 1675                     | 2349                        |
| 50  | 73294              | 2116                     | 2958                        |
| 55  | 113227             | 2554                     | 3593                        |
| 60  | 143268             | 3094                     | 4300                        |
| 65  | 186087             | 3658                     | 5073                        |
| 70  | 245637             | 4306                     | 5975                        |
| 75  | 317582             | 5017                     | 7109                        |
| 80  | 411608             | 5949                     | 8195                        |
| 85  | 510693             | 6806                     | 9317                        |

</center>

<figure align="center">
    <img src="img/ballwalk.png" />
    <figcaption>Fig.1 - 3000 uniform samples by Ball Walk.
</figure>

### **RDHR Walk algorithm**

RDHR Walk mainly uses `line_intersect` functions. I observed in experiments that the `line_intersect` function of `Spectrahedron` class is much slower than the ones of `CorrelationSpectrahedron` and `CorrelationSpectrahedron`. This bottleneck comes from the difference in implementations `line_intersect` in those classes:
* In `Spectrahedron`, it uses two eigenvalue problems that construct the matrices twice and compute respectively the positive and negative distances. Moreover, I corrected a bug which leads to wrong sign of negative distances (a `-` missing).
* In the new classes, I use `Spectra` to compute in one function call both distances.

Deeper reasons that lead to large difference in performance may depend on how `Spectra` is implemented and are unclear for me at the moment.

<center>

Timings (in milliseconds) for RDHR Walk to sample $10000$ $n\times n$ correlation matrices

| $n$ | Old implementation | CorrelationSpectrahedron | CorrelationSpectrahedron_MT |
|-----|--------------------|--------------------------|-----------------------------|
| 5   | 89                 | 29                       | 31                          |
| 10  | 190                | 195                      | 203                         |
| 15  | 467                | 262                      | 282                         |
| 20  | 790                | 312                      | 335                         |
| 25  | 1400               | 386                      | 415                         |
| 30  | 2248               | 448                      | 492                         |
| 35  | 3500               | 561                      | 634                         |
| 40  | 5500               | 624                      | 706                         |
| 45  | 9320               | 761                      | 860                         |
| 50  | 17150              | 829                      | 957                         |
| 55  | 33910              | 984                      | 1143                        |
| 60  | 42530              | 1070                     | 1242                        |
| 65  | 63800              | 1262                     | 1476                        |
| 70  | 88000              | 1348                     | 1596                        |

</center>

<figure align="center">
    <img src="img/rdhrwalk.png" />
    <figcaption>Fig.2 - 3000 uniform samples by Hit-and-Run Walk.
</figure>

### **Billiard Walk algorithm**

Billiard Walk mainly uses the `intersection` and `reflection` oracles. The new implementations of these oracles in `CorrelationSpectrahedron` and `CorrelationSpectrahedron_MT` accelerate Billiard Walk in `volesti` for sampling correlation matrices.

The main acceleration comes from the new `compute_reflection` function which is specialized for correlation matrices. Since Billiard Walk uses intensively **reflection**, this results in an impressive performance illustrated in the table below.

<center>

Timings (in milliseconds) for Billiard Walk to sample $1000$ $n\times n$ correlation matrices

| $n$ | Old implementation | CorrelationSpectrahedron | CorrelationSpectrahedron_MT |
|-----|--------------------|--------------------------|-----------------------------|
| 5   | 40                 | 0                        | 0                           |
| 10  | 630                | 70                       | 70                          |
| 15  | 4110               | 360                      | 350                         |
| 20  | 12650              | 920                      | 900                         |
| 25  | 37330              | 2150                     | 2280                        |
| 30  | 121150             | 4780                     | 4640                        |
| 35  | 350640             | 8660                     | 8790                        |
| 40  | 868060             | 14200                    | 14690                       |

</center>

<figure align="center">
    <img src="img/billiardwalk.png" />
    <figcaption>Fig.3 - 3000 uniform samples by Billiard Walk.
</figure>

### **ReHMC Walk algorithm**

<figure align="center">
    <img src="img/gaussianReHMC.png" />
    <figcaption>Fig.4 - 3000 Gaussian samples by ReHMC Walk.
    </figcaption>
</figure>

<figure align="center">
    <img src="img/exponentialReHMC.png" />
    <figcaption>Fig.5 - 3000 exponential samples by ReHMC Walk.
    </figcaption>
</figure>

## Synthesis of experiments

For Ball Walk, RDHR Walk and Billiard Walk, the new classes provide good performance for generating correlation matrices comparing with `Spectrahedron` class.

From the tables of timing above, we see that each step of Ball Walk and RDHR Walk is faster than Billiard Walk (which is obvious since each step of Billiard Walk is more sophisticated).

However, we also need to consider the mixing time, i.e., the number of steps to achieve the convergence of Markov random walks. This mixing time can be measured by PSRF which, by experiments, is much smaller for Billiard Walk. For instance, in high dimension, e.g., $n = 100$, we cannot obtain PSRF < 1.1 by Ball Walk or RDHR Walk even with $10000$ points while Billiard Walk reaches this PSRF with only $1000$ points.

## V - Further plan

A possible follow-up after the final submission of Google Summer of Code 2022 consists of:
* Reimplement **membership** and **intersection** oracles in `Spectrahedron` (e.g. using `Eigen::LDLT`)
* Comparing with other software for generating correlation matrices
* A general implementation of ReHMC walk that is templated by a functor `f(x)`
* Further optimization of member functions of `CorrelationSpectrahedron_MT` and `CorreMatrix` classes
* Simplify user functions (e.g., fusion `_MT` and non `MT` user functions)

### Doctest bug in `volesti`

While building `volesti`'s test on an Ubuntu 22.04 system, I encountered an error:

```
doctest.h:4189:47: error: size of array 'altStackMem' is not an integral constant-expression
 4189 |         static char             altStackMem[4 * SIGSTKSZ];
 ```

This comes from `doctest.h` v1.2.9 being used in `test/`. A proposed fix is to use `doctest.h` v2.4.8 found at [https://github.com/doctest/doctest/blob/v2.4.8/doctest/doctest.h](https://github.com/doctest/doctest/blob/v2.4.8/doctest/doctest.h).
