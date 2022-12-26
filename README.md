# smallexpimv
Fortran routines and a python wrapper to apply the imaginary exponential of a matrix (real tridiagonal or complex). More precisely,

$$
\beta\mathrm{e}^{\pm\mathrm{i} t X} u
$$

for a matrix $X\in\mathbb{C}^{n\times n}$, a time-step $t\in\mathbb{R}$, a scaling factor $\beta\in\mathbb{R}$, a vector $u\in\mathbb{C}^n$ and the imaginary number $\mathrm{i}^2 = -1 $. This matrix exponential is approximated by an adaptive and restarted Taylor approximation. We also consider a tridiagonal case, namely

$$
\beta\mathrm{e}^{\pm\mathrm{i} t T} e_1
$$

for a tridiagonal symmetric matrix $T\in\mathbb{R}^{n\times n}$ and the first unit vector $e_1 = (1,0,\ldots,0)^\ast \in\mathbb{R}^{n}$.  For the tridiagonal case we use a lapack eigendecomposition.

The matrix exponentials above have some relevance when applying Lanczos or Arnoldi propagators of Schrödinger-type problems. In this context, $X$ or $T$ correspond to small dimensional Krylov representations of a large problem.

The main code is written in Fortran, and can be used in Python using f2py.

To generate the main python module smallexpimv.so run
> python3 -m numpy.f2py -c src/F_smallexpimv.F90 -m smallexpimv -lblas

optional, generate f2py header files:
> python3 -m numpy.f2py -h smexp.pyf -m smallexpimv src/F_smallexpimv.F90

To handle working memory in Python (most relevant when applying the matrix exponential multiple times)
we also provide the module smallexpimv_pyclass.py. Keep in mind that using this module requires smallexpimv.so

Makefile: use
> make test

to run python and fortran tests and
> make libs

to generate smallexpimv.so and smallexpimv_pyclass.pyc. To use the Makefile, first check if the compiler commands and the f2py filename extension are defined correctly therein for your system!
