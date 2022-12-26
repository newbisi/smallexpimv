# smallexpimv
A Python class to apply the imaginary exponential of a matrix (real tridiagonal or complex) using Fortran code; these functions apply well to solve the small matrix exponential which occurs for Lanczos or Arnoldi propagator of Schrödinger-type problems 

generate smallexpimv.so using f2py which can be imported as a module in python3:
> python3 -m numpy.f2py -c src/F_smallexpimv.F90 -m smallexpimv -lblas

optional, generate f2py header files:
> python3 -m numpy.f2py -h smexp.pyf -m smallexpimv src/F_smallexpimv.F90

to handle working memory for the fortran code (most relevant when applying the matrix exponential multiple times)
we provide the code smallexpimv_pyclass.py. Keep in mind that using this module requires smallexpimv.so

Makefile: use
> make test

to run python and fortran tests and
> make libs

to generate smallexpimv.so and smallexpimv_pyclass.pyc.
To run the makefile, first define the compiler commands and the f2py filename extension!
