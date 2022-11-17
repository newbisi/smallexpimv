# python_smallexpimv
A Python class to apply the imaginary exponential of a matrix (real tridiagonal or complex) using Fortran code; these functions apply well to solve the small matrix exponential which occurs for Lanczos or Arnoldi propagator of Schrödinger-type problems 

compile python code:
> python3 -m numpy.f2py -c F_smallexpimv.F90 -m smallexpimv -lblas

optional, generate f2py header files:
> python3 -m numpy.f2py -h smexp.pyf -m smallexpimv F_smallexpimv.F90

run example:
> python3 example.py
