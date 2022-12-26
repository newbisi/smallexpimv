# define f2py command and filename extension for .pyc and .so
F2PY = python3 -m numpy.f2py
PYCENDING = cpython-36
F2PYENDING = cpython-36m-x86_64-linux-gnu

# fortran compiler is only needed for fortran example
CFORTRAN = gfortran

PATHTOFORTRAN = ./src
PATHTOPYTHON = ./src
LIBS = -llapack -lblas

SOURCEF = $(PATHTOFORTRAN)/F_smallexpimv.F90
OBJECTS = F_smallexpimv.o

.PHONY: libs clean test ftest

%.so: $(PATHTOFORTRAN)/F_%.F90
	$(F2PY) -c $< -m $* $(LIBS)
	mv $*.$(F2PYENDING).so $*.so

%.out: %.o
	$(CFORTRAN) $(OBJECTS) $< -o $@ $(LIBS)

%.o: %.f90 $(SOURCEF)
	$(CFORTRAN) -c $(SOURCEF) $<

libs: smallexpimv.so $(PATHTOPYTHON)/smallexpimv_pyclasses.py
	$(MAKE) -C $(PATHTOPYTHON) smallexpimv_pyclasses.pyc PYCENDING=$(PYCENDING)
	mv $(PATHTOPYTHON)/smallexpimv_pyclasses.pyc ./smallexpimv_pyclasses.pyc

clean:
	rm -f *.o *.mod *.so *.pyc *.out *.pyf
	$(MAKE) -C $(PATHTOPYTHON) clean

test: smallexpimv.so $(PATHTOPYTHON)/smallexpimv_pyclasses.py example.py example.out
	$(info --- run Python and Fortran tests ---)
	$(MAKE) -C $(PATHTOPYTHON) smallexpimv_pyclasses.pyc PYCENDING=$(PYCENDING)
	mv $(PATHTOPYTHON)/smallexpimv_pyclasses.pyc ./smallexpimv_pyclasses.pyc
	python3 example.py
	./example.out

#	(cd src; ${MAKE} smallexpimv_pyclasses.pyc);
