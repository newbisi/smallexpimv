#PYCENDING = cpython-36

.PHONY: clean

%.pyc: %.py
	python3 -m compileall $*.py
	mv ./__pycache__/$*.$(PYCENDING).pyc $*.pyc
	rmdir ./__pycache__

clean:
	rm -f *.o *.mod *.so *.pyc *.out
