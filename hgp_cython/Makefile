pyinc = -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/ -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/

hgp_cython.o: hgp_cython.c
	gcc -dynamic -pipe -fPIC -fwrapv -fno-common -fno-strict-aliasing $(pyinc) -c $^ -o $@

hgp_cython.so: hgp_cython.o
	gcc -isysroot / -L/opt/local/lib -bundle -undefined dynamic_lookup -L/opt/local/lib $^ -o $@

hgp_cython.c: hgp_cython.pyx
	cython $^
