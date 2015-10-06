# Sparse

Handle Sparse Matrices (and regular ones).

This library is used with lsqr (C. C. Paige and M. A. Saunders) to solve *very large* system of linear equations (see [lsqrsolve](http://github.com/marcopovitch/lsqrsolve)).


# Beware

Still usable but written 10 years ago !  Use it only if you know what you are doing ... no support will be provided !

# Compilation

To make all the autotools files, run :

`glibtoolize -i`

`./autogen.sh`

To configure and install the sparse library in your installation directory `$INSTALL_PATH`, run :

`.configure --prefix=$INSTALL_PATH`

`make install`
 
 

# Links

* [http://web.stanford.edu/group/SOL/software/lsqr/](http://web.stanford.edu/group/SOL/software/lsqr/)
* [ http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c](http://stanford.edu/group/SOL/software/lsqr/c/lsqr.c)
