lib_LTLIBRARIES = libsparse.la

libsparse_la_SOURCES = \
	matrice.h matrice.c \
	sparse.h sparse.c

LIBRARY_VERSION=0:1:0
libsparse_la_LDFLAGS= -version-info $(LIBRARY_VERSION)

library_includedir=$(includedir)/sparse
library_include_HEADERS = matrice.h sparse.h 
