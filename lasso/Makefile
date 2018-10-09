include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in

FFTW3_INC = -I/soft/fftw/3.3.4/include
FFTW3_LIB = -L/soft/fftw/3.3.4/lib -lfftw3

INCLUDE := $(MIROOTLIB_INC) $(MIBLASLIB_INC) $(MILIB_INC) $(ROOT_INC) $(FFTW3_INC) $(MT_INC)
LDLIBS  := $(MIROOTLIB_LIB) $(MIBLASLIB_LIB) $(MILIB_LIB) $(ROOT_ALLLIBS) $(FFTW3_LIB) $(BLAS_LIB) $(MT_LIB) $(LIB) -lgfortran

TARGET = lasso

all: $(TARGET) 

lasso: lasso.o  arg_lasso.o  sub.o
	$(CXXN) $@.o arg_$@.o sub.o ${LDLIBS} -o $@

.SUFFIXES:
.SUFFIXES: .o .c .f .F .cc

.c.o:
	$(CCN) -c $<

.F.o : 
	$(FCN) -c $<

.f.o:
	$(FCN) -c $<

.cc.o:
	${CXXN} -c $<

clean:
	-rm -f *.o

cleaner: clean
	-rm -f ${TARGET} *% *~ core
