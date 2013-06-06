
include ./Make.inc

# ----------------------------------------------------
# set CXX_FLAGS for release / debug
# ----------------------------------------------------
ifeq ($(USE_DEBUG), 1)
	CXX_FLAGS = $(CXX_DEBUG) $(CXX_WARNINGS)
else
	CXX_FLAGS = $(CXX_RELEASE) $(CXX_WARNINGS)
endif


# ----------------------------------------------------
# check whether we link static or dynamically
# ----------------------------------------------------
ifeq ($(LINK_SHARED), 1)
	CXX_FLAGS += -fpic
	F77_FLAGS += -fpic
	CXX_LIBS += $(CXX_LIBS_SHARED)
	TDKPSHELL = libtdkp.so tdkpshell 
	ifeq ($(USE_REMOTE), 0)
		CXX_INC += -DNOREMOTESOLVER
	else
		TDKPSHELL += atze
	endif
else
	CXX_LIBS += $(CXX_LIBS_STATIC)
	TDKPSHELL = tdkpshell-static
	ifeq ($(USE_REMOTE), 0)
		CXX_INC += -DNOREMOTESOLVER
	else
		TDKPSHELL += atze-static
	endif
endif

# --------------------------------------------------
# all objects which are required to link atze statically
# --------------------------------------------------
ATZEOBJS = bin/SocketIO.o bin/all.o bin/Logger.o \
	bin/ArpackSIPardisoSolver.o bin/CSRMatrix.o bin/CSRSparseMatrixInterface.o \
	bin/EigenSolver.o bin/Configuration.o 

all:  bin/Version.h $(TDKPSHELL)

depend: bin/Version.h
	rm -vf dependencies
	./make_depend.sh dependencies tdkp/common tdkp/geometry tdkp/io tdkp/kpmatrices tdkp/main tdkp/povray tdkp/solvers tdkp/utilities tdkp/probdefs tdkp/coulomb tdkp/clc 

test:
	(cd unittests; pwd; make)

clean:
	rm -f bin/*.o bin/*.cpp bin/*.h test.bin rm utilities/*

atze: tdkp/atze/atze_remote_gevp.cpp tdkp/utilities/SocketIO.h libtdkp.so
	$(CXX) $(CXX_FLAGS) -o utilities/atze_remote_gevp $^ $(CXX_INC) -L. -ltdkp $(CXX_LIBS)

atze-static: tdkp/atze/atze_remote_gevp.cpp tdkp/utilities/SocketIO.h tdkp/common/all.h \
	tdkp/solvers/ArpackSIPardisoSolver.h tdkp/main/CSRMatrix.h libtdkp.a
	$(CXX) $(CXX_FLAGS) -o utilities/atze_remote_gevp $^ $(CXX_INC) libtdkp.a -L. $(CXX_LIBS)

dox:
	doxygen doc/tdkpdoxyfile

utilities: bin/write_mat_to_matlab_readable

include ./dependencies

libtdkp.so: bin/ModifiedLapackRoutines.o $(OBJS)
	$(CXX) -shared -Wl,-soname,libtdkp.so -o libtdkp.so \
    $(OBJS) bin/ModifiedLapackRoutines.o $(CXX_LIBS)

libtdkp.a: bin/ModifiedLapackRoutines.o $(OBJS)
	ar rcs $@ $^

tdkpshell: libtdkp.so bin/tdkpshell.o bin/tdkpshell_wrap.o
	$(CXX) $(CXX_FLAGS) $(CXX_INC)	-o tdkpshell.bin \
	bin/tdkpshell.o \
	bin/tdkpshell_wrap.o \
	-Wl,-rpath,$(PWD) -L$(PWD) -ltdkp $(CXX_LIBS)

# -------------------------------------
# don't use shared lib
# -------------------------------------
tdkpshell-static: libtdkp.a bin/tdkpshell.o bin/tdkpshell_wrap.o
	$(CXX) $(CXX_FLAGS) $(CXX_INC) -o tdkpshell.bin \
	bin/tdkpshell.o bin/tdkpshell_wrap.o \
	libtdkp.a $(CXX_LIBS)

bin/sp.o: sp.c tdkp/utilities/SchroedingerPoisson.h
	$(CXX) -c $(CXX_FLAGS) -o $@ sp.c

bin/ModifiedLapackRoutines.o: tdkp/solvers/ModifiedLapackRoutines.f
	$(F77) -c $(F77_FLAGS) -o $@ $<

bin/tdkpshell.o: tdkp/tcl/tdkpshell.cpp $(HEADERS)
	$(CXX) -c $(CXX_FLAGS) $(CXX_INC) -o $@ $<

bin/tdkpshell_wrap.cpp: tdkp/tcl/tdkpshell.i $(HEADERS)
	$(SWIG) -Wall -c++ -tcl $(CXX_INC) -o $@ $<

bin/tdkpshell_wrap.o: bin/tdkpshell_wrap.cpp
	$(CXX) -c $(CXX_FLAGS) $(CXX_INC) -o $@ $<


ifneq ($(SVNSRCPATH),)
bin/Version.h:
	rm -f bin/Version.h
	touch bin/Version.h
	echo "#define TDKP_REVISION "`cd $(SVNSRCPATH);svn info | grep Revision | sed -e 's/Revision://'` >> bin/Version.h
	echo "#define MACHINE \"`hostname`\"" >> bin/Version.h
else
bin/Version.h:
	rm -f bin/Version.h
	touch bin/Version.h
	echo "#define TDKP_REVISION \"0.0\"" >> bin/Version.h
	echo "#define MACHINE \"`hostname`\"" >> bin/Version.h
endif

bin/write_mat_to_matlab_readable: tdkp/utilities/write_mat_to_matlab_readable.cc
	$(CXX) $(CXX_FLAGS) -o $@ $^ $(CXX_INC) -ltdkp $(CXX_LIBS)

