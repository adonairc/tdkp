# SYSTEM-SPECIFIC PART OF TDKP MAKEFILE

# blanc.cscs.ch: note that tdkpshell.bin cannot be done --> use make libtdkp.a instead


# -----------------------------------------------------------------------------------
# external include directories and libraries
# -----------------------------------------------------------------------------------
LIB_PATH       = /users/steiger/src/release/gnu/lib

# for omp.h
GCC_INC        = -I/opt/ibmcmp/xlsmp/1.7/include/
GCC_LIB_PATH   = 

BOOST_INC      = -I../boost_1_35_0/ 
BOOST_LIB      = # no lib necessary, only template stuff is used

DFISE_INC      = -I../DF-ISE/c++/
DFISE_LIB      = -L$(LIB_PATH) -lDF-ISE++ -lDF-ISE

LINALG_INC     = 
LINALG_LIB     = /apps/lapack/3.1.1/XL/lib/lapack_ppc64.a /apps/lapack/3.1.1/XL/lib/blas_ppc64.a

FLENS_INC      = -I../FLENS-lite_0608/ -I../cblas/src/
FLENS_LIB      = # no lib necessary, only template stuff is used

UMFPACK_INC    = -I../umfpack/UMFPACK/Include/ -I../umfpack/UFconfig/ -I../umfpack/AMD/Include/
UMFPACK_LIB    = -L$(LIB_PATH) -lumfpack
UMFPACK_DEF    = -DLINSOLV_INCLUDE_UMFPACK

PARDISO_INC    = 
PARDISO_LIB    = 
PARDISO_DEF    = -DLINSOLV_INCLUDE_PARDISO

SUPERLU_INC    = 
SUPERLU_LIB    = 
SUPERLU_DEF    = -DLINSOLV_INCLUDE_SUPERLU

AZTEC00_INC    = 
AZTEC00_LIB    = 
AZTEC00_DEF    = -DLINSOLV_INCLUDE_AZTECOO 

VEPREK_PATH      = 
ATZE_MPI_INC     = 
ATZE_MPI_LIBS 	 = 
ATZE_SERIAL_INC  =
ATZE_SERIAL_LIBS = 

METIS_INC      = 
METIS_LIB      = 

ILS_INC        = 
ILS_LIB        = 
ILS_DEF        = -DLINSOLV_INCLUDE_ILS

TCL_INC        = 
TCL_LIB        = 

ZLIB_INC       = -I/apps/zlib/1.2.3/XL/include
ZLIB_LIB       = -L/apps/zlib/1.2.3/XL/lib -lz

ARPACK_INC     = 
ARPACK_LIB     = 
ARPACK_DEF     = -DEIGNSOLV_INCLUDE_ARPACK

JAC_DAV_DEF    = -DEIGNSOLV_INCLUDE_JACOBI_DAVIDSON
JDQZ_DEF       = -DEIGNSOLV_INCLUDE_JDQZ

# -----------------------------------------------------------------------------------
# exclude unwanted solvers in the next section
# -----------------------------------------------------------------------------------

EIGSOLVERS    = $(JDQZ_DEF)
# not defined: $(ARPACK_DEF) $(JAC_DAV_DEF) 

LINSOLVERS     = $(UMFPACK_DEF) 
# not defined: $(ILS_DEF) $(PARDISO_DEF) $(SUPERLU_DEF) $(AZTEC00_DEF)

CXX_INC        = $(UMFPACK_INC) $(METIS_INC) $(DFISE_INC) $(ZLIB_INC) $(BOOST_INC) $(FLENS_INC) $(LINALG_INC) $(GCC_INC) -I.
# not included: $(ILS_INC) $(PARDISO_INC) $(SUPERLU_INC) $(AZTEC00_INC)

CXX_LOCAL_LIBS = $(UMFPACK_LIB) $(METIS_LIB) $(DFISE_LIB) $(ZLIB_LIB) $(ARPACK_LIB) -ljdqz $(LINALG_LIB)
# not linked: $(ILS_LIB) $(PARDISO_LIB) $(SUPERLU_LIB) $(AZTEC00_LIB)

CXX_LIBS       = -nodefaultlibs $(CXX_LOCAL_LIBS) -L$(GCC_LIB_PATH) -lstdc++ -lm -lnuma -lc 
CXX            = mpCC -compiler gcc
F77            = mpfort -compiler gcc
SWIG           = 
MPICC          = 

DEFINES        = -DNDEBUG -DNOACML -DNOMETIS -DNOUNDERSCORE -DNOREMOTESOLVER $(LINSOLVERS) $(EIGSOLVERS)
CXX_FLAGS      = -q64 -Wall -Wno-unknown-pragmas -fno-strict-aliasing -O3 $(DEFINES)
F77_FLAGS      = -q64 -O3 -fno-underscoring 
