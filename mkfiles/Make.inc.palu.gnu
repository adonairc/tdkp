# SYSTEM-SPECIFIC PART OF TDKP MAKEFILE
# palu.cscs.ch: note that tdkpshell.bin cannot be done --> use make libtdkp.a instead

# -----------------------------------------------------------------------------------
# external include directories and libraries
# -----------------------------------------------------------------------------------
LIB_PATH       = /users/steiger/src/release/gnu/lib

# GCC include is used in make_depend.sh
GCC_INC        =
GCC_LIB_PATH   = 

BOOST_INC      = -I../boost_1_35_0/ 
BOOST_LIB      = # no lib necessary, only template stuff is used

DFISE_INC      = -I../DF-ISE/c++/
DFISE_LIB      = -L$(LIB_PATH) -lDF-ISE++ -lDF-ISE

LINALG_INC     = -I/opt/acml/4.0.1a/gfortran64/include/
#LINALG_LIB     = -lacml -lacml_mv
LINALG_LIB     = -lacml_mv

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

METIS_INC      = -I../metis/
METIS_LIB      = -L$(LIB_PATH) -lmetis

ILS_INC        = 
ILS_LIB        = 
ILS_DEF        = -DLINSOLV_INCLUDE_ILS

TCL_INC        = 
TCL_LIB        = 

ZLIB_INC       = -I../zlib-1.2.3/
ZLIB_LIB       = -L$(LIB_PATH) -lz

ARPACK_INC     = 
ARPACK_LIB     = -larpack 
ARPACK_DEF     = -DEIGNSOLV_INCLUDE_ARPACK

JAC_DAV_DEF    = -DEIGNSOLV_INCLUDE_JACOBI_DAVIDSON
JDQZ_DEF       = -DEIGNSOLV_INCLUDE_JDQZ

# -----------------------------------------------------------------------------------
# exclude unwanted solvers in the next section
# -----------------------------------------------------------------------------------

# -DNOREMOTESOLVER to exclude RemoteGEVPSolver functionality
# -DNOMETIS disables compilation of GraphReordering class
# -DNOUNDERSCORE switches to fortran routine names w/o trailing underscores
# -DNOACML disables ACML functionality

EIGSOLVERS    = $(JDQZ_DEF)
# not defined: $(ARPACK_DEF) $(JAC_DAV_DEF) 

LINSOLVERS     = $(UMFPACK_DEF) 
# not defined: $(ILS_DEF) $(PARDISO_DEF) $(SUPERLU_DEF) $(AZTEC00_DEF)

DEFINES        = -DNDEBUG -DNOREMOTESOLVER -DNOMETIS $(LINSOLVERS) $(EIGSOLVERS)

CXX_INC        = $(UMFPACK_INC) $(METIS_INC) $(DFISE_INC) $(ZLIB_INC) $(BOOST_INC) $(FLENS_INC) $(LINALG_INC) -I.
# not included: $(ILS_INC) $(PARDISO_INC) $(SUPERLU_INC) $(AZTEC00_INC)

CXX_LOCAL_LIBS = $(UMFPACK_LIB) $(METIS_LIB) $(DFISE_LIB) $(ZLIB_LIB) $(ARPACK_LIB) -ljdqz $(LINALG_LIB)
# not linked: $(ILS_LIB) $(PARDISO_LIB) $(SUPERLU_LIB) $(AZTEC00_LIB)

CXX_LIBS       = -nodefaultlibs $(CXX_LOCAL_LIBS) -L$(GCC_LIB_PATH) -lgfortran -lstdc++ -lm -lnuma -lc 
CXX            = CC
F77            = ftn
SWIG           = 
MPICC          = 

CXX_FLAGS      = -Wall -Wno-unknown-pragmas -fno-strict-aliasing -O3 $(DEFINES)
F77_FLAGS      = -O3 
