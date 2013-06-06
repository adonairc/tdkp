
include ../Make.inc

#CXX  = g++-4.0.2
TESTGEN = ./cxxtest/cxxtestgen.pl

#CXX_WARNINGS = -Wall
#CXX_DEBUG    = -g
#CXX_LIBS     =  -L../../DF-ISE/c++/ -L../../DF-ISE/src/ -lz -lgfortran
#CXX_INC_LOCAL =  -Icxxtest/ -I. -I../../boost -I../../DF-ISE/c++ -I../../DF-ISE/include
CXX_INC_LOCAL =  -Icxxtest/ -I. -I../../sebise/sebise $(CXX_INC)


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
	TDKPSHELL = libtdkp.so tdkpshell atze
	CXX_LIBS += $(CXX_LIBS_SHARED)
else
	TDKPSHELL = tdkpshell-static atze-static
	CXX_LIBS += $(CXX_LIBS_STATIC)
endif




all: clean_output make_tests

clean_output:
	rm -f *.log

clean:
	rm -f bin/* *.log */runner.cpp *.dat *.txt *.bin *.out *.mat

make_tests: bin/02_Elements bin/03_Geometry bin/04_PropertyContainer bin/05_MaterialDatabase \
		   bin/06_SYMCSRMatrix bin/07_EigenSolvers bin/08_KPMatrices bin/09_Dot3D bin/10_DataIOTest \
		   bin/11_Bandstructure bin/12_BulkBandstructure bin/13_LinearSolvers \
		   bin/15_DensityOfStates bin/16_DomainIntegrator bin/17_MatrixElements bin/18_CoulombMatrix \
		   bin/19_IsoSurface bin/20_Splines bin/21_ElementBoundary bin/23_SLCMBShift

tests: make_tests
	./bin/03_Geometry
	./bin/04_PropertyContainer
	./bin/05_MaterialDatabase
	./bin/06_SYMCSRMatrix
	./bin/07_EigenSolvers
	./bin/08_KPMatrices
	./bin/09_Dot3D
	./bin/10_DataIOTest
	./bin/11_Bandstructure
	./bin/12_BulkBandstructure
	./bin/13_LinearSolvers
	./bin/15_DensityOfStates
	./bin/16_DomainIntegrator
	./bin/17_MatrixElements
	./bin/18_CoulombMatrix
	./bin/19_IsoSurface
	./bin/20_Splines
	./bin/21_ElementBoundary
	./bin/23_SLCMBShift      
	./bin/02_Elements




########################
# unit 01
#######################
bin/01_CSCtest: 01_CSCtest/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I01_CSCtest/ -I.. -o $@ $^

01_CSCtest/runner.cpp: 01_CSCtest/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 02
#######################
bin/02_Elements: 02_Elements/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I02_Elements/ -I.. -o $@ $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

02_Elements/runner.cpp: 02_Elements/definition.h 02_Elements/validator.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 03
#######################
bin/03_Geometry: 03_Geometry/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I02_Elements/ -I.. -o $@ \
	       $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

03_Geometry/runner.cpp: 03_Geometry/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^


########################
# unit 04
#######################
bin/04_PropertyContainer: 04_PropertyContainer/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I04_PropertyContainer/ -I.. -I../../boost -o $@ $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

04_PropertyContainer/runner.cpp: 04_PropertyContainer/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^



########################
# unit 05
#######################
bin/05_MaterialDatabase: 05_MaterialDatabase/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I05_MaterialDatabase/ -I.. -I../../boost -o $@ $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

05_MaterialDatabase/runner.cpp: 05_MaterialDatabase/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 06
#######################
bin/06_SYMCSRMatrix: 06_SYMCSRMatrix/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I06_SYMCSRMatrix/ -I.. -o $@ $^ \
	$(CXX_LOCAL_LIBS) $(CXX_LIBS)

06_SYMCSRMatrix/runner.cpp: 06_SYMCSRMatrix/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 07
#######################
bin/07_EigenSolvers: 07_EigenSolvers/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) $(CXX_DEFINES) -I07_EigenSolvers/ -I.. -o $@ $^ \
		$(CXX_LOCAL_LIBS) $(CXX_LIBS)

07_EigenSolvers/runner.cpp: 07_EigenSolvers/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 08 - kp matrices
#######################
bin/08_KPMatrices: 08_KPMatrices/runner.cpp ../libtdkp.so

	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I08_KPMatrices/ -I.. \
	           -I../src/geometry -Icommon/ -I../src/io -o $@ $^ \
				$(CXX_LOCAL_LIBS) $(CXX_LIBS)


08_KPMatrices/runner.cpp: 08_KPMatrices/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^


########################
# unit 09 - 3d dot
#######################
bin/09_Dot3D: 09_Dot3D/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I08_KPMatrices/ -I.. -I../../DF-ISE/c++/ -I../../boost/ -o $@ $^ \
	    $(CXX_LOCAL_LIBS) $(CXX_LIBS)

09_Dot3D/runner.cpp: 09_Dot3D/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 10 - data io
#######################
bin/10_DataIOTest: 10_DataIOTest/runner.cpp ../libtdkp.so
	$(CXX)  $(CXX_FLAGS) $(CXX_INC_LOCAL) -I10_DataIOTest/ -I.. -o $@ $^ \
	    $(CXX_LOCAL_LIBS) $(CXX_LIBS) -lz

10_DataIOTest/runner.cpp: 10_DataIOTest/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

#######################
# unit 11 - bandstructure
#######################
bin/11_Bandstructure: 11_Bandstructure/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I11_Bandstructure/ -I.. -o $@ $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

11_Bandstructure/runner.cpp: 11_Bandstructure/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^


########################
# unit 12 - bulk bandstructure
########################
bin/12_BulkBandstructure: 12_BulkBandstructure/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I12_BulkBandstructure/ -I.. -Icommon/ -o $@ $^ $(CXX_LOCAL_LIBS) $(CXX_LIBS)

12_BulkBandstructure/runner.cpp: 12_BulkBandstructure/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

###########################
# unit 13 Linear solveres
###########################
bin/13_LinearSolvers: 13_LinearSolvers/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) $(CXX_DEFINES) -I../../SuperLU_MT/SRC -I../../ils-2.0/src -Icommon -I13_LinearSolvers/ -I.. -o $@ $^ \
		$(CXX_LOCAL_LIBS) $(CXX_LIBS)

13_LinearSolvers/runner.cpp: 13_LinearSolvers/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

###########################
# unit 14 aqua interface
###########################
bin/14_AquaInterface: 14_AquaInterface/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL)  -I../../ils-2.0/src -Icommon -I.. -o $@ $^ \
		$(CXX_LOCAL_LIBS) $(CXX_LIB_FLENS) $(CXX_LIBS)

14_AquaInterface/runner.cpp: 14_AquaInterface/definition.h common/MatrixCreation.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

############################
# unit 15 density of states
############################

bin/15_DensityOfStates: 15_DensityOfStates/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I15_DensityOfStates/ -I.. -Icommon -o $@ $^ \
	$(CXX_LOCAL_LIBS) $(CXX_LIBS)


15_DensityOfStates/runner.cpp: 15_DensityOfStates/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 16 domain integrator
#######################
bin/16_DomainIntegrator: 16_DomainIntegrator/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I16_DomainIntegrator/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

16_DomainIntegrator/runner.cpp: 16_DomainIntegrator/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 17 matrix elements
#######################
bin/17_MatrixElements: 17_MatrixElements/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I17_MatrixElements/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

17_MatrixElements/runner.cpp: 17_MatrixElements/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 18 CoulombMatrix
#######################
bin/18_CoulombMatrix: 18_CoulombMatrix/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I18_CoulombMatrix/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

18_CoulombMatrix/runner.cpp: 18_CoulombMatrix/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 19 IsoSurface
#######################
bin/19_IsoSurface: 19_IsoSurface/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I19_IsoSurface/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

19_IsoSurface/runner.cpp: 19_IsoSurface/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 20 Splines
#######################
bin/20_Splines: 20_Splines/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I20_Splines/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

20_Splines/runner.cpp: 20_Splines/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 21 Element Boundary
#######################
bin/21_ElementBoundary: 21_ElementBoundary/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I21_ElementBoundary/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

21_ElementBoundary/runner.cpp: 21_ElementBoundary/definition.h 21_ElementBoundary/moc_grid_reader.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^


########################
# unit 22 aqua poisson solver
#######################
bin/22_PoissonAQUA: 22_PoissonAQUA/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I22_PoissonAQUA/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

22_PoissonAQUA/runner.cpp: 22_PoissonAQUA/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

########################
# unit 23 slc mb shifter
#######################
bin/23_SLCMBShift: 23_SLCMBShift/runner.cpp ../libtdkp.so
	$(CXX) $(CXX_FLAGS) $(CXX_INC_LOCAL) -I23_SLCMBShift/ -I.. -o $@ $^  $(CXX_LOCAL_LIBS) $(CXX_LIBS)

23_SLCMBShift/runner.cpp: 23_SLCMBShift/definition.h ../libtdkp.so
	$(TESTGEN) -o $@ --error-printer $^

