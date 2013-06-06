#ifndef DEFINITION_H_
#define DEFINITION_H_

#define TDKP_UNITTEST

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"


// ------------------------------------------
// class includes
// ------------------------------------------
#include <math.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>
#include "tdkp/utilities/MatrixElements.h"
#include "tdkp/main/EigenSolution.h"

using namespace tdkp;

namespace tdkp {

/** faky momentum operator with diagonal tensor */	
class MomentumOperatorBulkMoc : public MomentumOperator {
public:		
	MomentumOperatorBulkMoc(bool k_dependent_, unsigned int basis_size,  unsigned int solution_length);
	virtual ~MomentumOperatorBulkMoc();		
	virtual void lock();
	virtual void release();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual bool is_k_dependent() const { return k_dependent; }
	virtual void set_k_value(const DomainPoint& point);
	virtual void eval(vector<cplx>& result, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 2; }
private:
	bool k_dependent;	
	unsigned int basis_size;
	unsigned int solution_length;
	
};
	
}


class MatrixElementsTest : public CxxTest::TestSuite {

public:			
	MatrixElementsTest() { 
		fout.open("17_MatrixElementsTest_output.log", ios::app | ios::out);		
		Logger::get_instance()->add_listener(&fout);
		//Logger::get_instance()->del_listener(&std::cout);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);	
	}
	
	~MatrixElementsTest() {
		Logger::get_instance()->del_listener(&fout);
		fout.close();	
	}
			
	std::fstream fout;
		
	void setUp() {}				
	void tearDown() {}
	
	 
	
	void test_moc_matrix_elements();
	
};

MomentumOperatorBulkMoc::MomentumOperatorBulkMoc(bool k_dependent_, unsigned int basis_size_, unsigned int solution_length_) 
: k_dependent(k_dependent_),
  basis_size(basis_size_),
  solution_length(solution_length_)  
{
	
}

MomentumOperatorBulkMoc::~MomentumOperatorBulkMoc() {
	
}		
void MomentumOperatorBulkMoc::lock() {

}
void MomentumOperatorBulkMoc::release() {
	
}
bool MomentumOperatorBulkMoc::compatible(short cb_basis_size_, int cb_length, short vb_basis_size_, int vb_length) const {
	
	ostringstream sout;
	if((signed)basis_size != cb_basis_size_ || (signed)basis_size != vb_basis_size_) {		
		sout << "basis_size: " << basis_size << " cb_basis_size: " << cb_basis_size_ 
		     << "  vb_basis_size: " << vb_basis_size_ << endl;
		Logger::get_instance()->emit(LOG_WARN, sout.str());		     
		return false;	
	} 
	if(cb_length != (signed)solution_length || vb_length != (signed)solution_length) {
		sout << "solution_length:  " << solution_length << " cb length: " << cb_length << " "
		     << "vb_length: " << vb_length << endl;
		Logger::get_instance()->emit(LOG_WARN, sout.str());				
		return false;	
	}
	return true;	
}

void MomentumOperatorBulkMoc::set_k_value(const DomainPoint& point) {
	
}
void MomentumOperatorBulkMoc::eval(vector<cplx>& result, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {
	
	TDKP_ASSERT(cb_wave.size() == vb_wave.size(), "cb_wave.size() == vb_wave.size()");
	TDKP_ASSERT(cb_wave.size() == solution_length * basis_size, "cb_wave.size() (" << cb_wave.size() << ") == solution_length (" << solution_length << ") * basis_size (" << basis_size << ")");  
	
	cplx ret(0.0);
	for(unsigned int ii = 0; ii < solution_length * basis_size; ii++) {
		ret += conj(cb_wave[ii]) * vb_wave[ii];		
	}	
	result[0] = abs(ret * (static_cast<double>(get_operator_direction()) + 1.0));
	result[1] = ret * (static_cast<double>(get_operator_direction()) + 1.0);
}


void MatrixElementsTest::test_moc_matrix_elements() {
	
	DomainMaster domain;
	create_3D_domain_radial(domain, Vector3D(1.0, 0.0, 0.0), 0.0, 1.2, 8);
	
	// basis size == 4, num subbands == 2, solution length == 2 (* basis_size)
	BandstructureDomain<cplx> cb_bands(4, 2, 2, domain);
	// basis size == 4, num subbands == 4, solution length == 2	(* basis_size)
	BandstructureDomain<cplx> vb_bands(4, 4, 2, domain);
		
	// ----------------------------------------------------
	// build bandstructures
	// matrix elements between cb <-> vb
	// <cb0|M|vb0> == 1
	// <cb0|M|vb1> == <cb1|M|vb0> == 0
	// <cb1|M|vb1> == 1
	// <cb0|M|vb2/3> == <cb1|M|vb2/3> == 0
	// respectively, the momentum operator retrusn <cb0Mvb0> * (dir_p + 1)
	// ----------------------------------------------------	
	vector<cplx> tmp(8);
	tmp[0] = 1.0;
	tmp[4] = 1.0;
	for(unsigned int kk = 0; kk < 8; kk++) {
		cb_bands.add_eigensolution(kk, 0, new EigenSolution<cplx>(0.1, &tmp[0], 2, 4));
		vb_bands.add_eigensolution(kk, 0, new EigenSolution<cplx>(0.1, &tmp[0], 2, 4));		
	}
	tmp.assign(8, 0.0);
	tmp[1] = 1.0;
	tmp[5] = 1.0;
	for(unsigned int kk = 0; kk < 8; kk++) {
		cb_bands.add_eigensolution(kk, 1, new EigenSolution<cplx>(0.1, &tmp[0], 2, 4));
		vb_bands.add_eigensolution(kk, 1, new EigenSolution<cplx>(0.1, &tmp[0], 2, 4));		
	}
	tmp.assign(8, 0.0);
	for(unsigned int kk = 0; kk < 8; kk++) {
		vb_bands.add_eigensolution(kk, 2, new EigenSolution<cplx>(0.2, &tmp[0], 2, 4));
		vb_bands.add_eigensolution(kk, 3, new EigenSolution<cplx>(0.2, &tmp[0], 2, 4));		
	}	
	
	try {
		MomentumOperatorBulkMoc m_op(true, 4, 2); 
		MatrixElements melem(m_op);	
		melem.calculate(cb_bands, vb_bands);
		
		TS_ASSERT(melem.get_num_cb_bands() == 2);
		TS_ASSERT(melem.get_num_vb_bands() == 4);
		TS_ASSERT(melem.get_num_k_values() == 8);
		
		for(unsigned int dd = 0; dd < 3; dd++) {
			for(unsigned int kk = 0; kk < 8; kk++) {				
				TS_ASSERT(melem.get_abs_square(0,0,dd,kk) == 2.0 * (1.0 + dd));
				TS_ASSERT(melem.get_abs_square(1,1,dd,kk) == 2.0 * (1.0 + dd));
				TS_ASSERT(melem.get_abs_square(0,1,dd,kk) == 0.0);
				TS_ASSERT(melem.get_abs_square(1,0,dd,kk) == 0.0);
				TS_ASSERT(melem.get_abs_square(0,2,dd,kk) == 0.0);
				TS_ASSERT(melem.get_abs_square(0,3,dd,kk) == 0.0);
				TS_ASSERT(melem.get_abs_square(1,2,dd,kk) == 0.0);
				TS_ASSERT(melem.get_abs_square(1,3,dd,kk) == 0.0);				
			}
		} 

			
	} catch(Exception *e) {
		TS_FAIL(e->get_reason());
	}
	

	
}


#endif /*DEFINITION_H_*/
