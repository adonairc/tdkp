// ------------------------------------------------------------
//
// This file is part of tdkp, a simulation tool for nanostrutctures
// of optoelectronics developed at ETH Zurich
//
// (C) 2005-2009 Ratko G. Veprek, ETH Zurich, veprek@iis.ee.ethz.ch
//
// 1) As of 18.6.2009 this code is property of ETH Zurich and must not be
// transferred, modified or used by third parties without appropriate
// licenses issued by authorized agents of ETH Zurich.
//
// 2) Violation of this will result in judicial action according to civil
// and penal law.
//
// 3) Any claim of authorship other than by the author himself is
// strictly forbidden.
//
// 4) The source code must retain the copyright notice, this list
// of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
// BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
// IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------------------------------------

#ifndef COULOMBMATRIXELEMENT_H_
#define COULOMBMATRIXELEMENT_H_

#include "tdkp/common/all.h"
#include "tdkp/coulomb/CoulombIntegrator.h"
#include "tdkp/coulomb/CoulombFunction.h"
#include "tdkp/coulomb/CoulombBlochStrategy.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/geometry/Geometry.h"

namespace tdkp {



// --------------------------------------------------
// the next classes should be inlined classes in the
// coulomb matrix element object. unfortunately i had
// to move them outside as swig can not handle it.
// --------------------------------------------------

/** coulomb matrix element index holder */
class CMEIdx {
public:	
	CMEIdx(unsigned int n1_, unsigned int n2_, unsigned int n3_, unsigned int n4_)
	: n1(n1_),n2(n2_),n3(n3_),n4(n4_) {}
	unsigned int n1, n2, n3, n4;
	bool operator==(const CMEIdx& rhs) const;		
};	

/** data container for coulomb matrix elements */
class CoulombMatrixElementData : public XYData<double> {
	friend class CoulombMatrixElement;
	friend class CoulombMatrixElementQuantized;
	friend class CoulombMatrixElementBulk;
	
public:
	CoulombMatrixElementData(const char* filename);
	void write_binary(ostream& stream) const;
	void write_binary(const char* filename) const;
	void read_binary(istream& stream);
	const vector<double>& get_q_values() const { return q_values; }
	const vector<cplx>& get_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4) const;
	
	// ---------------------------------------------
	// wavefunction info functions
	// ---------------------------------------------
	const string& get_wavefunction_type(unsigned int ni) const;
	unsigned int  get_num_wavefunctions() const;
		
	// ---------------------------------------------
	// XYData functions
	// ---------------------------------------------	
	virtual int    get_x_length()                const;
	virtual int    get_num_y_sets()              const;
	virtual int    get_num_x_sets()              const; 
	virtual void   get_x(int xidx, vector<double> &x) const;
	virtual void   get_y(int yidx, vector<double>& y) const; 
	virtual string get_x_identifier(int xidx)    const;
	virtual string get_y_identifier(int yidx) 	 const;				
protected:
	CoulombMatrixElementData() {} // protected, may empty only be created via coulomb matrix element
	vector<string>		  wavefunction_types;	
	vector<double>        q_values;
	vector<CMEIdx>        requests;
	vector<vector<cplx> > coulomb_matrix_elements;
private:		
	static const int magic;		
};	

/** coulomb matrix element evaluator 
 *
 * the coulomb matrix element is defined by
 * <n1 n2 | v | n3 n4> = int dr dr' g*(n1;r) g*(n2;r') v(r,r') g(n4;r) g(n3;r') 
 * 
 * be careful, we return for bulk, wells and wires the fourier transformed
 * coulomb potential
 * 
 * bulk:
 *   that would be 4pi/q^2 <g1g2g3g4>. to prevent numerical difficulties, we omit the 1/q^2 
 *   (SO ITS YOUR JOB TO ADD THAT APPROPRIATELY)
 * 
 * well:
 *   as in the bulk case, we return 2pi*exp(-q|z-z') and OMIT 1/q
 * 
 * wire:
 *   here, we really return the full coulomb matrix element 
 *  
 */
class CoulombMatrixElement {
public:
	CoulombMatrixElement(unsigned int wavefunction_length);
	virtual ~CoulombMatrixElement();
	// ---------------------------------------------
	// setup functions
	// ---------------------------------------------
	void clear_all();
	void clear_wavefunctions();
	void clear_requests_and_results();
	void set_q_range(double q_min, double q_max, unsigned int num_q_values, double progression = 1.0);
	void set_q_range(vector<double>& q_values);	
	unsigned int add_wavefunction(const EigenSolution<complex<double> >* wave, const char* name = "ud");
	// ---------------------------------------------
	// calculation request functions
	// ---------------------------------------------	
	void request_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4);
	virtual void calculate() = 0; 
	// ---------------------------------------------
	// result access
	// ---------------------------------------------
	const vector<double>& get_q_values() const { return container.get_q_values(); }
	const vector<cplx>& get_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4);
		
	// -----------------------------------------------
	// read and saveable data class
	// -----------------------------------------------
	const CoulombMatrixElementData& get_data_object() const { return container; } 
		
protected:
	const EigenSolution<complex<double> >* get_wavefunction(unsigned int idx) const;
	
	/** the integration over z index holder (so for n1/n4) */
	class TwoInt {
	public:
		TwoInt(unsigned int n1_, unsigned int n4_) : n1(n1_), n4(n4_) {}
		unsigned int n1, n4;
		bool operator==(const TwoInt& rhs) const;
	};
	
	vector<const EigenSolution<cplx>*> wavefunctions;

	CoulombMatrixElementData	container;
private:
	unsigned int wavefunction_length;	
			
};

/** coulomb matrix element evaluator for quantized systems 
 *
 */
class CoulombMatrixElementQuantized : public CoulombMatrixElement {
		
public:
	CoulombMatrixElementQuantized(
		const Geometry& geometry, 
		MaterialDatabase& material_database,
		unsigned int kp_basis_size				
	);
	virtual ~CoulombMatrixElementQuantized();
	
	void calculate();
	 
private:
	
	const Geometry&                    geometry;
	MaterialDatabase&                  material_database;
	CoulombBlochStrategy*              bloch_strategy;

};

/** coulomb matrix element evaluator for bulk */
class CoulombMatrixElementBulk : public CoulombMatrixElement {		
public:
	CoulombMatrixElementBulk();
	virtual ~CoulombMatrixElementBulk();	
	void calculate();


};

}

#endif /*COULOMBMATRIXELEMENT_H_*/
