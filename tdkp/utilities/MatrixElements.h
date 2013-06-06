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

#ifndef MATRIXELEMENTS_H_
#define MATRIXELEMENTS_H_

#include <vector>
#include <string>
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/Domain.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/main/FEMSolverLE.h"
#include "tdkp/probdefs/LinearProblem.h"
#include "tdkp/kpmatrices/KPMatrixBase.h"
#include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix6x6Wurtzite.h"
#include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"

namespace tdkp {
	
/** momentum operator definition for matrix element calculation 
 * 
 * the aim of this class is to provide the matrix element calculation
 * a strategy for the passed bandstructure.
 * 
 * the tensor takes two wavefunctions |c> and |v> and returns
 * |<c|p|v>|^2 and <c|p|v>
 * 
 * if |c> is effmass, it returns 3 values,
 * |<c spin up|p|v>|^2 + |<c spin down|p|v>|^2, <c spin up|p|v>, <c spin down|p|v>  
 * 
 */
class MomentumOperator {
public:
	vector<cplx>::const_iterator data_const_iterator;

	MomentumOperator();
	virtual ~MomentumOperator();
	/** prepares class and geometry object. lock is set bevore eval is used and released after last evaluation. */
	virtual void lock() = 0;
	/** releases class data and resets geometry object */
	virtual void release() = 0;
	/** test whether the momentum operator is compatible to our bands */ 
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const = 0;
	/** returns true if the momentum operator depends on k */
	virtual bool is_k_dependent() const = 0;
	/** set current point in k space (throws exception if operator is not k dependent and if class is locked) */
	virtual void set_k_value(const DomainPoint& point) = 0;
	/** set momentum operator direction */
	void set_operator_direction(unsigned short dir) { TDKP_ASSERT(dir < 3, "dir < 3"); momentum_operator_direction = dir; }
	/** return operator direction */	
	unsigned short get_operator_direction() const { return momentum_operator_direction; }  
	/** evaluate tensor */
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const = 0;
	/** return length of results we will return in eval */
	virtual int eval_results_length() const = 0;
	/** return a string giving a hint what kind of operator this is */
	virtual string get_description() const { return string("the unknown momentum operator"); }
	
				
private:
	unsigned short    momentum_operator_direction;	 
};

/** base class to calculate momentum matrix elements
 *  
 *  we are primarily interested to calculate 
 *  |<c|p|v>|^2 from the kp solutions for slc, and
 *  <c|p|v> for clc
 */  
class MatrixElements : public XYData<double> {

public:
	MatrixElements(MomentumOperator& momentum_operator_);
		
	// --------------------------------------------------
	// the value iterators
	// --------------------------------------------------
	typedef std::vector<double>::iterator 		value_iterator;
	typedef std::vector<double>::const_iterator const_value_iterator;

	// --------------------------------------------------
	// constructors
	// --------------------------------------------------	
	virtual ~MatrixElements();
	
	// --------------------------------------------------
	// solving
	// --------------------------------------------------
	void calculate(unsigned int max_cb_bands, unsigned int max_vb_bands, const BandstructureDomain<complex<double> >& cb_bands, const BandstructureDomain<complex<double> >& vb_bands);   
	void calculate(const BandstructureDomain<complex<double> >& cb_bands, const BandstructureDomain<complex<double> >& vb_bands);
	
	// --------------------------------------------------
	// special solvers ...
	// --------------------------------------------------
	void calculate_bulk_effmass(const BandstructureDomain<complex<double> >& vb_bands);
	
	// --------------------------------------------------
	// create dispersive effective mass matrix elements
	// --------------------------------------------------
	MatrixElements get_disp_matrix_elements(const DomainMaster& new_domain);
		  	
	// --------------------------------------------------
	// solution properties
	// --------------------------------------------------		
	unsigned int get_num_cb_bands() const;
	unsigned int get_num_vb_bands() const; 
	unsigned int get_num_k_values() const;
	
	// --------------------------------------------------
	// raw value access
	// --------------------------------------------------
	const DomainMaster&    get_domain() const { return domain; }
	/** return the square absolute value of the matrix element |<c|p|v>|^2 */
	const double& get_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk) const;
	/** return the square absolute value of the matrix element |<c|p|v>|^2 */
	double&       get_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk);
	
	/** return the matrix element <c|p|v> */
	const complex<double>& get(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx) const;
	/** return the matrix element <c|p|v> */
	complex<double>&       get(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx);
			
	/** set the matrix element |<c|p|v>|^2 */		
	void set_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, const double& value);
	/** set the matrix element <c|p|v> (idx value_idx)*/
	void set(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx, const complex<double>& value);
	

	void read_binary(const char* filename);
	void read_binary(istream& in);
	void write_binary(const char* filename) const;
	void write_binary(ostream& out) const;
		
	// -----------------------------------------------
	// interface XYData (returns square matrix elements)
	// -----------------------------------------------
	virtual int         get_x_length()                const;
	virtual int         get_num_y_sets()              const;
	virtual int         get_num_x_sets()              const {  return this->domain.get_dimension(); }
	virtual void        get_x(int xidx, std::vector<double> &x) const;
	virtual void        get_y(int yidx, std::vector<double>& y) const; 
	virtual std::string get_x_identifier(int xidx)   const;
	virtual std::string get_y_identifier(int yidx) 	 const;	

private:	
	unsigned int get_values_idx(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk) const;
		
	MomentumOperator& momentum_operator;
	DomainMaster      domain;

	unsigned int num_cb_bands; 
	unsigned int num_vb_bands;
						
	std::vector<double> values_square;
	std::vector<complex<double> > values_raw;
		
};


/** base class for bulk matrix elements 
 *
 * for the bulk matrix elements, we skip any strain dependence except for
 * the 8x8 models, where strain dependence can be included via setting 
 * a specific strain to the kp matrix.
 */
class MomentumOperatorBulk : public MomentumOperator {
public:
	virtual ~MomentumOperatorBulk() {};	
	virtual void set_k_value(const DomainPoint& point);	
	virtual void release();
	virtual bool is_k_dependent() const { return false; }	
protected:
	MomentumOperatorBulk();	
	bool locked;	
	vector<cplx> tensor;
	double kvalues[3];	
};

/* momentum operator for the bulk 8x8 matrix element (wurtzite or zincblende) */
class MomentumOperatorBulk8x8 : public MomentumOperatorBulk {
public:
	MomentumOperatorBulk8x8(const KPMatrix8x8EndersForeman& matrix);
	MomentumOperatorBulk8x8(const KPMatrix8x8Wurtzite& matrix);
	virtual ~MomentumOperatorBulk8x8() {};
	virtual void lock();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	/** set whether we use the simple gamma point approximation */
	void set_k_dependent(bool k_dependent_) { TDKP_ASSERT(!locked, "object is locked! so you can not set anything to it!"); k_dependent = k_dependent_; }				
	virtual bool is_k_dependent() const { return k_dependent; }
	virtual int eval_results_length() const { return 2; }	// no effmass
	
private:	
	const KPMatrixBase& matrix;
	bool k_dependent;

};


/** momentum operator for bulk 6x6 matrix element */
class MomentumOperatorBulk6x6 : public MomentumOperatorBulk {
public:
	MomentumOperatorBulk6x6(const double& optical_matrix_parameter);
	MomentumOperatorBulk6x6(const Material& material);
	virtual ~MomentumOperatorBulk6x6() {};
	virtual void lock();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 3; }	
private:
	double optical_matrix_param;	
};


/** momentum operator for bulk 6x6 wurtzite matrix element */
class MomentumOperatorBulk6x6WZ : public MomentumOperatorBulk {
public:
	MomentumOperatorBulk6x6WZ(const Material& material);
	virtual ~MomentumOperatorBulk6x6WZ() {};
	virtual void lock();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 3; }	
private:
	/** optical momentum matrix params */
	double P1, P2;
};

/** momentum operator for bulk 4x4 matrix element */
class MomentumOperatorBulk4x4 : public MomentumOperatorBulk {
public:
	MomentumOperatorBulk4x4(const double& optical_matrix_parameter);
	MomentumOperatorBulk4x4(const Material& material);
	virtual ~MomentumOperatorBulk4x4() {};
	virtual void lock();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 3; }
private:
	double optical_matrix_param;
};

/** bulk momentum operator used for effective mass calculation */
class MomentumOperatorBulkEffectiveMass : public MomentumOperatorBulk {
public:
	MomentumOperatorBulkEffectiveMass(const double& optical_matrix_parameter);
	virtual ~MomentumOperatorBulkEffectiveMass() {};
	virtual void lock();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 2; }
private:
	double optical_matrix_param;
};


/** base momentum operator for pure effective mass kp 4x4 and 6x6 calculations */
class MomentumOperatorQuantized : public MomentumOperator, public LinearProblem<complex<double> > {
public:
	MomentumOperatorQuantized(int num_equation_per_node, short cb_basis_size, short vb_basis_size, Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperatorQuantized();

	// -------------------------------------------------
	// inherited from momentum operator
	// -------------------------------------------------	
	virtual void lock();
	virtual void release();
	virtual bool compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const;
	virtual bool is_k_dependent() const { return false; }
	virtual void set_k_value(const DomainPoint& point);
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const = 0;
	virtual int  eval_results_length() const { return 3; }	

	// -------------------------------------------------------
	// inherited from linear problem, used for matrix assembly
	// -------------------------------------------------------
	virtual void prepare() = 0;
	virtual void display_solution_info() const {};	
		
	// -------------------------------------------------------
	// factory!
	// -------------------------------------------------------
	template<class CBPC, class VBPC>
	static MomentumOperatorQuantized* factory(Geometry& geometry, MaterialDatabase& material_database);		

	/** interoperation with FEMSolverLE (to assemble momentum matrix) */
	virtual void calculate_element_matrices(const Element* elem, cplx* lhs, cplx* rhs, int* node_internal_indices, int &n) const = 0;

protected:
	FEMSolverLE<cplx> fem_assembler;
	double k_values[3];
	bool   locked;
	short  obj_cb_basis_size;
	short  obj_vb_basis_size;
		
}; 

/** simple kp XxX momentum operator for quantized systems 
 * 
 * here, we assume p to be constant over k 
 */
class MomentumOperatorSimpleXxX : public MomentumOperatorQuantized {
public:	
	MomentumOperatorSimpleXxX(int num_equation_per_node, short cb_basis_size, short vb_basis_size, Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperatorSimpleXxX() {};		
	virtual void calculate_element_matrices(const Element* elem, cplx* lhs,cplx *rhs, int* node_internal_indices, int &n) const;			
protected:	
	vector<vector<cplx> > momentum_tensors;
				
};

/** kp 6x6 momentum operator for quantized systems
 * 
 * [rv] checked against gebas results 2008-04-30 
 *  
 */
class MomentumOperator6x6 : public MomentumOperatorSimpleXxX {
public:
	MomentumOperator6x6(Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperator6x6();
	virtual void prepare();
	virtual const int* get_node_sparsity_pattern(int&) const;	
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual string get_unique_identifier() const { return string("MatrixElement6x6"); }	
private:
	static const int sparsity_pattern[12];	
}; 

/** kp 6x6 wurtzite momentum operator for quantized systems
 *
 * [rv] checked against gebas results 2008-04-30 
 */
class MomentumOperator6x6WZ : public MomentumOperator6x6 {
public:
	MomentumOperator6x6WZ(Geometry& geometry, MaterialDatabase& material_database);
	virtual void prepare();
	virtual string get_unique_identifier() const { return string("MatrixElement6x6WZ"); }
};


/** simple momentum operator for kp 4x4 equations 
 *
 * [rv] checked against gebas results 2008-04-30 
 */
class MomentumOperator4x4 : public MomentumOperatorSimpleXxX {
public:
	MomentumOperator4x4(Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperator4x4();
	virtual void prepare();
	virtual const int* get_node_sparsity_pattern(int&) const;	
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual string get_unique_identifier() const { return string("MatrixElement4x4"); }	
private:
	static const int sparsity_pattern[16];	
};

/** basic routines for simple momentum operator for kp 8x8 equations */
class MomentumOperator8x8Base : public MomentumOperatorSimpleXxX {
public:
	MomentumOperator8x8Base(Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperator8x8Base();
	virtual void prepare();
	virtual const int* get_node_sparsity_pattern(int&) const;	
	virtual void eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const;
	virtual int  eval_results_length() const { return 2; }

protected:
	virtual KPMatrixBase* get_kp_matrix() const = 0;
	virtual KPMatrixBase* get_kp_matrix(const Material* material) const = 0;
	void    create_matrix_structures();				
private:
	vector<int> sparsity_pattern;
	vector<KPMatrixBase*> kp_matrices;	
};

/** zincblende 8x8 momentum operator */
class MomentumOperator8x8 : public MomentumOperator8x8Base {
public:
	MomentumOperator8x8(Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperator8x8() {}
	virtual KPMatrixBase* get_kp_matrix() const;
	virtual KPMatrixBase* get_kp_matrix(const Material* material) const;
	virtual string get_unique_identifier() const { return string("MatrixElement8x8"); }	
};

/** wurtzite 8x8 momentum operator
 * 
 * [rv] checked against gebas results 2008-04-30 
 */
class MomentumOperator8x8WZ : public MomentumOperator8x8Base {
public:
	MomentumOperator8x8WZ(Geometry& geometry, MaterialDatabase& material_database);
	virtual ~MomentumOperator8x8WZ() {}
	virtual KPMatrixBase* get_kp_matrix() const;
	virtual KPMatrixBase* get_kp_matrix(const Material* material) const;
	virtual string get_unique_identifier() const { return string("MatrixElement8x8WZ"); }	
};

}

#endif /*MATRIXELEMENTS_H_*/
