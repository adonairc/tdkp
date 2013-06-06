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

#ifndef OVERLAP_H_
#define OVERLAP_H_

#include "tdkp/common/all.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/probdefs/LinearProblem.h"
#include "tdkp/main/FEMSolverME.h"


namespace tdkp {

/** calculates <v1 | v2>
 * 
 * basically <i|j> = \int sum_{a} f^{*}_i,a(z) f_{j},a(z) dz
 */
class Overlap : public LinearProblem<double>{
public:
	Overlap(Geometry& geometry, MaterialDatabase& matdb, SparseMatrixProperties matrix_type = symmetric_matrix);
	virtual ~Overlap();
	
	complex<double> evaluate(
		const EigenSolution<complex<double> >& v1, 
		const EigenSolution<complex<double> >& v2
	) const;
	 	
	/** interoperation with FEMSolverLE (to assemble overlap matrix) */
	virtual void calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* node_internal_indices, int &n) const;
	virtual const int* get_node_sparsity_pattern(int&) const;
	virtual string get_unique_identifier() const { return string("ClassOverlap"); }
	const CSRMatrix<double>& get_csr_matrix() { return fem_assembler->get_csr_matrix(); }
private:
	FEMSolverME<double>* fem_assembler;
	static const int sparsity_pattern[2];		
};


/** regroup, resort, reorthogonalize bands -> warning, object sets bc to geometry */
class DeCross {
public:	
	DeCross(Geometry& geometry, MaterialDatabase& matdb);	
	DeCross();
	~DeCross();
	void resort(Bandstructure<complex<double> >& bands) const;		
private:
/*	void swap_bands(
		Bandstructure<complex<double> >& bands, 
		unsigned int v1idx, 
		unsigned int v2idx, 
		unsigned int kidx
	);
	void swap_bands(
		Bandstructure<complex<double> >& bands, 
		unsigned int v1idx, 
		unsigned int v2idx, 
		unsigned int ksidx,
		unsigned int keidx
	);
	*/	 
	cplx get_overlap(const EigenSolution<cplx>& lhs, const EigenSolution<cplx>& rhs) const; 
	Overlap* overlap;

};


}

#endif /*OVERLAP_H_*/
