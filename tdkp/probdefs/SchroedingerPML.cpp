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

#include "tdkp/probdefs/SchroedingerPML.h"

namespace tdkp {

PMLStretch::PMLStretch(
	const double& pml_start_, 
	const double& pml_end_,
	const double& alpha_, 
	const double& beta_,
	unsigned int pml_power_
) : pml_start(pml_start_),
	pml_end(pml_end_),
	fact(alpha_, beta_),
	pml_power(pml_power_)
{
	TDKP_ASSERT(pml_start != pml_end, "");
}

complex<double> PMLStretch::evaluate(const double& x) const {

	double tmp = (x-pml_start)/(pml_end - pml_start);
	TDKP_BOUNDS_ASSERT(tmp >= -1.0e-5, "tmp >= -1.0e-5, [" << pml_start << ", " << pml_end << "] for x = " << x);
	TDKP_BOUNDS_ASSERT(tmp <= 1.0001, "");
	double val = 1.0;
	for(unsigned int ii = 0; ii < pml_power; ii++) {
		val *= tmp;
	}
	return (1.0 + fact*val);	
	
}


template<>
void SchroedingerPML<EffectiveMass>::calculate_element_matrices(const Element* elem, cplx* lhs, cplx* rhs, int* node_internal_indices, int &n) const { 

	
	vector<double> crhs(num_equations_per_node*num_equations_per_node*Element::max_num_nodes*Element::max_num_nodes);
	vector<double> clhs(num_equations_per_node*num_equations_per_node*Element::max_num_nodes*Element::max_num_nodes);
	// -----------------------------------------
	// use pml wrapper element if we are in pml region
	// ----------------------------------------- 
	if(elem->enabled() && pml_regions[elem->get_region().get_index_global()]) {
		// -----------------------------------------
		// init wrapper with my element
		// -----------------------------------------
		TDKP_ASSERT(pml_element_wrappers.size() > elem->get_element_unique_type_key(), "");
		TDKP_ASSERT(pml_element_wrappers[elem->get_element_unique_type_key()] != 0, "");		
		ElementPML copy(*pml_element_wrappers[elem->get_element_unique_type_key()], elem);
		for(unsigned int dd = 0; dd < get_geometry().get_dimension(); dd++) {			
			copy.set_pml_stretch(dd,stretch[dd][elem->get_region().get_index_global()]);
		}
		copy.prepare();
		copy.return_real_part();
		base_problem.calculate_element_matrices(&copy, &clhs[0], &crhs[0], node_internal_indices, n);
		// copy real values from decorated class
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] = crhs[ii];
			lhs[ii] = clhs[ii];				
		}
		copy.return_imaginary_part();
		base_problem.calculate_element_matrices(&copy, &clhs[0], &crhs[0], node_internal_indices, n);
		// add imaginary parts
		cplx i(0.0, 1.0);
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] += i * crhs[ii];
			lhs[ii] += i * clhs[ii];				
		}							
	} else {
		// -----------------------------------------
		// no pml, use standard method
		// -----------------------------------------
		base_problem.calculate_element_matrices(elem, &clhs[0], &crhs[0], node_internal_indices, n);
		// copy values from decorated class
		for(unsigned int ii = 0; ii < crhs.size(); ii++) {
			rhs[ii] = crhs[ii];
			lhs[ii] = clhs[ii];				
		}		
	}

}

}
