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

#ifndef CLCINGREDIENTS_H_
#define CLCINGREDIENTS_H_

#include "tdkp/common/all.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/utilities/MatrixElements.h"
#include "tdkp/coulomb/CoulombMatrixElement.h"
#include "tdkp/clc/Splines.h"
#include "tdkp/clc/CLCCoulombLambda.h"

namespace tdkp {

/** ingredients used for the clc optics calculation */
class CLCIngredients {
public:
	virtual ~CLCIngredients() {}
	
	// --------------------------------------------
	// data access functions
	// --------------------------------------------
	const DomainMaster& get_domain() const;
	const BandstructureDomain<complex<double> >& get_cb_bands() const; 
	const BandstructureDomain<complex<double> >& get_vb_bands() const;
	const MatrixElements& get_matrix_elements() const;
	const CoulombMatrixElementData& get_coulomb_matrix_elements() const;
	const double& get_homogeneous_broadening() const;		
	const double& get_optical_refractive_index() const;
	const double& get_static_dielectric_constant() const;
	const double& get_geometry_output_division_factor() const;	
	
protected:
	CLCIngredients(
		const BandstructureDomain<complex<double> >& cb_bands, 
		const BandstructureDomain<complex<double> >& vb_bands,
		const MatrixElements& matrix_elements,
		const CoulombMatrixElementData& coulomb_matrix_elements,
		const double& homogeneous_broadening,		
		const double& optical_refractive_index,
		const double& static_dielectric_constant,
		const double& geometry_output_division_factor
	);
private:		
	const BandstructureDomain<complex<double> >& cb_bands; 
	const BandstructureDomain<complex<double> >& vb_bands;
	const MatrixElements& matrix_elements;
	const CoulombMatrixElementData& coulomb_matrix_elements;
	const double homogeneous_broadening;		
	const double optical_refractive_index;
	const double static_dielectric_constant;
	const double geometry_output_division_factor;
		
};

/** ingredients for the radial bulk / well and the standard wire case (with 1D splines) */
class CLCIngredients1DR : public CLCIngredients {
	
public:
	CLCIngredients1DR(
		const BandstructureDomain<complex<double> >& cb_bands, 
		const BandstructureDomain<complex<double> >& vb_bands,
		const MatrixElements& matrix_elements,
		const CoulombMatrixElementData& coulomb_matrix_elements,
		const double& homogeneous_broadening,		
		const double& optical_refractive_index,
		const double& static_dielectric_constant,
		const double& geometry_output_division_factor
	);
	virtual ~CLCIngredients1DR();
	
	const Spline1D& get_cb_band(unsigned int cb_idx) const;
	const Spline1D& get_vb_band(unsigned int vb_idx) const;	
	const CLCLinearCurve& get_momentum_matrix_element(
		unsigned int cb_idx, 
		unsigned int vb_idx, 
		unsigned int dir, 
		unsigned int value_index
	) const;
	// Lambda_mm, radial classes also contain the radial integration
	const CLCRadialCoulombLambda& get_coulomb_matrix_element_cb_cb(
		unsigned int cb_idx
	) const;
	// Lambda_nm, radial classes also contain the radial integration	
	const CLCRadialCoulombLambda& get_coulomb_matrix_element_vb_cb(
		unsigned int vb_idx,
		unsigned int cb_idx
	) const;
	// Lambda_nn, radial classes also contain the radial integration
	const CLCRadialCoulombLambda& get_coulomb_matrix_element_vb_vb(
		unsigned int vb_idx
	) const;	
	
private:

	unsigned int get_storage_idx(		
		unsigned int cb_idx, 
		unsigned int vb_idx, 
		unsigned int dir, 
		unsigned int value_index
	) const;

	vector<Spline1D*>         curve_cb_bands;
	vector<Spline1D*>         curve_vb_bands;	
	vector<CLCLinearCurve*>   curve_matrix_elements;
	vector<CLCRadialCoulombLambda* > cb_cb_coulomb_matrix_elements;
	vector<CLCRadialCoulombLambda* > vb_vb_coulomb_matrix_elements;
	vector<CLCRadialCoulombLambda* > vb_cb_coulomb_matrix_elements;
		
	CLCRadialCoulombLambda* create_clc_coulomb_lambda_object(
		const vector<double>& q_values, 
		const vector<cplx>& data, 
		unsigned int kspace_dimension
	) const;		
		
};



}

#endif /*CLCINGREDIENTS_H_*/
