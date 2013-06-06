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

#include "GridInterpolator.h"

namespace tdkp {



/** initialize grid interpolator for mapping data from source grid to target grid
 * 
 * a rotation matrix and a translation vector can be passed which will be applied
 * to the coordinates of the target geometry. thereby, lower dimensional 
 * grids can efficently be placed into higher dimensional ones.
 * 
 * @param source_ source grid where the source data is defined on
 * @param target_ target grid where the target data should be defined on
 * @param target_rotation_matrix rotation matrix that should be applied to the target coordinates BEFORE translation
 * @param target_translation translation matrix that should be applied to the target coordinates AFTER rotation
 */
GridInterpolator::GridInterpolator(
	const Geometry& source_, const Geometry& target_, 
	const RMatrix<double>& target_rotation_matrix, const Vector3D& target_translation,
	const int interpolation_flag_
) : source(source_), target(target_), interpolation_flag(interpolation_flag_) 
{
	this->init(target_rotation_matrix, target_translation);
}
GridInterpolator::~GridInterpolator() {
	
}
	

void GridInterpolator::init(
	const RMatrix<double>& target_rotation_matrix, const Vector3D& target_translation
) {
	
	// ---------------------------------------------
	// reject mappings where source is lower dimensional than target
	// ---------------------------------------------
	TDKP_ASSERT(source.get_dimension() >= target.get_dimension(), "source.get_dimension() >= target.get_dimension() failed!");

	// ----------------------------------------------
	// allocate space for target coordinates
	// ----------------------------------------------	
	vector<Vector3D> target_coordinates;
	vector<bool> disabled_element(target.get_num_elements(), false);
					
	// ---------------------------------------------
	// affine mapping of target coordinates
	// ---------------------------------------------
	unsigned int element_offset = 0; 
	unsigned int node_threshold = target.get_num_nodes();
	if(interpolation_flag == both_interpolations) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "GridInterpolator: affine mapping of target coordinates (nodes and\nelement midpoints) Rx + t where t = " << target_translation << " and R = " << target_rotation_matrix);
		target_coordinates.resize(target.get_num_nodes() + target.get_num_elements());
		element_offset = target.get_num_nodes(); // enable element mid points after nodes 
	} else if(interpolation_flag == nodal_interpolations) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "GridInterpolator: affine mapping of target coordinates (nodes only)\nRx + t where t = " << target_translation << " and R = " << target_rotation_matrix);
		target_coordinates.resize(target.get_num_nodes());		
	} else if(interpolation_flag == element_interpolations) {
		TDKP_LOGMSG(LOG_INFO_DEVEL2, "GridInterpolator: affine mapping of target coordinates (elements only)\nRx + t where t = " << target_translation << " and R = " << target_rotation_matrix);
		target_coordinates.resize(target.get_num_elements());
		node_threshold = 0; // disable nodes
	} else {
		TDKP_GENERAL_EXCEPTION("unexpected interpolation flag passed");	
	}	
	

	
	unsigned int disabled = 0;	
	for(unsigned int ii = 0; ii < target_coordinates.size(); ii++) {
		Vector3D x;
		if(ii < node_threshold) {
			x = Vector3D(target.get_node(ii).get_coords());
		} else {
			unsigned int ee = ii - element_offset;
			const Element& elem = target.get_element(ee);
			x = Vector3D(0.0, 0.0, 0.0);
			if(elem.enabled()) {
				// calculate element middle point				
				for(unsigned int nn = 0; nn < elem.get_num_corners(); nn++) {
					x = x + Vector3D(elem.get_corner_node(nn).get_coords());	
				}
				x = x * (1.0 / static_cast<double>(elem.get_num_corners()));
			} else {
				// remember if this is a disabled element
				disabled_element[ee] = true;
				disabled++;	
			} 				
		}
		target_coordinates[ii] = (target_rotation_matrix * x) + target_translation; 	
	}
	
	// ---------------------------------------------
	// for target every node 
	// ---------------------------------------------
	node_contributions.resize(target.get_num_nodes());
	element_contributions.resize(target.get_num_elements());
	unsigned int inside = 0;
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "GridInterpolator: looking for mapped target nodes in source grid");
	#pragma omp parallel 	
	{			
		unsigned int my_inside = 0;
		#pragma omp for schedule(dynamic,10)			
		for(int ii = 0; ii < (signed)target_coordinates.size(); ii++) {						
			// ------------------------------------------------
			// only search for points which are not disabled
			// ------------------------------------------------
			if(ii < (signed)target.get_num_nodes() || !disabled_element[ii - target.get_num_nodes()]) {										
				// ------------------------------------------------
				// loop over all ??? source elements
				// ------------------------------------------------
				const Vector3D& vecii = target_coordinates[ii];
				for(unsigned int ss = 0; ss < source.get_num_elements(); ss++) {
					// ------------------------------------------------
					// check if point is inside element
					// ------------------------------------------------
					const Element& elem = source.get_element(ss);
					if(elem.enabled() && elem.inside_element(vecii(0), vecii(1), vecii(2))) {
						// ------------------------------------------------
						// evaluate form function at that point
						// this gives me the local contribution
						// ------------------------------------------------
						for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
							Contribution contrib;
							contrib.contribution = elem.evaluate_form_function_global(nn, vecii(0), vecii(1), vecii(2));
							if(ii < (signed)target.get_num_nodes()) {
								// node wise data needs to know which node
								contrib.source_index_global = elem.get_node(nn).get_index_global();
								node_contributions[ii].push_back(contrib);
							} else {
								// element wise data needs to know which element (local node is clear by indexing in vector)
								contrib.source_index_global = elem.get_index_global();
								unsigned int ee = ii - target.get_num_nodes();
								element_contributions[ee].push_back(contrib);	
							}						
						}
						my_inside++;
						break;
					}	
				}
			}		
		}
		#pragma omp atomic
		inside += my_inside;		
	}		
	TDKP_LOGMSG(LOG_INFO_DEVEL2, "GridInterpolator: " << inside << " out of " << target_coordinates.size() - disabled << " target points (" << disabled << " are disabled) are inside the source grid");
			
}

void GridInterpolator::interpolate(const StrainField& source_data, StrainField& target_data) const {
	
	TDKP_ASSERT(source_data.get_length() == (signed)source.get_num_elements(), "the source strain field has not the same number of elements as the source geometry has.");
	TDKP_ASSERT(source_data.get_length() == (signed)source.get_num_elements(), "the target strain field has not the same number of elements as the target geometry has.");
	TDKP_ASSERT(interpolation_flag == both_interpolations || interpolation_flag == element_interpolations,"");
	
	if(!target_data.constant_strain()) {
		TDKP_LOGMSG(LOG_WARN, "GridInterpolator: sorry, but the target will have only element-wise constant strain.");
	} 	

	// -----------------------------------------
	// first, get the element data correctly (mid point of target grid)
	// -----------------------------------------
	for(unsigned int ii = 0; ii < element_contributions.size(); ii++) {			
		list<Contribution>::const_iterator it  = element_contributions[ii].begin();		
		list<Contribution>::const_iterator end = element_contributions[ii].end();
		// is there anything for that element?
		if(it != end) {
			// get associated element in source
			const Element& source_elem = source.get_element((*it).source_index_global);
			TDKP_ASSERT((source_elem.get_element_order() == 1 && source_data.constant_strain()) || (source_elem.get_element_order() > 1 && !source_data.constant_strain()), "");
			// source data is constant in the given element
			if(source_elem.get_element_order() == 1) {
				target_data.get(ii) = source_data.get(source_elem.get_index_global());
			} else {										
				// source data is not constant, so interpolate 
				unsigned int nn = 0;
				// zero strain tensor
				StrainTensor ipol_strain;					
				// sum over all strain tensors on nodes in source element
				while(it != end) {			
					ipol_strain += source_data.get_nodal_strain((*it).source_index_global, nn)
					             * (*it).contribution;				
					it++;
					nn++;
				}
				target_data.get(ii) = ipol_strain;
			}
		}
	}	
	// -------------------------------------------------
	// now copy the constant strain to the element nodes
	// -------------------------------------------------
	if(!target_data.constant_strain()) {
		for(unsigned int ee = 0; ee < target.get_num_elements(); ee++) {
			const Element& elem = target.get_element(ee);
			if(elem.enabled()) {
				for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
					target_data.get_nodal_strain(ee, nn) = target_data.get(ee);
				}
			}		
		}
	}		
}


} // end of namespace
