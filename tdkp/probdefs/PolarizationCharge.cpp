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

#include "tdkp/probdefs/PolarizationCharge.h"

namespace tdkp {

PolarizationCharge::PolarizationCharge(
	const Geometry& geometry_,
	MaterialDatabase& material_database_
) : geometry(geometry_),
 	material_database(material_database_),
 	spont_polarizations(material_database.get_num_materials()),
 	ignored_materials(material_database.get_num_materials(), false),
 	element_nodal_polarizations(geometry.get_num_elements()),
 	surface_charges(Element::max_num_nodes, geometry_.get_num_boundaries()),
 	volume_charges(1, geometry_.get_num_elements()),
 	strains(0),
 	rotation_matrix(identity_matrix)
{
	// allocate storage for nodal polarization vector
	Vector3D nopol(0.0, 0.0, 0.0);
	for(unsigned int ii = 0; ii < element_nodal_polarizations.size(); ii++) {
		element_nodal_polarizations[ii].assign(
			Element::max_num_nodes, Vector3D(0.0, 0.0, 0.0) 
		);	
	}
	// init spont polarizations (per material)
	this->init_spont_polarizations();
} 	
 	
PolarizationCharge::~PolarizationCharge() {
	
	for(unsigned int ii = 0; ii < piezo_tensors.size(); ii++) {
		delete piezo_tensors[ii];	
	}
	
}

void PolarizationCharge::set_strains(const StrainField* strain_field) {
	strains = strain_field;	
}

void PolarizationCharge::init_spont_polarizations() {
	
	const string psp_keys[] = {
		string("spontaneous_polarization_Psp_xx"),
		string("spontaneous_polarization_Psp_yy"),
		string("spontaneous_polarization_Psp_zz")	
	};
	
	for(int ii = 0; ii < material_database.get_num_materials(); ii++) {
		Vector3D tmp(0.0, 0.0, 0.0);
		// calculate new psp
		for(unsigned int pp = 0; pp < 3; pp++) {
			const Material& mat = *material_database.get_material(ii); 
			if(mat.valid_key(psp_keys[pp].c_str()) && mat.is_set(psp_keys[pp].c_str())) {
				tmp.set(pp, mat.get(psp_keys[pp].c_str()));	
			}	
		}
		// and rotate!
		spont_polarizations[ii] = rotation_matrix * tmp;					
		// check if its gas
		if(material_database.get_material_name(ii) == string("Gas")) {
			ignored_materials[ii] = true;	
		}
	}
}

/** calculate surface charge */
void PolarizationCharge::calculate() {

	ostringstream sout;
	sout << "calculating polarization charge. using rotation matrix\n "
	     << rotation_matrix 
		 << "for a back- and forth transformation. spontaneous polarizations\n"
		 << "and piezo tensors ";
	if(strains == 0) {
		sout << "(skipped, no strains set)";
	} 
	sout << "are given by:\n";
	unsigned int str_length = 0;
	for(int ii = 0; ii < material_database.get_num_materials(); ii++) {
		str_length = max(str_length, static_cast<unsigned int>(material_database.get_material_name(ii).size()));
		if(ignored_materials[ii]) {
			str_length += 10;	
		} 
	}
	sout.setf(ios::left);			
	for(int ii = 0; ii < material_database.get_num_materials(); ii++) {
		ostringstream mname;
		mname << material_database.get_material_name(ii);
		if(ignored_materials[ii]) {
			mname << " (ignored)";	
		}
		sout << setw(str_length) << mname.str() << ": " << spont_polarizations[ii] << "\n";
		if(strains != 0) {
		     sout << *piezo_tensors[ii];
		}
	}
	
	TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
			 	     

	this->calculate_element_polarizations();
	this->calculate_surface_charges();
	this->calculate_volume_charges();
		
} 

/** calculates for each element the nodal polarization
 * 
 * polarization jumps from element to element. we basically have
 * P = Pspont + Ppiezo
 */
void PolarizationCharge::calculate_element_polarizations() {
	
	// --------------------------------------------
	// loop over elements
	// --------------------------------------------
	for(unsigned int ii = 0; ii < geometry.get_num_elements(); ii++) {
		
		const Element& elem = geometry.get_element(ii);		
		if(elem.enabled()) {
			TDKP_BOUNDS_ASSERT(element_nodal_polarizations[ii].size() >= elem.get_num_nodes(), "");
			// ----------------------------------------------
			// initialize nodal polarization with spontaneous 
			// polarizations
			// ----------------------------------------------			
			for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
				TDKP_BOUNDS_ASSERT(spont_polarizations.size() > elem.get_material().get_id(), "");							
				element_nodal_polarizations[ii][nn] = spont_polarizations[elem.get_material().get_id()];
			}
			
			// ----------------------------------------------
			// add piezo polarizations
			// ----------------------------------------------
			// TODO: do this on a per node basis and don't assume constant strain over the element
			// TODO: 1. interpolate strain onto nodes
			//       2. calculate piezo polarization from strain and piezo tensor
			//       3. add this to nodal polarization
			// TODO: REMOVE ABOVE STUFF, I IMPLEMENTED IT ...			
			if(strains != 0) {
				// get piezoelectric tensor
				const PiezoElectricTensor& piezo_tensor = *this->piezo_tensors[elem.get_material().get_id()];
				for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
					Vector3D elem_pol(0.0, 0.0, 0.0);
					// get nodal strain tensor
					const StrainTensor& strain = this->strains->get_nodal_strain(elem.get_index_global(), nn);					
					piezo_tensor.evaluate_polarization(strain, elem_pol);				 
					element_nodal_polarizations[ii][nn] = element_nodal_polarizations[ii][nn] + elem_pol;
				}
			}
		}				
	}	 	
}

void PolarizationCharge::calculate_surface_charges() {
	// -------------------------------------------------
	// loop over all element boundaries
	// check if polarization vector changes
	// calculate drop of polarization between elements
	// and store them into result vector
	// -------------------------------------------------
	for(unsigned int ii = 0; ii < geometry.get_num_boundaries(); ii++) {
		const ElementBoundary& boundary = geometry.get_element_boundary(ii);
		if(boundary.get_num_elements() > 1) {
			 			
			TDKP_ASSERT(boundary.get_num_nodes() > 0, "detected a boundary (gid: " << ii << ") without any node! did you call prepare_boundaries() on the geometry object before?"); 			 			
			 			    
			// ------------------------------------------------
			// get both elements
			// ------------------------------------------------
			const Element& elem0 = boundary.get_element(0);
			const Element& elem1 = boundary.get_element(1);

			// ------------------------------------------------
			// are we ignoring the surface?
			// ------------------------------------------------
			if(ignored_materials[elem0.get_material().get_id()] || ignored_materials[elem1.get_material().get_id()]) {
				for(unsigned int nn = 0; nn < boundary.get_num_nodes(); nn++) {
					surface_charges.set_node_value(ii, nn, 0.0);
				}
			} else {
			
	/*			ostringstream sout;
				sout << "bnd " << ii << " with " << boundary.get_num_nodes() << " nodes has two adjacent elements:\n"
				     << "  elem0: " << elem0.get_index_global() << ", mat0 " 
				     << elem0.get_material().get_id() << "\n"
					 << "  elem1: " << elem1.get_index_global() << ", mat1 " 
				     << elem1.get_material().get_id() << "\n";
		*/		      
				     
				// ------------------------------------------------
				// for every node on the boundary
				// ------------------------------------------------
				for(unsigned int nn = 0; nn < boundary.get_num_nodes(); nn++) {
					double charge = 0.0;
					unsigned int e0_node_index = get_local_node_idx_from_element(elem0, boundary.get_node(nn).get_index_global());
					unsigned int e1_node_index = get_local_node_idx_from_element(elem1, boundary.get_node(nn).get_index_global());				
					// -----------------------------------------------
					// and every dimension
					// -----------------------------------------------
					for(unsigned int dd = 0; dd < elem0.get_dimension(); dd++) {
			/*			if(elem0.get_material().get_id() != elem1.get_material().get_id()) {
							sout << "charge due node " << nn << ", dim " << dd << ": "
							     << ( (element_nodal_polarizations[elem1.get_index_global()][e1_node_index](dd) 
						              - element_nodal_polarizations[elem0.get_index_global()][e0_node_index](dd))
						              * boundary.get_normal()(dd)
						            )
						         << ", normal = " << boundary.get_normal()(dd) << "\n";  																	        
						}
			*/
						// ------------------------------------------------
						// calculate delta P_dd * n_dd
						// ------------------------------------------------										
						charge -= (element_nodal_polarizations[elem1.get_index_global()][e1_node_index](dd) 
						           - element_nodal_polarizations[elem0.get_index_global()][e0_node_index](dd))
						        * boundary.get_normal()(dd);
					}
					// -----------------------------------------------
					// set node value
					// -----------------------------------------------
					surface_charges.set_node_value(ii, nn, charge);
				}
				//TDKP_LOGMSG(LOG_INFO_DEVEL2, sout.str());					
			}			
		}		
	}	
}

/** volume charge is given by - div P */
void PolarizationCharge::calculate_volume_charges() {

	// ------------------------------------------------
	// calculate volume charge 
	// WARNING: WE ASSUME THAT P CHANGES LINEAR
	// AND THEREFORE THE CHARGE DENSITY IS CONSTANT OVER 
	// THE ELEMENT
	// ------------------------------------------------
	vector<double> values(Element::max_num_nodes, 0.0);
	for(unsigned int ii = 0; ii < geometry.get_num_elements(); ii++) {
		const Element& elem = geometry.get_element(ii);		
		if(elem.enabled()) {
			// --------------------------------------------
			// complain if somebody wants to use higher order elements ...
                        // --------------------------------------------
			TDKP_ASSERT(elem.get_element_order() <= 2, "wooohaaaa! i assume that i only have elements of up to second order and therefore use a constant element charge resulting from polarizations. so you either fix me or you comment this assert.");
			// --------------------------------------------
			// calculate - div P 
			// --------------------------------------------
			// for every dim (of P)
			double charge = 0.0;
			for(unsigned int dd = 0; dd < elem.get_dimension(); dd++) {
				// ----------------------------------------------
				// copy nodal P_dd into values
				// ----------------------------------------------
				for(unsigned int nn = 0; nn < elem.get_num_nodes(); nn++) {
					values[nn] = element_nodal_polarizations[elem.get_index_global()][nn](dd);
				}
				// ----------------------------------------------
				// calculate - d/d dd P_dd
				// ----------------------------------------------
				charge -= elem.differentiate_solution(dd, &values[0]);									
			}	
			volume_charges.set_element_value(elem.get_index_global(), 0, charge);		
		}  			
	} 		
}

void PolarizationCharge::map_surface_charge_to_volume_charge(
	const NodeData<double>& passed_surface_charge, ElementData<double>& target_volume_charge
) const {
	
	TDKP_ASSERT(passed_surface_charge.get_length() == (signed)geometry.get_num_boundaries(), "");
	
	// -----------------------------------------
	// reset target volume charge
	// -----------------------------------------
	target_volume_charge.set_length(geometry.get_num_elements(), 1);
	for(unsigned int ii = 0; ii < geometry.get_num_elements(); ii++) {
		target_volume_charge.set_element_value(ii, 0, 0.0);	
	}
	// -----------------------------------------
	// loop over all element boundaries
	// -----------------------------------------	
	for(unsigned int ii = 0; ii < geometry.get_num_boundaries(); ii++) {
		const ElementBoundary& boundary = geometry.get_element_boundary(ii);
		if(boundary.get_num_elements() > 1) {
			TDKP_ASSERT(boundary.get_num_nodes() > 0, "boundary.get_num_nodes() > 0! did you call prepare_boundaries on the geometry object?");
			double charge = 0.0;
			// -----------------------------------------	
			// calculate average nodal surface charge
			// -----------------------------------------	
			for(unsigned int nn = 0; nn < boundary.get_num_nodes(); nn++) {
				charge += passed_surface_charge.get_node_value(ii,nn);	
			}
			charge /= static_cast<double>(boundary.get_num_nodes());
			// -----------------------------------------	
			// integrate over surface
			// -----------------------------------------	
			charge *= boundary.get_surface_volume();
			// -----------------------------------------
			// distribute to elements
			// -----------------------------------------
			for(unsigned int ee = 0; ee < 2; ee++) {
				target_volume_charge.get_element_value(
					boundary.get_element(ee).get_index_global(), 0
				) += charge / 2.0 / boundary.get_element(ee).get_volume();
			}	
		}
	}	
}

/** return calculated surface charge
 * 
 * @return surface charge vector[global_boundary_idx] -> surface charge
 */
const StdNodeData<double>& PolarizationCharge::get_surface_charge() const {
	return surface_charges;		
}

/** return calculated volume charge density */
const StdElementData<double>& PolarizationCharge::get_volume_charge() const {
	return volume_charges;	
}

unsigned int PolarizationCharge::get_local_node_idx_from_element(
	const Element& elem, unsigned int node_index_global
) const {
	for(unsigned int ii = 0; ii < elem.get_num_nodes(); ii++) {
		if(static_cast<unsigned int>(elem.get_node(ii).get_index_global()) == node_index_global) {
			return ii;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("could not find node " << node_index_global << " in element!");
}

 
void PolarizationCharge::set_rotation(const RMatrix<double>& rotation_matrix_) {
		
	this->rotation_matrix = rotation_matrix_;
	
	// ---------------------------------------
	// set rotation matrix to all piezo tensors
	// --------------------------------------- 
	for(unsigned int ii = 0; ii < piezo_tensors.size(); ii++) {
		piezo_tensors[ii]->set_rotation(rotation_matrix);	
	}
	this->init_spont_polarizations();
	
}

/** initialize piezo tensors for every material */
void PolarizationCharge::init_piezo_tensors() {
	
	TDKP_ASSERT(piezo_tensors.size() == 0, "init_piezo_tensors must only be called once!");	
	// ---------------------------------------
	// set rotation matrix to all piezo tensors
	// --------------------------------------- 
	for(int ii = 0; ii < material_database.get_num_materials(); ii++) {
		piezo_tensors.push_back(
			this->create_piezo_tensor(*material_database.get_material(ii))
		);
	}		
}

PolarizationChargeZincBlende::PolarizationChargeZincBlende(
	const Geometry& geometry_, MaterialDatabase& material_database_
) 
: PolarizationCharge(geometry_, material_database_)
{
	// init tensors (must be called in derived class as create piezo tensor is virtual 
	// and exists only after the base constructor is called
	this->init_piezo_tensors();
}

PiezoElectricTensor* PolarizationChargeZincBlende::create_piezo_tensor(const Material& material) const {
	return new PiezoElectricTensorTdZincBlende(material);
}	



PolarizationChargeWurtzite::PolarizationChargeWurtzite(
	const Geometry& geometry_, MaterialDatabase& material_database_
) 
: PolarizationCharge(geometry_, material_database_)
{
	this->init_piezo_tensors();
}
PiezoElectricTensor* PolarizationChargeWurtzite::create_piezo_tensor(const Material& material) const {
	return new PiezoElectricTensorHexagonalWurtzite(material);
}	


} // end of namespace
