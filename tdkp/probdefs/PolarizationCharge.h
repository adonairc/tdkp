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

#ifndef POLARIZATIONCHARGE_H_
#define POLARIZATIONCHARGE_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/geometry/ElementBoundary.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/main/Fields.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/utilities/PiezoElectricTensor.h"

namespace tdkp {

/** surface and volume charges on element boundaries and inside due to changing polarization
 *  
 * the class is intended to calculate the surface charge between elements
 * caused by a change of the polarization vector.
 *  
 * the system does not know anything about crystal orientation. therefore 
 * it is on the user to supply the spontaneous polarization with respect
 * of system orientation!
 * 
 * further, 
 *   
 */  
class PolarizationCharge {
public:	
	virtual ~PolarizationCharge();
	
	/** set rotation matrix (rotates crystal axes on grid axes) */
	void set_rotation(const RMatrix<double>& rotation_matrix_);
	/** calculate surface charge */
	void calculate(); 
	/** set strains */
	void set_strains(const StrainField* strain_field); 
	/** return calculated surface charge density
	 * 
	 * @return surface charge [global_boundary_index][boundary_node_idx]   
	 */
	const StdNodeData<double>& get_surface_charge() const;
	/** return calculated volume charge density (charge for every node in element) */
	const StdElementData<double>& get_volume_charge() const;
	/** map surface charge to volume charge 
	 * 
	 * the surface charge density is multiplied by the surface volume
	 * divided by two, distributed among the elements and divided by element volume
	 * 
	 * @param target_volume_charge std element data to store volume charge 
	 */
	void map_surface_charge_to_volume_charge(const NodeData<double>& surface_charge, ElementData<double>& target_volume_charge) const;
	

protected:
	PolarizationCharge(const Geometry& geometry, MaterialDatabase& material_database);
	virtual PiezoElectricTensor* create_piezo_tensor(const Material& material) const = 0;
	void init_piezo_tensors();
	
private:
	void calculate_element_polarizations();
	void calculate_surface_charges();
	void calculate_volume_charges();
	void init_spont_polarizations();
	unsigned int get_local_node_idx_from_element(const Element& elem, unsigned int node_index_global) const;

	const Geometry&              geometry;
	MaterialDatabase&            material_database;
	vector<Vector3D>             spont_polarizations;
	vector<bool>                 ignored_materials;
	vector<PiezoElectricTensor*> piezo_tensors;
	vector<vector<Vector3D> >    element_nodal_polarizations; // polarization in each element on every node
	StdNodeData<double>          surface_charges;
	StdElementData<double>       volume_charges;
	const StrainField*           strains;	
	RMatrix<double> 			 rotation_matrix;
};

class PolarizationChargeZincBlende : public PolarizationCharge {
public:	
	PolarizationChargeZincBlende(const Geometry& geometry, MaterialDatabase& material_database);
private:	
	virtual PiezoElectricTensor* create_piezo_tensor(const Material& material) const;	
};

class PolarizationChargeWurtzite : public PolarizationCharge {
public:	
	PolarizationChargeWurtzite(const Geometry& geometry, MaterialDatabase& material_database);
private:	
	virtual PiezoElectricTensor* create_piezo_tensor(const Material& material) const;	
};

} // end of namespace

#endif /*POLARIZATIONCHARGE_H_*/
