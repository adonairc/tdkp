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

#ifndef FIELDS_H_
#define FIELDS_H_

#include <vector>

#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/geometry/Element.h"
#include "tdkp/common/Vector3D.h"


namespace tdkp {

/** indexes for accessing strain components in strain tensor */
typedef short StrainComponentIndex ;
const StrainComponentIndex iexx = 0;
const StrainComponentIndex ieyy = 1;
const StrainComponentIndex iezz = 2;
const StrainComponentIndex iexy = 3;
const StrainComponentIndex iexz = 4;
const StrainComponentIndex ieyz = 5;
const StrainComponentIndex ieyx = 3;
const StrainComponentIndex iezx = 4;
const StrainComponentIndex iezy = 5;

class KPMatrixBase;

/** strain tensor definition 
*
* note, the strain is NOT saved in Voigt's notation. i just save all components 
* of the upper triangle.
*/
class StrainTensor {
	friend class KPMatrixBase;
	friend class StrainField;
public:
	StrainTensor() { 
		for(StrainComponentIndex ii = 0; ii < 6; ii++) data[ii] = 0.0e0; 
	}
	
	StrainTensor(const double data_[]) { 
		for(StrainComponentIndex ii = 0; ii < 6; ii++) this->data[ii] = data_[ii];
	}
	StrainTensor(const StrainTensor& copy);
	double get(StrainComponentIndex ieij) const { return this->data[ieij]; }
	double get(short ii, short jj) const;
	void set(StrainComponentIndex ieij, double value) { this->data[ieij] = value; }
	void set_all(double data_[]) {
		for(StrainComponentIndex ii = 0; ii < 6; ii++) this->data[ii] = data_[ii];
	}
	void rotate(const RMatrix<double>& rot_matrix);
	const StrainTensor& operator+=(const StrainTensor& rhs);
	StrainTensor operator/(const double& rhs) const;
	StrainTensor operator*(const double& rhs) const; 	

protected:

	double data[6];	//!< strain tensor, see indices above. NO engineering strains (i.e. offdiag components are not saved with the additional factor of two.
};


/** potential energy field */
class PotentialEnergyField : public StdNodeData<double> {		
public:
	PotentialEnergyField(unsigned int nnodes);	
	virtual ~PotentialEnergyField() {}
	/** set potential energy */
	void set(const vector<double>& nodal_potential);
	
};

/** strain field, defined constant on an element */
class StrainField : public ElementData<double> {
	
public:	
	StrainField(const Geometry& geometry);
	StrainField(const char* binary_file);
 	
	virtual ~StrainField() {}
		
	// ---------------------------------------------------------
	// functions for acessing the strain averaged over element
	// ---------------------------------------------------------		
	StrainTensor& operator[](unsigned int elem_idx);
	const StrainTensor& get(unsigned int elem_idx) const;
	StrainTensor& get(unsigned int elem_idx);
	
	// ---------------------------------------------------------
	// functions for acessing the strain on the elements nodes
	// ---------------------------------------------------------
	/** return true if strain is constant over element */
	bool constant_strain() const { return nodal_strains.size() == 0; }
	/** return strain at node in designated element */
	StrainTensor& get_nodal_strain(unsigned int elem_idx, unsigned int node_idx);
	/** return strain at node in designated element */
	const StrainTensor& get_nodal_strain(unsigned int elem_idx, unsigned int node_idx) const;
	
	// ---------------------------------------------------------
	// functions for object manipulation 
	// ---------------------------------------------------------
	void rotate(const RMatrix<double>& rot_matrix);	
	void write_binary(const char* filename) const;
	void read_binary(const char* filename); 
	
	// ---------------------------------------------------------
	// inhertited virtual functions from ElementData, used for 
	// plotting output, DISABLED FOR READING
	// ---------------------------------------------------------
	const double& get_element_value(unsigned int eidx, unsigned int eqidx) const;
	double&       get_element_value(unsigned int eidx, unsigned int eqidx);
	void   		  set_element_value(unsigned int eidx, unsigned int eqidx, double value);
	int           get_length() const { return this->element_strains.size(); }
	void          set_length(int length, int num_data_per_element);
	string        get_identifier(unsigned int eqidx) const;
			
private:

	static const string binary_header;
	static const string binary_footer;

	vector<StrainTensor>          element_strains;
	vector<vector<StrainTensor> > nodal_strains;		
};




} // end namespace tdkp

#endif /*FIELDS_H_*/
