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

#include "Fields.h"

namespace tdkp {

StrainTensor::StrainTensor(const StrainTensor& copy) {
	for(unsigned int ii = 0; ii < 6; ii++) {
		data[ii] = copy.data[ii];	
	}	
}

double StrainTensor::get(short ii, short jj) const {
	TDKP_BOUNDS_ASSERT(ii >= 0 && jj >= 0 && ii < 3 && jj < 3, "ii >= 0 && jj >= 0 && ii < 3 && jj < 3");
	if(ii == jj) {
		return this->data[ii];
	} else { 
		return this->data[2 + ii + jj];
	}
}

/** rotation of the strain
 * 
 * given by e' = R e R^t
 */
void StrainTensor::rotate(const RMatrix<double>& rot_matrix) {
		
	RMatrix<double> strain_tmp(3,3);
	RMatrix<double> strain_res(3,3);
	for(unsigned int ii = 0; ii < 3; ii++) {
		for(unsigned int jj = 0; jj < 3; jj++) {
			strain_tmp(ii,jj) = this->get(ii,jj);	
		}	
	}
	TDKP_BOUNDS_ASSERT(rot_matrix.cols() == rot_matrix.rows() && rot_matrix.rows() == 3, "rot_matrix.cols() == rot_matrix.rows() && rot_matrix.rows() == 3");
	strain_res = rot_matrix * strain_tmp * rot_matrix.get_transpose();
	strain_tmp = strain_res; 
	this->set(iexx, strain_tmp(0,0));
	this->set(ieyy, strain_tmp(1,1));
	this->set(iezz, strain_tmp(2,2));
	this->set(iexy, strain_tmp(0,1));
	this->set(iexz, strain_tmp(0,2));
	this->set(ieyz, strain_tmp(1,2));
	TDKP_BOUNDS_ASSERT(tdkp_math::abs(strain_tmp(0,1) - strain_tmp(1,0)) < 1.0e-10, "strain_tmp(0,1) (" << strain_tmp(0,1) << ") == strain_tmp(1,0) (" << strain_tmp(1,0) << ")");
	TDKP_BOUNDS_ASSERT(tdkp_math::abs(strain_tmp(0,2) - strain_tmp(2,0)) < 1.0e-10, "strain_tmp(0,2) (" << strain_tmp(2,0) << ") == strain_tmp(2,0) (" << strain_tmp(2,0) << ")");
	TDKP_BOUNDS_ASSERT(tdkp_math::abs(strain_tmp(1,2) - strain_tmp(2,1)) < 1.0e-10, "strain_tmp(1,2) (" << strain_tmp(1,2) << ") == strain_tmp(2,1) (" << strain_tmp(2,1) << ")");		 	
		
}


 const StrainTensor& StrainTensor::operator+=(const StrainTensor& rhs) {	
	for(unsigned int ii = 0; ii < 6; ii++) {
		this->data[ii] += rhs.data[ii]; 	
	}
	return *this;
}
StrainTensor StrainTensor::operator/(const double& rhs) const {
	StrainTensor copy(*this);
	for(unsigned int ii = 0; ii < 6; ii++) {
		copy.data[ii] /= rhs; 	
	}
	return copy;	
}
StrainTensor StrainTensor::operator*(const double& rhs) const {
	StrainTensor copy(*this);
	for(unsigned int ii = 0; ii < 6; ii++) {
		copy.data[ii] *= rhs; 	
	}
	return copy;	
}


/** field constructor
 * 
 * @param nnodes number of nodes = number of values in field
 */
PotentialEnergyField::PotentialEnergyField(unsigned int nnodes) { 
	this->set_length(nnodes, 1);
	this->set_identifier(0, "Potential");
}

/** set potential from vector
 * 
 * vector must have length nnodes or we throw an exception
 */
void PotentialEnergyField::set(const vector<double>& nodal_potential) {
	TDKP_ASSERT((signed)nodal_potential.size() == get_length(), "");
	for(unsigned int ii = 0; ii < nodal_potential.size(); ii++) {
		this->get_node_value(ii,0) = nodal_potential[ii];	
	}	
}

/** string used in binary file for checking if data was read consistently */
const string StrainField::binary_header("StrainFieldBegin");
/** string used in binary file for checking if data was read consistently */
const string StrainField::binary_footer("StrainFieldEnd");

/** strain field constructor 
 * 
 * @param geometry geometry object of associated strain field 
 */
StrainField::StrainField(const Geometry& geometry) 
: element_strains(geometry.get_num_elements())
{
	this->num_data_per_element = 6;
	// -------------------------------------------
	// initialize nodal strains if we have higher order elements
	// -------------------------------------------
	unsigned int order = 0;
	for(unsigned int ii = 0; ii < geometry.get_num_elements(); ii++) {
		if(geometry.get_element(ii).get_element_order() > order) {
			order = geometry.get_element(ii).get_element_order(); 	
		}	
	}
	TDKP_ASSERT(order > 0, "");
	if(order > 1) {
		nodal_strains.resize(geometry.get_num_elements());
		for(unsigned int ii = 0; ii < geometry.get_num_elements(); ii++) {
			const Element& elem = geometry.get_element(ii);
			if(elem.enabled()) {
				nodal_strains[ii].resize(elem.get_num_nodes());	
			}	
		}
	}	
} 
  

/** strain field constructor (reads strain data from binary file) */
StrainField::StrainField(const char* binary_file) {
	this->num_data_per_element = 6;
	this->read_binary(binary_file);	
}


/** rotate strain field */
void StrainField::rotate(const RMatrix<double>& rot_matrix) {
	
	// -----------------------------------------
	// rotate element wise strain field
	// -----------------------------------------
	#pragma omp parallel for
	for(int ii = 0; ii < (signed)element_strains.size(); ii++) {
		element_strains[ii].rotate(rot_matrix);		
	}
	// -----------------------------------------
	// rotate nodal strain field
	// -----------------------------------------
	if(!constant_strain()) {
		#pragma omp parallel for
		for(int ii = 0; ii < (signed)nodal_strains.size(); ii++) {
			for(unsigned int jj = 0; jj < nodal_strains[ii].size(); jj++) {
				nodal_strains[ii][jj].rotate(rot_matrix);
			}		
		}		
	}	
} 

/** return strain at node in designated element (excpetion for disabled elements) */
StrainTensor& StrainField::get_nodal_strain(unsigned int elem_idx, unsigned int node_idx) {
	return const_cast<StrainTensor&>(static_cast<const StrainField&>(*this).get_nodal_strain(elem_idx, node_idx));
}
/** return strain at node in designated element */
const StrainTensor& StrainField::get_nodal_strain(unsigned int elem_idx, unsigned int node_idx) const {
	// if strain is constant, return element strain
	if(constant_strain()) {
		return this->get(elem_idx);
	} else {		
		TDKP_BOUNDS_ASSERT(elem_idx < nodal_strains.size(), "");
		TDKP_BOUNDS_ASSERT(node_idx < nodal_strains[elem_idx].size(), "");
		return nodal_strains[elem_idx][node_idx];
	}	
}

/** store strain field in binary file */
void StrainField::write_binary(const char* filename) const {
	
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "StrainField: writing binary strain data to file " << filename);
	
	// -------------------------------------------
	// create binary file
	// -------------------------------------------
	ofstream fout(filename, ios::binary);
	if(!fout) {
		TDKP_GENERAL_EXCEPTION("can not write to binary file " << filename);	
	}
	// ---------------------------------------------
	// write header
	// ---------------------------------------------
	fout.write(binary_header.c_str(), binary_header.size());
	// ---------------------------------------------
	// write element data
	// ---------------------------------------------	
	int size     	 	  = element_strains.size();	
	fout.write((char*)(&size),  sizeof(int));
	for(int ii = 0; ii < size; ii++) {
		fout.write((char*)(&element_strains[ii].data), sizeof(double) * 6);		
	} 
	// ---------------------------------------------
	// write nodal strains
	// ---------------------------------------------
	size = nodal_strains.size();
	fout.write((char*)(&size),  sizeof(int));
	int elem_size;
	for(int ii = 0; ii < size; ii++) {
		// write number of nodes in element
		elem_size = nodal_strains[ii].size();		
		fout.write((char*)(&elem_size),  sizeof(int));
		// write strains
		for(int jj = 0; jj < elem_size; jj++) {
			fout.write((char*)(&nodal_strains[ii][jj].data), sizeof(double) * 6);
		}		
	} 
	// ---------------------------------------------
	// write footer
	// ---------------------------------------------	
	fout.write(binary_footer.c_str(), binary_footer.size());
	fout.close();
	
}

void StrainField::read_binary(const char* filename) {
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "StrainField: reading binary strain data from file " << filename);

	// -------------------------------------------
	// open binary file
	// -------------------------------------------
	ifstream fin(filename, ios::binary);
	if(!fin) {
		TDKP_GENERAL_EXCEPTION("can not read from binary file " << filename);	
	}
	// ---------------------------------------------
	// read header
	// ---------------------------------------------
	char buf[binary_header.size() + 1]; buf[binary_header.size()] = '\0';
	fin.read(buf, binary_header.size());
	TDKP_ASSERT(binary_header == buf, "this is not strain field data! data starts with " << buf);
	
	// ---------------------------------------------
	// read element data
	// ---------------------------------------------	
	int size;	
	fin.read((char*)(&size), sizeof(int));
	TDKP_ASSERT(size > 0, "");
	element_strains.resize(size);
	for(int ii = 0; ii < size; ii++) {
		fin.read((char*)(&element_strains[ii].data), sizeof(double) * 6);		
	} 
	// ---------------------------------------------
	// read nodal strains
	// ---------------------------------------------	
	fin.read((char*)(&size),  sizeof(int));
	if(size > 0) {
		nodal_strains.resize(size);
		int elem_size;
		for(int ii = 0; ii < size; ii++) {
			// read number of nodes in element		
			fin.read((char*)(&elem_size),  sizeof(int));
			nodal_strains[ii].resize(elem_size);
			// read strains
			for(int jj = 0; jj < elem_size; jj++) {
				fin.read((char*)(&nodal_strains[ii][jj].data), sizeof(double) * 6);
			}		
		}
	} 
	// ---------------------------------------------
	// read footer
	// ---------------------------------------------
	char buf2[binary_footer.size() + 1]; buf2[binary_footer.size()] = '\0';	
	fin.read(buf2, binary_footer.size());
	TDKP_ASSERT(binary_footer == buf2, "strange strain field data does not end with correct footer!");
	fin.close();	
		
} 


string StrainField::get_identifier(unsigned int eqidx) const {
	ostringstream sout;
	sout << "ElasticStrainEL";
	switch(eqidx) {
		case iexx:
			sout << "XX";
			break;
		case ieyy:
			sout << "YY";
			break;
		case iezz:
			sout << "ZZ";
			break;
		case iexy:
			sout << "XY";
			break;
		case iexz:
			sout << "XZ";
			break;
		case ieyz:
			sout << "YZ";
			break;			
	}
	return sout.str();
}

StrainTensor& StrainField::operator[](unsigned int elem_idx) {
	TDKP_BOUNDS_ASSERT(elem_idx < this->element_strains.size(),"");
	return this->element_strains[elem_idx];	
}

const StrainTensor& StrainField::get(unsigned int elem_idx) const {
	TDKP_BOUNDS_ASSERT(elem_idx < this->element_strains.size(),"");
	return this->element_strains[elem_idx];	
}

StrainTensor& StrainField::get(unsigned int elem_idx) {
	return const_cast<StrainTensor&>(
		static_cast<const StrainField&>(*this).get(elem_idx)
	);
}


const double& StrainField::get_element_value(unsigned int eidx, unsigned int eqidx) const {
	TDKP_BOUNDS_ASSERT(eidx >= 0 && eidx < this->element_strains.size() && eqidx >= 0 && eqidx < 6, "");
	return this->element_strains[eidx].data[eqidx];	
}

double& StrainField::get_element_value(unsigned int eidx, unsigned int eqidx) {
	return const_cast<double&>(static_cast<const StrainField&>(*this).get_element_value(eidx,eqidx));		
}

void StrainField::set_element_value(unsigned int eidx, unsigned int eqidx, double value) {
	TDKP_GENERAL_EXCEPTION("StrainField::set_element_value is disabled! sorry, but IO for strain field object must be done using the objects write/read binary functions. element data is only interfaced for plotting purposes");
}

void StrainField::set_length(int size, int num_data_per_element) {
	TDKP_GENERAL_EXCEPTION("StrainField::set_length is disabled! sorry, but IO for strain field object must be done using the objects write/read binary functions. element data is only interfaced for plotting purposes");
}


} // end namespace tdkp
