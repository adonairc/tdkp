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

#ifndef GRIDINTERPOLATOR_H_
#define GRIDINTERPOLATOR_H_

#include "tdkp/common/Logger.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/main/Fields.h"

namespace tdkp {

/** interpolator to interpolate from one grid to the another */
class GridInterpolator {
public:

	static const int both_interpolations    = 0;
	static const int nodal_interpolations   = 1;
	static const int element_interpolations = 2;

	GridInterpolator(
		const Geometry& source_, 
		const Geometry& target_, 
		const RMatrix<double>& target_rotation_matrix, 
		const Vector3D& target_translation,
		const int interpolation_flag_ = both_interpolations
	);
	virtual ~GridInterpolator();
	
	template<class T>
	void interpolate(const NodeData<T>& source_data, NodeData<T>& target_data) const;
	template<class T>
	void interpolate(const ElementData<T>& source_data, ElementData<T>& target_data) const;	
	void interpolate(const StrainField& strain_field, StrainField& target_data) const; 

private:
	void init(const RMatrix<double>& target_rotation_matrix, const Vector3D& target_translation);
	
	const Geometry& source;
	const Geometry& target;
	
	struct Contribution {				
		unsigned int source_index_global;
		double contribution;
	};
	
	vector<list<Contribution> > node_contributions;
	vector<list<Contribution> > element_contributions;
	const int interpolation_flag; 
			
};

template<class T>
void GridInterpolator::interpolate(
	const NodeData<T>& source_data, NodeData<T>& target_data
) const {
	
	TDKP_ASSERT(source_data.get_length() == (signed)source.get_num_nodes(), "the data you provided for interpolation has length " << source_data.get_length() << " while the geometry has " << source.get_num_nodes() << " nodes. this will not work ...");
	TDKP_ASSERT(interpolation_flag == both_interpolations || interpolation_flag == nodal_interpolations,"");
	TDKP_BOUNDS_ASSERT(node_contributions.size() == target.get_num_nodes(), "");
	
	// init target
	target_data.set_length(node_contributions.size(), source_data.get_num_data_per_node());
	// for every node in target grid
	for(unsigned int ii = 0; ii < node_contributions.size(); ii++) {			
		list<Contribution>::const_iterator it  = node_contributions[ii].begin();
		list<Contribution>::const_iterator end = node_contributions[ii].end();		 
		// for every dataset
		while(it != end) {
			for(int xx = 0; xx < source_data.get_num_data_per_node(); xx++) {
				double val = target_data.get_node_value(ii, xx);	
				val += (*it).contribution * source_data.get_node_value((*it).source_index_global, xx);
				target_data.set_node_value(ii, xx, val);
			}
			it++;
		}
	}
		
}

template<class T>
void GridInterpolator::interpolate(
	const ElementData<T>& source_data, ElementData<T>& target_data
) const {
	
	TDKP_ASSERT(source_data.get_length() == (signed)source.get_num_elements(), "the data you provided for interpolation has length " << source_data.get_length() << " while the geometry has " << source.get_num_elements() << " elements. this will not work ...");
	TDKP_ASSERT(interpolation_flag == both_interpolations || interpolation_flag == element_interpolations,"");
	TDKP_BOUNDS_ASSERT(element_contributions.size() == target.get_num_elements(), "");
	
	// init target
	target_data.set_length(element_contributions.size(), source_data.get_num_data_per_element());
	// for every node in target grid
	for(unsigned int ii = 0; ii < element_contributions.size(); ii++) {			
		list<Contribution>::const_iterator it  = element_contributions[ii].begin();
		list<Contribution>::const_iterator end = element_contributions[ii].end();		 
		// for every dataset
		while(it != end) {
			for(int xx = 0; xx < source_data.get_num_data_per_element(); xx++) {
				double val = target_data.get_element_value(ii, xx);	
				val += (*it).contribution * source_data.get_element_value((*it).source_index_global, xx);
				target_data.set_element_value(ii, xx, val);
			}
			it++;
		}
	}
		
}


}

#endif /*GRIDINTERPOLATOR_H_*/
