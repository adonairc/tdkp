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

#include "tdkp/coulomb/CoulombIntegrator.h"

namespace tdkp {

/** element middle point constructor */
ElementMiddlePoints::ElementMiddlePoints(const Geometry& geometry) 
: coordinates(geometry.get_dimension() * geometry.get_num_elements()),
  dimension(geometry.get_dimension())
{	
	
	// -------------------------------------------
	// calculate mid points of elements
	// -------------------------------------------
	Geometry::element_const_iterator it, end;
	end = geometry.elements_end();
	double tmp[3];
	for(it  = geometry.elements_begin(); it != end; it++) {
		tmp[0] = tmp[1] = tmp[2] = 0.0e0;
		for(unsigned int vv = 0; vv < (*it)->get_num_nodes(); vv++) {
			for(unsigned int dd = 0; dd < dimension; dd++) {
				tmp[dd] += (*it)->get_node(vv).get_coord(dd);
			}	
		}
		for(unsigned int dd = 0; dd < dimension; dd++) {
			coordinates[(*it)->get_index_global() * geometry.get_dimension() + dd] = tmp[dd] / static_cast<double>((*it)->get_num_nodes());
		}			
	}
		
}
template<>
const int CoulombIntegrator<cplx>::sparsity_pattern[]   = {0,0};
template<>
const int CoulombIntegrator<double>::sparsity_pattern[] = {0,0};



}
