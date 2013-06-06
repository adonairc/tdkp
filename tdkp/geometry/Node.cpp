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

#include "tdkp/geometry/Node.h"
#include "tdkp/common/Exception.h"
#include "tdkp/common/Logger.h"


using namespace tdkp;

namespace tdkp {

ostream& operator<<(ostream& stream, const Node & node) {
	const double* coords = node.get_coords(); 
	stream << "node " << node.get_index_global() << ", vidx = " 
	       << node.get_index_vertex() << ",  (internal: "
	       << node.get_index_internal() << ")" << " <";
	for(int ii = 0; ii < node.get_dimension(); ii++) {
		if(ii > 0) { stream << ", "; }
		stream << coords[ii];	
	}
	stream << "> ";
	return stream;		
}

void Node::print() const {
	cout << *this;	
}

}

/** 3d node constructor
 * 
 * creates a 3d node
 * 
 * @param x x - coordinate
 * @param y y - coordinate
 * @param z z - coordinate
 */ 
Node::Node(double x_,double y_, double z_)
: index_internal(0),
  index_global(-1),
  index_vertex(-1),
  dimension(3),
  location(location_undefined),
  node_value_type(FunctionValue)
{
	this->coord[0]       = x_;
	this->coord[1]       = y_;
	this->coord[2]       = z_; 
}
/** 3d node constructor */
Node::Node(unsigned int index_global_, double x_,double y_, double z_) 
: index_internal(0),
  index_global(index_global_),
  index_vertex(-1),
  dimension(3),
  location(location_undefined),
  node_value_type(FunctionValue) 
{	
	this->coord[0] = x_;
	this->coord[1] = y_;
	this->coord[2] = z_;		
}

/** 2d node constructor */
Node::Node(double x, double y) 
: index_internal(0),
  index_global(-1),
  index_vertex(-1),
  dimension(2),
  location(location_undefined),
  node_value_type(FunctionValue)
{
	this->coord[0]       = x;
	this->coord[1]       = y;
	this->coord[2]       = 0;			
}
/** 2d node constructor */
Node::Node(unsigned int index_global_, double x, double y) 
: index_internal(0),
  index_global(index_global_),
  index_vertex(-1),
  dimension(2),
  location(location_undefined),
  node_value_type(FunctionValue) 
{
	this->coord[0]       = x;
	this->coord[1]       = y;	
	this->coord[2]       = 0;		
}
/** 1d node constructor */
Node::Node(double x)
: index_internal(0),
  index_global(-1),
  index_vertex(-1),
  dimension(1),
  location(location_undefined),
  node_value_type(FunctionValue)
{
	this->coord[0]       = x;
	this->coord[1]       = 0;
	this->coord[2]       = 0;
}
/** 1d node constructor */
Node::Node(unsigned int index_global_, double x)
: index_internal(0),
  index_global(index_global_),
  index_vertex(-1),
  dimension(1),
  location(location_undefined),
  node_value_type(FunctionValue) 
{
	this->coord[0]       = x;
	this->coord[1]       = 0;
	this->coord[2]       = 0;
}

Node::Node(const Node& copy) 
: index_internal(copy.index_internal),
  index_global(copy.index_global),
  index_vertex(copy.index_vertex),
  dimension(copy.dimension),
  location(copy.location),
  node_value_type(copy.node_value_type)
{
	for(unsigned int ii = 0; ii < 3; ii++) {
		coord[ii] = copy.coord[ii];	
	}
}

Node::~Node() {	
	// rien du tout		
}

void Node::set_location(NodeLocation loc) {
	this->location = loc;	 	
}

NodeLocation Node::get_location() const {
	return this->location;	 	
}

/** set primary node status
 * 
 * primary nodes are nodes corresponding to the value of the solution,
 * while other nodes correspond to some derivatives 
 */
void Node::set_value_type(NodeValueType node_value_type_) {
	node_value_type = node_value_type_;	
}

void Node::set_index_internal(int index_internal_) {
	index_internal = index_internal_;
}

void Node::set_index_global(unsigned int index_global_) {
	index_global = index_global_; 	
}

/** set vertex index of node (or -1 if it isnt a vertex) */
void Node::set_index_vertex(int index_vertex_) {
	index_vertex = index_vertex_;	
}

double Node::distance(const Node &a, const Node &b) {
	double dist = 0;
	for(short ii = 0; ii < a.dimension; ii++) {
		dist += (a.coord[ii] - b.coord[ii]) * (a.coord[ii] - b.coord[ii]);	
	}	
	return sqrt(dist);
}

unsigned int Node::get_num_contributions() const {
	if(get_index_vertex() == -1) {
		return vertex_contributions.size();
	} else {
		return 1;	
	}
}

void Node::get_contribution(unsigned int ii, int& contrib_index_vertex, double& contrib) const {
	// its a vertex
	if(get_index_vertex() != -1) {
		TDKP_BOUNDS_ASSERT(ii == 0, "");
		contrib_index_vertex = get_index_vertex();
		contrib = 1.0;
	} else {
		// its a node
		TDKP_BOUNDS_ASSERT(ii < vertex_contributions.size(), "");
		contrib_index_vertex = vertex_contributions[ii].index_vertex;
		contrib = vertex_contributions[ii].contrib;
	}	
}

void Node::set_contribution(const int contrib_index_vertex, const double& contrib) {
	TDKP_ASSERT(get_index_vertex() == -1, "you can not set contributions to a vertex node!");
	Contribution tmp;
	tmp.index_vertex = contrib_index_vertex;
	tmp.contrib      = contrib;
	for(unsigned int ii = 0; ii < vertex_contributions.size(); ii++) {
		TDKP_ASSERT(tmp.index_vertex != vertex_contributions[ii].index_vertex, "vertex contribution of vertex " << tmp.index_vertex << " is already set!");	
	}
	vertex_contributions.push_back(tmp);
}
