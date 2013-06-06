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

#include "tdkp/common/Domain.h"
#include "tdkp/common/Logger.h"

namespace tdkp {

/** constructor for 2D domain points */
DomainPoint::DomainPoint(const double& x, const double& y) 
: weight(0.0), 
  index(-1) 
{ 
	coords[0] = x; 
	coords[1] = y; 
	coords[2] = 0.0e0; 
}

/** constructor for 3D domain points */
DomainPoint::DomainPoint(const double& x, const double& y, const double& z) 
: weight(0.0), 
  index(-1) 
{ 
	coords[0] = x; 
	coords[1] = y; 
	coords[2] = z; 
}

/** constructor for 1D domain points */
DomainPoint::DomainPoint(const double& x) 
: weight(0.0), 
  index(-1) 
{ 
	coords[0] = x; 
	coords[1] = 0.0e0; 
	coords[2] = 0.0e0; 
}

ostream& operator<<(ostream& out, const DomainPoint& point) {
	out << "coords = " << point.get_coord(0) << ", " << point.get_coord(1) 
	    << ", " << point.get_coord(2)
	    <<", weight = " << point.get_weight() << ", index = " 
	    << point.get_index();
	return out;	    	
}


/** protected domain node constructor */
DomainNode::DomainNode(unsigned short dimension_)
: point(0),
  weight(0),
  dimension(dimension_)
{
	//TDKP_TRACE("created DomainNode: " << this);
}

/** copy constructor to perform deep copy */
DomainNode::DomainNode(const DomainNode& copy) 
: point(0),
  weight(copy.weight),
  dimension(copy.dimension)
{
	
	//TDKP_TRACE("copy DomainNode: " << this << " from " << &copy);
	// copy point if necessary
	if(copy.point != 0) {
		point = new DomainPoint(*copy.point);	
	} 
	// copy children
	for(unsigned int ii = 0; ii < copy.children.size(); ii++) {
		children.push_back(copy.children[ii]->clone());
	}
}

/** destructor, deletes children and own point */
DomainNode::~DomainNode() {
//	TDKP_TRACE("deleting DomainNode: " << this);
	for(unsigned int ii = 0; ii < children.size(); ii++) {
		delete children[ii];	
	}
	if(point != 0) {
		delete point; point = 0;	
	}	 	
}

/** pop child from children stack
 * 
 * returns the pointer to the domain node which you then have to delete when 
 * you don't need it anymore
 */
DomainNode* DomainNode::pop_child() {
	if(this->children.size() > 0) {
		DomainNode* tmp = this->children.back();
		this->children.pop_back();
		return tmp; 	
	}	
	return 0;
}




// -----------------------------------------------------------
// Domain Line implementation
// -----------------------------------------------------------
/** public constructor for the 1D line domain node */
DomainNodeLine::DomainNodeLine(const double& lx, const double& rx)
: DomainNode(1)
{	
	double coords_[] = {lx, rx};
	point = new DomainPoint((coords_[0] + coords_[1]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeLine::DomainNodeLine(const DomainNodeLine& copy)
: DomainNode(copy) 
{
	for(unsigned short ii = 0; ii < 2; ii++) {
		coords[ii] = copy.coords[ii];
	}		
}

/** creates a deep clone of my node */
DomainNode* DomainNodeLine::clone() const {
	return new DomainNodeLine(*this);	
}

DomainNodeLine::~DomainNodeLine() {
	
}	

/** protected constructor for the line node that gets the domain point if its parent */
DomainNodeLine::DomainNodeLine(DomainPoint* point_, const double coords_[2])
: DomainNode(1) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the 1D line node that creates its own new constructor */
DomainNodeLine::DomainNodeLine(const double coords_[2])
: DomainNode(1) 
{
	point = new DomainPoint((coords_[0] + coords_[1]) / 2.0);
	this->init(coords_);
}	

void DomainNodeLine::init(const double coords_[2]) {	
	for(unsigned short ee = 0; ee < 2; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = tdkp_math::abs((coords[1] - coords[0]));		     
	point->set_weight(weight); 		
}

/** split the domain node into three new nodes */
void DomainNodeLine::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 3 line elements
		double new_coords[2];
		double seg_x;
		for(unsigned short dd = 0; dd < 3; dd++) {
			seg_x = double(dd) / 3.0;
			// calculate new rectangle boundaries
			new_coords[0] = coords[0] * (1.0 - seg_x)
				          + coords[1] * seg_x; 	
			new_coords[1] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
				          + coords[1] * (seg_x + 1.0 / 3.0);					          
			// ----------------------------------------
			// passing my point
			// ----------------------------------------
			if(dd == 1) {
				children.push_back(
					new DomainNodeLine(
						point,
						new_coords
					)
				);					
				point = 0; 
			} else {
				children.push_back(
					new DomainNodeLine(						
						new_coords
					)
				);
			}		
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}
	
}

// -----------------------------------------------------------
// Domain Rectangle implementation
// -----------------------------------------------------------
/** public constructor for the rectangular domain node */
DomainNodeRectangle::DomainNodeRectangle(const double& llx, const double& lly, const double& urx, const double& ury)
: DomainNode(2)
{	
	double coords_[] = {llx, lly, urx, ury};
	point = new DomainPoint((coords_[0] + coords_[2]) / 2.0, (coords_[1] + coords_[3]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeRectangle::DomainNodeRectangle(const DomainNodeRectangle& copy)
: DomainNode(copy) 
{
	for(unsigned short ii = 0; ii < 4; ii++) {
		coords[ii] = copy.coords[ii];
	}		
}

/** creates a deep clone of my node */
DomainNode* DomainNodeRectangle::clone() const {
	return new DomainNodeRectangle(*this);	
}

DomainNodeRectangle::~DomainNodeRectangle() {
	
}	

/** protected constructor for the rectangular node that gets the domain point if its parent */
DomainNodeRectangle::DomainNodeRectangle(DomainPoint* point_, const double coords_[4])
: DomainNode(2) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the rectangular node that creates its own new constructor */
DomainNodeRectangle::DomainNodeRectangle(const double coords_[4])
: DomainNode(2) 
{
	point = new DomainPoint((coords_[0] + coords_[2]) / 2.0, (coords_[1] + coords_[3]) / 2.0);
	this->init(coords_);
}	

void DomainNodeRectangle::init(const double coords_[4]) {	
	for(unsigned short ee = 0; ee < 4; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = tdkp_math::abs((coords[2] - coords[0]) * (coords[3] - coords[1]));	     
	point->set_weight(weight); 		
}

/** split the domain node into nine new nodes */
void DomainNodeRectangle::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 9 rectangles
		double new_coords[4];
		double seg_x, seg_y;
		for(unsigned short dd = 0; dd < 3; dd++) {
			for(unsigned short ee = 0; ee < 3; ee++) {
				seg_x = double(dd) / 3.0;
				seg_y = double(ee) / 3.0;
				// calculate new rectangle boundaries
				new_coords[0] = coords[0] * (1.0 - seg_x)
					          + coords[2] * seg_x; 	
				new_coords[1] = coords[1] * (1.0 - seg_y)
					          + coords[3] * seg_y;
				new_coords[2] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
					          + coords[2] * (seg_x + 1.0 / 3.0);					          
				new_coords[3] = coords[1] * (1.0 - seg_y - 1.0 / 3.0)
					          + coords[3] * (seg_y + 1.0 / 3.0);
				// ----------------------------------------
				// passing my point
				// ----------------------------------------
				if(dd == 1 && ee == 1) {
					children.push_back(
						new DomainNodeRectangle(
							point,
							new_coords
						)
					);					
					point = 0; 
				} else {
					children.push_back(
						new DomainNodeRectangle(						
							new_coords
						)
					);
				}	
			}	
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}	
}




// -----------------------------------------------------------
// implementation of Radial Approximation of Plane Domain
// -----------------------------------------------------------
/** public constructor for the rectangular domain node */
DomainNodeRadialPlane::DomainNodeRadialPlane(const double& dir_x_, const double& dir_y_, const double& r0_, const double& r1_)
: DomainNode(2),
  dir_x(dir_x_),
  dir_y(dir_y_),
  r0(r0_),
  r1(r1_)
{	
	this->init();	
	point = new DomainPoint(((r1 + r0) / 2.0) * dir_x, ((r1 + r0) / 2.0) * dir_y);
	point->set_weight(weight); 	
}

/** copy constructor for the node */
DomainNodeRadialPlane::DomainNodeRadialPlane(const DomainNodeRadialPlane& copy)
: DomainNode(copy),
  dir_x(copy.dir_x),
  dir_y(copy.dir_y),
  r0(copy.r0),
  r1(copy.r1)
{
	
}

/** creates a deep clone of my node */
DomainNode* DomainNodeRadialPlane::clone() const {
	return new DomainNodeRadialPlane(*this);	
}

DomainNodeRadialPlane::~DomainNodeRadialPlane() {
	
}	

/** protected constructor for the radial node that gets the domain point from its parent */
DomainNodeRadialPlane::DomainNodeRadialPlane(DomainPoint* point_, const double& dir_x_, const double& dir_y_, const double& r0_, const double& r1_)
: DomainNode(2),
  dir_x(dir_x_),
  dir_y(dir_y_),
  r0(r0_),
  r1(r1_)
{
	this->init();
	point = point_;	
	point->set_weight(weight);
}

/** calculates the radial weight and rescales the direction */
void DomainNodeRadialPlane::init() {
		
	// 2pi ri dr
	weight = constants::pi * ((r1 * r1) - (r0 * r0));
	
	double normalize = sqrt(dir_x * dir_x + dir_y * dir_y);
	dir_x /= normalize;
	dir_y /= normalize;
		 	
}

/** split the domain node into three new nodes */
void DomainNodeRadialPlane::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 3 line elements
		double new_coords[2];
		double seg_x;
		for(unsigned short dd = 0; dd < 3; dd++) {
			seg_x = double(dd) / 3.0;
			// calculate new rectangle boundaries
			new_coords[0] = r0 * (1.0 - seg_x)
				          + r1 * seg_x; 	
			new_coords[1] = r0 * (1.0 - seg_x - 1.0 / 3.0)
				          + r1 * (seg_x + 1.0 / 3.0);					          
			// ----------------------------------------
			// passing my point
			// ----------------------------------------
			if(dd == 1) {
				children.push_back(
					new DomainNodeRadialPlane(
						point,
						dir_x,
						dir_y,
						new_coords[0],
						new_coords[1]
					)
				);					
				point = 0; 
			} else {
				children.push_back(
					new DomainNodeRadialPlane(						
						dir_x,
						dir_y,
						new_coords[0],
						new_coords[1]
					)
				);
			}		
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}
}

// --------------------------------------------------------
// singular point implementation
// --------------------------------------------------------
/** initialize a singular point for a singular point (0D system) */
DomainNodeSingularPoint::DomainNodeSingularPoint(const double& weight_) 
: DomainNode(0),
  is_radial(false)
{
	weight = weight_;
	point  = new DomainPoint(0.0);
	point->set_weight(weight);
}

/** initialize a singular point for a line (probably for outputting reason, as else you could e.g. miss the zero point ...) */ 
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x)
: DomainNode(1),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x);
	point->set_weight(weight);	
}
/** initialize a singular point for a plane (probably for outputting reason, as else you could e.g. miss the zero point ...) */
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x, const double& y) 
: DomainNode(2),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x,y);
	point->set_weight(weight);	
}

/** initialize a singular point for a space (probably for outputting reason, as else you could e.g. miss the zero point ...) */
DomainNodeSingularPoint::DomainNodeSingularPoint(bool radial_, const double& weight_, const double& x, const double& y, const double& z)
: DomainNode(3),
  is_radial(radial_)
{
	weight = weight_;
	point = new DomainPoint(x,y,z);
	point->set_weight(weight);	
}

/** copy constructor */
DomainNodeSingularPoint::DomainNodeSingularPoint(const DomainNodeSingularPoint& copy)
: DomainNode(copy),
  is_radial(copy.is_radial)
{
}

/** empty destructur */
DomainNodeSingularPoint::~DomainNodeSingularPoint() {}

/** no refinement happens on a singular point */	
void DomainNodeSingularPoint::refine() {}

/** clone node */
DomainNode* DomainNodeSingularPoint::clone() const {
	return new DomainNodeSingularPoint(*this);	
}

// --------------------------------------------------------
// 3D cuboid implementation
// --------------------------------------------------------	
/** public cuboid constructor */
DomainNodeCuboid::DomainNodeCuboid(const double& llx, const double& lly, const double& llz, const double& urx, const double& ury, const double& urz)
: DomainNode(3)
{	
	double coords_[] = {llx, lly, llz, urx, ury, urz};
	point = new DomainPoint((coords_[0] + coords_[3]) / 2.0, (coords_[1] + coords_[4]) / 2.0, (coords_[2] + coords_[5]) / 2.0);
	this->init(coords_);	
}

/** copy constructor for the node */
DomainNodeCuboid::DomainNodeCuboid(const DomainNodeCuboid& copy)
: DomainNode(copy)
{
	for(unsigned short ii = 0; ii < 6; ii++) {
		coords[ii] = copy.coords[ii];
	}
}

/** creates a deep clone of my node */
DomainNode* DomainNodeCuboid::clone() const {
	return new DomainNodeCuboid(*this);	
}

DomainNodeCuboid::~DomainNodeCuboid() {
	
}	

/** protected constructor for the cuboid node that gets the domain point if its parent */
DomainNodeCuboid::DomainNodeCuboid(DomainPoint* point_, const double coords_[6])
: DomainNode(3) 
{
	point = point_;
	this->init(coords_);
}

/** protected constructor for the cuboid node that creates its own new constructor */
DomainNodeCuboid::DomainNodeCuboid(const double coords_[6])
: DomainNode(3) 
{
	point = new DomainPoint((coords_[0] + coords_[3]) / 2.0, (coords_[1] + coords_[4]) / 2.0, (coords_[2] + coords_[5]) / 2.0);
	this->init(coords_);
}	

void DomainNodeCuboid::init(const double coords_[6]) {	
	for(unsigned short ee = 0; ee < 6; ee++) {
		coords[ee] = coords_[ee];
	}
	weight = tdkp_math::abs((coords[3] - coords[0]) * (coords[4] - coords[1]) * (coords[5] - coords[2]));	     
	point->set_weight(weight); 		
}

/** split the domain node into 27 new nodes */
void DomainNodeCuboid::refine() {

	// --------------------------------------
	// if we are a leaf, we refine ourselves
	// --------------------------------------
	if(get_number_of_children() == 0) {
		// subdivide into 27 cubes
		double new_coords[6];
		double seg_x, seg_y, seg_z;
		for(unsigned short dd = 0; dd < 3; dd++) {
			for(unsigned short ee = 0; ee < 3; ee++) {
				for(unsigned short ff = 0; ff < 3; ff++) {
					seg_x = double(dd) / 3.0;
					seg_y = double(ee) / 3.0;
					seg_z = double(ff) / 3.0;
					// calculate new cube boundaries
					new_coords[0] = coords[0] * (1.0 - seg_x)
						          + coords[3] * seg_x; 	
					new_coords[1] = coords[1] * (1.0 - seg_y)				
						          + coords[4] * seg_y;
					new_coords[2] = coords[2] * (1.0 - seg_z)
						          + coords[5] * seg_z;					          
					new_coords[3] = coords[0] * (1.0 - seg_x - 1.0 / 3.0)
						          + coords[3] * (seg_x + 1.0 / 3.0);					          
					new_coords[4] = coords[1] * (1.0 - seg_y - 1.0 / 3.0)
						          + coords[4] * (seg_y + 1.0 / 3.0);
					new_coords[5] = coords[2] * (1.0 - seg_z - 1.0 / 3.0)
						          + coords[5] * (seg_z + 1.0 / 3.0);					          
					// ----------------------------------------
					// passing my point
					// ----------------------------------------
					if(dd == 1 && ee == 1 && ff == 1) {
						children.push_back(
							new DomainNodeCuboid(
								point,
								new_coords
							)
						);					
						point = 0; 
					} else {
						children.push_back(
							new DomainNodeCuboid(						
								new_coords
							)
						);
					}
				}	
			}	
		}
	} else {
		// -----------------------------------------------
		// we are a tree node, so we refine our children
		// -----------------------------------------------
		for(unsigned int ii = 0; ii < children.size(); ii++) {
			children[ii]->refine();
		}	
	}	
}




// --------------------------------------------------
// implementation of freezed node
// --------------------------------------------------
/** take a domain node, steal its point and its weight */
DomainNodeFreeze::DomainNodeFreeze(DomainNode& take_over) 
: DomainNode(take_over.dimension),
  is_radial(false)
{
	TDKP_ASSERT(take_over.leaf(), "Domain Freezing does only work on collapsed trees! every node must be a leaf and the node you passed is not a leaf!");
	TDKP_ASSERT(take_over.point != 0, "damn, the domain node i'm freezing should have a point");
	TDKP_ASSERT(take_over.children.size() == 0, "take_over.children.size() == 0");
	TDKP_ASSERT(tdkp_math::abs(take_over.weight - take_over.point->get_weight()) < 1.0e-12, "tdkp_math::abs(take_over.weight - take_over.point->get_weight()) < 1.0e-12");
	this->weight    = take_over.weight;
	this->point     = take_over.point;
	this->is_radial = take_over.radial();
	take_over.point = 0; // steal his point!
}

/** copy constructor, deep copy */
DomainNodeFreeze::DomainNodeFreeze(const DomainNodeFreeze& copy)
: DomainNode(copy.dimension) 
{
	this->weight    = copy.weight;
	this->point     = new DomainPoint(*copy.point);
	this->is_radial = copy.is_radial;
}

/** initialize node from binary stream */
DomainNodeFreeze::DomainNodeFreeze(istream& in)
: DomainNode(1) // dummy default 
{
	this->read_binary(in);
} 

/** empty destructor */
DomainNodeFreeze::~DomainNodeFreeze() {
	
}

/** no refinement on frozen nodes! */	
void DomainNodeFreeze::refine() {
	TDKP_GENERAL_EXCEPTION("no, you can not refine freezed nodes!");
}

/** clone node */
DomainNode* DomainNodeFreeze::clone() const {
	return new DomainNodeFreeze(*this);	
}

/** write node data to binary stream */
void DomainNodeFreeze::write_binary(ostream& out) const {
	
	TDKP_ASSERT(point != 0, "point != 0");
	TDKP_ASSERT(point->get_index() >= 0, "point->get_index() >= 0");
		
	struct serialdata {
		unsigned int dimension;
		double weight;
		double x, y, z;
		int point_idx;
		bool is_radial;	
	} serial;
	
	serial.dimension = dimension;
	serial.weight    = weight;
	serial.x         = this->point->get_coord(0);
	serial.y         = this->point->get_coord(1);
	serial.z         = this->point->get_coord(2);
	serial.point_idx = this->point->get_index();
	serial.is_radial = is_radial;
	
	// ------------------------------------------------
	// write struct and eigensolutions
	// ------------------------------------------------
	out.write((char*)&serial, sizeof(serialdata));

}
	
/** reads domain node data and creates a point */	
void DomainNodeFreeze::read_binary(istream& in) {
	
	// works only for fresh stuff
	TDKP_ASSERT(point == 0, "point == 0"); 
			
	struct serialdata {
		unsigned int dimension;
		double weight;
		double x, y, z;
		int point_idx;
		bool is_radial;	
	} serial;
	
	in.read((char*)&serial, sizeof(serialdata));
	
	dimension = serial.dimension;
	weight    = serial.weight;
    is_radial = serial.is_radial;
	point     = new DomainPoint(serial.x, serial.y, serial.z);
	point->set_index(serial.point_idx);
	point->set_weight(weight);
	
}	



// --------------------------------------------------------
// DomainMaster implementation
// --------------------------------------------------------
DomainMaster::DomainMaster(DomainNode* node) 
: node_counter(0),
  is_frozen(false)
{	
	this->add_node(node);		
}


/** empty domain master constructor */
DomainMaster::DomainMaster()
: node_counter(0),
  is_frozen(false)
{

} 




// define small comperator
bool index_cmp(DomainPoint* a, DomainPoint* b) {
	return a->get_index() < b->get_index();
}

DomainMaster::DomainMaster(const DomainMaster& copy) 
: root_nodes(copy.root_nodes.size()),
  domain_points(0),
  node_counter(copy.node_counter),  
  is_frozen(copy.is_frozen)
{	
	// clone root nodes and collect points
	for(unsigned int ii = 0; ii < copy.root_nodes.size(); ii++) {
		root_nodes[ii] = copy.root_nodes[ii]->clone();
		collect_all_indexed_points(*root_nodes[ii], domain_points);	
	}
			
	// and sort the points
	sort(domain_points.begin(), domain_points.end(), index_cmp);

	TDKP_ASSERT(domain_points.size() == copy.domain_points.size(), "domain_points.size() == copy.domain_points.size()");

#ifdef DEBUG	
	// validate	
	for(unsigned int ii = 0; ii < domain_points.size(); ii++) {
		TDKP_ASSERT(domain_points[ii]->get_index() == copy.domain_points[ii]->get_index(), "domain_points[ii]->get_index() (" << domain_points[ii]->get_index() << ") == copy.domain_points[ii]->get_index() (" << copy.domain_points[ii]->get_index() << ") at " << ii);
		TDKP_ASSERT(domain_points[ii]->get_weight() == copy.domain_points[ii]->get_weight(), "domain_points[ii]->get_weight() == copy.domain_points[ii]->get_weight() at " << ii);
		TDKP_ASSERT(domain_points[ii]->get_coord(0) == copy.domain_points[ii]->get_coord(0), "domain_points[ii]->get_coord(0) == copy.domain_points[ii]->get_coord(0) at " << ii);		
	} 	
#endif
			
}

/** deep assignment operator */
const DomainMaster& DomainMaster::operator=(const DomainMaster& rhs) {
	if(this == &rhs) {
		return *this;	
	}
	// purge old 
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			TDKP_GENERAL_EXCEPTION("somehow my root node is zero!");	
		}	
	}		
	root_nodes.resize(rhs.root_nodes.size());
  	domain_points.resize(0);
  	node_counter = rhs.node_counter;  
  	is_frozen    = rhs.is_frozen;
	// clone root nodes and collect points
	for(unsigned int ii = 0; ii < rhs.root_nodes.size(); ii++) {
		root_nodes[ii] = rhs.root_nodes[ii]->clone();
		collect_all_indexed_points(*root_nodes[ii], domain_points);	
	}			
	// and sort the points
	sort(domain_points.begin(), domain_points.end(), index_cmp);

	TDKP_ASSERT(domain_points.size() == rhs.domain_points.size(), "domain_points.size() == rhs.domain_points.size()");
	return *this;
}

DomainMaster::~DomainMaster() {
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
	//	TDKP_TRACE("master " << this << " deletes: " << root_nodes[ii]);
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			TDKP_GENERAL_EXCEPTION("somehow my root node is zero!");	
		}	
	}
}

/** reset domain master */
void DomainMaster::clean() {
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
	//	TDKP_TRACE("master " << this << " deletes: " << root_nodes[ii]);
		if(root_nodes[ii] != 0) { 
			delete root_nodes[ii];
		} else {
			TDKP_GENERAL_EXCEPTION("somehow my root node is zero!");	
		}	
	}
	root_nodes.resize(0);
	domain_points.resize(0); // not my pointers!
	node_counter = 0;
	is_frozen = false;	
}
	
void DomainMaster::refine() {
	
	TDKP_ASSERT(!this->frozen(), "can not refine frozen domain!");
	
	int num_points = domain_points.size();
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		root_nodes[ii]->refine();
	}
	update();
	ostringstream sout;	
	sout << "refined master domain from " 
	     << num_points << " to " << domain_points.size();
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	     
}	

void DomainMaster::update() {
	
	TDKP_ASSERT(!this->frozen(), "can not update frozen domain!");
	
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		this->update_domain_points(*root_nodes[ii]);		
	}

}

void DomainMaster::collapse() {
	
	TDKP_ASSERT(!this->frozen(), "can not collapse frozen domain!");
	
	vector<DomainNode*>::iterator it = root_nodes.begin();
	while(it != root_nodes.end()) {		 
		if(!(*it)->leaf()) {
			// store pointer (we will later delete the node)
			DomainNode* pnode = (*it);
			// remove node (this makes the iterator useless)	
			root_nodes.erase(it);
			// collapse its tree 
			this->collapse_node_tree(*pnode);
			// delete it
			delete pnode;
			// start fresh (collapse node tree possibly added new elements)
			it = root_nodes.begin(); 			
		} else {
			it++;
		}
	}
}

void DomainMaster::collapse_node_tree(DomainNode& node) {

	TDKP_ASSERT(!node.leaf(), "!node.leaf()");
	// process 
	while(node.get_number_of_children() > 0) {		
		unsigned int idx = node.get_number_of_children() - 1;
		// -------------------------------------------
		// if child its not a leaf, descend
		// -------------------------------------------		
		if(!node.get_child(idx).leaf()) {
			this->collapse_node_tree(node.get_child(idx));
			// and delete
			delete node.pop_child();
		} else {
			// -------------------------------------------
			// any leafs will be added to global tree
			// -------------------------------------------
			root_nodes.push_back(node.pop_child());			
		}
	}
	
}

/** descends tree and collects pointers of all points which are already indexed */
void DomainMaster::collect_all_indexed_points(DomainNode& node, vector<DomainPoint*>& points) {
	if(node.leaf()) {
		if(node.get_point().get_index() != -1) {
			points.push_back(&node.get_point());
		}
	} else {
		for(unsigned int ii = 0; ii < node.get_number_of_children(); ii++) {
			collect_all_indexed_points(node.get_child(ii), points);
		}
	}		
}

void DomainMaster::update_domain_points(DomainNode& node) {

	// ------------------------------------------------
	// add leafs points if i don't know them already
	// ------------------------------------------------
	if(node.leaf()) {
		DomainPoint& the_point = node.get_point();
		if(the_point.get_index() == -1) {
			the_point.set_index(node_counter++);
			this->domain_points.push_back(&the_point);
		}
	} else {
		for(unsigned int ii = 0; ii < node.get_number_of_children(); ii++) {
			this->update_domain_points(node.get_child(ii));		
		}
	}
	
}

double DomainMaster::get_total_weight() const {
	double weight = 0.0;
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		weight += root_nodes[ii]->get_weight();	
	}	
	return weight;
}

		
DomainMaster::point_iterator DomainMaster::begin() {
	return this->domain_points.begin();	
}
DomainMaster::point_const_iterator DomainMaster::begin() const {
	return this->domain_points.begin();	
}
	 
void DomainMaster::add_node(DomainNode* node) {
	
	TDKP_ASSERT(!this->frozen(), "can not add nodes to frozen domain!");		
	if(this->root_nodes.size() > 0) {
		TDKP_ASSERT(this->root_nodes.back()->get_dimension() == node->get_dimension(), "a domain master may only have domain of equal dimension");
		TDKP_ASSERT(this->root_nodes.back()->radial() == node->radial(), "you can not mix radial and nonradial nodes!");		
	}
	this->root_nodes.push_back(node);
	update();	
}

const DomainNode& DomainMaster::get_root_node(unsigned int idx) const {
	TDKP_ASSERT(idx < root_nodes.size(), "idx < root_nodes.size()");
	return *root_nodes[idx];	
}

/** freeze domain -> replace all nodes by frozen ones 
 * 
 * so we can be serialized!
 */
void DomainMaster::freeze() {
	TDKP_ASSERT(!this->frozen(), "domain master is already frozen!");
	this->update();
	this->collapse();
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		// make a freeze
		DomainNodeFreeze* frosty = new DomainNodeFreeze(*root_nodes[ii]);
		// delete current node
		delete root_nodes[ii];
		// and replace by frozen one
		root_nodes[ii] = frosty;		
	}
	this->is_frozen = true;
}

DomainMaster::DomainMaster(istream& in)
: node_counter(0),
  is_frozen(true)
{
	this->read_binary(in);	
}

/** writes frozen nodes to stream. expects domain master to be frozen! */
void DomainMaster::write_binary(ostream& out) const {
	
	TDKP_ASSERT(this->frozen(), "writing does only work on frozen domains!");
	TDKP_ASSERT(node_counter == (signed)domain_points.size(), "node_counter == domain_points.size()"); 
	
		
	int magic = 1818;
	out.write((char*)&magic, sizeof(int));
	
	// ------------------------------------------------
	// just store the number of nodes 
	// ------------------------------------------------
	out.write((char*)&node_counter, sizeof(int));
	
	// ------------------------------------------------
	// and let the nodes to their business
	// ------------------------------------------------
	for(unsigned int ii = 0; ii < root_nodes.size(); ii++) {
		DomainNodeFreeze* node = dynamic_cast<DomainNodeFreeze*>(root_nodes[ii]);
		TDKP_ASSERT(node != 0, "could not cast node!");
		node->write_binary(out);	
	}
	out.write((char*)&magic, sizeof(int));
				
}

/** read from stream, expects to be an empty class! */
void DomainMaster::read_binary(istream& in) {
	
	TDKP_ASSERT(root_nodes.size() == 0, "reading binary does only work for empty domains! root_nodes.size() == 0");
	node_counter = 0;
  	is_frozen    = true;
	// read magic
	int magic;
	in.read((char*)&magic, sizeof(int));
	TDKP_ASSERT(magic == 1818, "magic key at domain start is wrong!");
	
	// -----------------------------------------------
	// read number of nodes
	// -----------------------------------------------
	in.read((char*)&node_counter, sizeof(int));
	
	ostringstream sout;	
	sout << "DomainMaster: reading " << node_counter << " nodes from stream";	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	
	// -----------------------------------------------
	// read frozen nodes
	// -----------------------------------------------
	for(int ii = 0; ii < node_counter; ii++) {
		root_nodes.push_back(new DomainNodeFreeze(in));
		collect_all_indexed_points(*root_nodes.back(), domain_points);
	}
	
	// -----------------------------------------------
	// sort points
	// -----------------------------------------------
	sort(domain_points.begin(), domain_points.end(), index_cmp);
	
	TDKP_ASSERT(domain_points.size() == root_nodes.size(), "domain_points.size() == root_nodes.size()"); 

	in.read((char*)&magic, sizeof(int));
	TDKP_ASSERT(magic == 1818, "magic key at domain end is wrong!");
				
}

/** sorter according to coordinates */
bool abs_cmp(DomainPoint* a, DomainPoint* b) {
	return a->get_coord_abs() < b->get_coord_abs();
}

/** reindexes all points according to their absolute coordinate value */
void DomainMaster::reindex_points() {
	update();	
	sort(domain_points.begin(), domain_points.end(), abs_cmp);
	// ----------------------------------------------
	// visit all points and reindex
	// ----------------------------------------------
	int index = 0;
	for(unsigned int ii = 0; ii < domain_points.size(); ii++) {
		TDKP_ASSERT(domain_points[ii]->get_index() >= 0, "domain_points[ii]->get_index() >= 0");
		domain_points[ii]->set_index(index++);
	} 
	TDKP_ASSERT(index == node_counter, "index == node_counter");
}


/** checks if points are indexed equally and are close enought to each other to be accepted as the same grid */
bool DomainMaster::compare_points(const DomainMaster& other) const {
	
	if(other.get_number_of_points() != this->get_number_of_points()) {
		return false;	
	}	
	unsigned int dimension = this->get_dimension();
	double tmp;
	for(unsigned int ii = 0; ii < domain_points.size(); ii++) {
		if(domain_points[ii]->get_index() != other.domain_points[ii]->get_index()) {
			return false;	
		}	
		tmp = 0.0;
		for(unsigned int cc = 0; cc < dimension; cc++) {
			tmp += tdkp_math::abs(domain_points[ii]->get_coord(cc) - other.domain_points[ii]->get_coord(cc));			
		}
		if(tmp > 1.0e-10) {
			ostringstream sout;
			sout << "domain points " << ii << " differ by " << tmp;
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
			return false;	
		}
	}	
	return true;	
}

bool DomainMaster::operator==(const DomainMaster& rhs) const {
	return this->compare_points(rhs);	
}

/** dump domain content */
void DomainMaster::dump() const {
	ostringstream out;
	out << "DomainMaster of " << get_dimension() << "D domain.\n";
	if(radial()) {
		out << "We use the radial approximation.\n";	
	}
	out << "There are " << root_nodes.size() << " root nodes and "
	    << domain_points.size() << " indexed points.\n";
	
	for(unsigned int ii = 0; ii < get_number_of_points(); ii++) {
		out << *domain_points[ii] << "\n"; 	
	}		
	cout << out.str();
}



// --------------------------------------------
// domain map implementatiion
// --------------------------------------------
DomainMap::DomainMap(const DomainMaster& source_, const DomainMaster& target_) 
: source(source_), 
  target(target_),
  point_map(target.get_number_of_points())
{
	// -----------------------------------------------
	// check that both domains have the same dimension
	// -----------------------------------------------
	TDKP_ASSERT(source.get_dimension() == target.get_dimension(), "source.get_dimension() == target.get_dimension()");
	TDKP_ASSERT(source.get_number_of_points() > 1, "source.get_number_of_points() > 1");  	
	TDKP_ASSERT(target.get_number_of_points() > 1, "target.get_number_of_points() > 1");	

	// -----------------------------------------------
	// choose strategy
	// -----------------------------------------------
	if((source.radial() && target.radial()) || source.get_dimension() < 2) {
		vector<double> source_place(source.get_number_of_points());
		vector<double> target_place(target.get_number_of_points());
		// --------------------------------------------
		// check if domain points are ordered properly
		// --------------------------------------------
		unsigned int idx = 0;
		DomainMaster::point_const_iterator it = source.begin();
		source_place[idx++] = (*it)->get_coord_abs();		
		for(it++; it != source.end(); it++, idx++) {
			source_place[idx] = (*it)->get_coord_abs();			
			TDKP_ASSERT(source_place[idx - 1] < source_place[idx], "source_place[idx - 1] < source_place[idx]"); 			    
		}
		idx = 0;
		it = target.begin();
		target_place[idx++] = (*it)->get_coord_abs();		
		for(it++; it != target.end(); it++, idx++) {
			target_place[idx] = (*it)->get_coord_abs();			
			TDKP_ASSERT(target_place[idx - 1] < target_place[idx], "target_place[idx - 1] < target_place[idx]"); 			    
		}
		// --------------------------------------------
		// loop over target values and build map
		// --------------------------------------------
		// extrapolation 
		// f(x1 + dx) = f(N1) + (f(N2) - f(N1)) / (N2 - N1) dx
		//            = f(N1) * (1 - dx/(N2 - N1)) + f(N2) * (dx / (N2 - N1))
		idx = 0;
		double dx; 
		double dN2N1 = source_place[1] - source_place[0];
		while(target_place[idx] < source_place[0]) {
			dx = target_place[idx] - source_place[0];
			point_map[idx].push_back(
				DomainPointContribution(0, 1.0 - dx / dN2N1)
			);
			point_map[idx].push_back(
				DomainPointContribution(1, dx / dN2N1)
			);
			idx++;			
		}
		// ---------------------------------------------
		// interpolate in between
		// ---------------------------------------------
		unsigned int source_frame_idx = 0;
		double fraction;
		while(idx < target_place.size() && target_place[idx] <= source_place.back()) {
			
  
			
			// shift source frame if necessary
			while(source_frame_idx < source_place.size() - 1 &&
			      target_place[idx] >= source_place[source_frame_idx + 1]) {
				source_frame_idx++;
			}
			/*
			cout << "target: " << idx << "(" << target_place.size() << ") = " << target_place[idx] << " source: " 
			     << source_frame_idx << "( " << source_place.size() << ") = " << source_place[source_frame_idx] << " next = " << source_place[source_frame_idx + 1] << "\n";*/  
			// catch exact match
			if(source_place[source_frame_idx] == target_place[idx]) {
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx, 1.0)
				);
			} else {
				TDKP_ASSERT(source_frame_idx < source_place.size() - 1, "source_frame_idx < source_place.size() - 1"); 
				// in between
				fraction = (target_place[idx] - source_place[source_frame_idx]) 
				         / (source_place[source_frame_idx + 1] - source_place[source_frame_idx]); 
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx, 1.0 - fraction)
				);				
				point_map[idx].push_back(
					DomainPointContribution(source_frame_idx + 1, fraction)
				);	
			}
			idx++;			       			
		}				
		// same with end points
		// f(xn + dx) = f(xn) + (f(xn) - f(xn-1)) / dxnxn-1 * dx
		//            = f(xn) * (1 + dx / dxnxn-1) - f(xn-1) * dx/dxnxn-1
		idx = target_place.size() - 1;
		dN2N1 = source_place[source_place.size() - 1] - source_place[source_place.size() - 2];
		while(target_place[idx] > source_place.back()) {
			dx = target_place[idx] - source_place.back();
			point_map[idx].push_back(
				DomainPointContribution(source_place.size() - 1, 1.0 + dx / dN2N1)
			);
			point_map[idx].push_back(
				DomainPointContribution(source_place.size() - 2, - dx / dN2N1)
			);
			idx--;							
		}
	} else {
		TDKP_GENERAL_EXCEPTION("sorry, mapping is currently only implemented for radial domains");	
	}
	
}


// --------------------------------------------------------
// radial domain creators
// --------------------------------------------------------
/** 1d bandstructure domain for quantum wires (weight factor multiplied by 2, x0 & x1 > 0) */
void create_1D_domain_wire_bands(DomainMaster& domain, const double& x0, const double& x1, const unsigned int num) {
	TDKP_ASSERT(x0 >= 0.0, "");
	TDKP_ASSERT(x1 >   x0, "");
	TDKP_ASSERT(num > 1, "num > 1");
	domain.clean(); 
	
	const double dx  = (x1 - x0) / static_cast<double>(num - 1);
	const double wdx = tdkp_math::abs(dx);
	
	for(unsigned int ii = 0; ii < num; ii++) {
		double weight = (ii > 0 && ii < num - 1) ? wdx : wdx / 2.0;
		// weight multiplied by a factor of 2 to include summation over -k
		weight *= 2.0;
		// node is marked as radial. radial in 1D means positiv and negativ part is mapped to that		
		domain.add_node(new DomainNodeSingularPoint(true, weight, x0 + dx * ii)); 
	}
	domain.update(); 	
}


/** build trapezoidal 1D domain */
void create_1D_domain_trapezoidal(DomainMaster& domain, const double& x0, const double& x1, const unsigned int num) {
	
	TDKP_ASSERT(num > 1, "num > 1");
	domain.clean(); 
	
	const double dx  = (x1 - x0) / static_cast<double>(num - 1);
	const double wdx = tdkp_math::abs(dx);
	
	for(unsigned int ii = 0; ii < num; ii++) {
		double weight = (ii > 0 && ii < num - 1) ? wdx : wdx / 2.0;
		domain.add_node(new DomainNodeSingularPoint(false, weight, x0 + dx * ii)); 
	}
	domain.update(); 		
}

void create_1D_domain_simpson(DomainMaster& domain, const double& x0, const double& x1, unsigned int num) {

	// -------------------------------------------
	// ensure that num is even
	// -------------------------------------------
	num = num + (num % 2);
	
	// -------------------------------------------
	// build simpson integration 
	// -------------------------------------------
	// f0 + 4f1 + 2f2 + 4f3 + ... + 4fn-1 + fn
	// -------------------------------------------
	double dx = (x1 - x0) / static_cast<double>(num);
	double wdx = tdkp_math::abs(dx);
	// first point
	domain.add_node(new DomainNodeSingularPoint(false, wdx / 3.0, x0));	
	for(unsigned int ii = 0; ii < num - 2; ii += 2) {		
		domain.add_node(new DomainNodeSingularPoint(false, 4.0 * wdx / 3.0, x0 + dx * (ii + 1)));
		domain.add_node(new DomainNodeSingularPoint(false, 2.0 * wdx / 3.0, x0 + dx * (ii + 2)));  			
	}
	domain.add_node(new DomainNodeSingularPoint(false, 4.0 * wdx / 3.0, x1 - dx));
	domain.add_node(new DomainNodeSingularPoint(false, wdx / 3.0, x1));
	domain.update();
}

/** build radial 2D domain in [10] direction */
void create_2D_domain_radial(DomainMaster& domain, const double& r0, const double& r1, const unsigned int num) {

	TDKP_ASSERT(r1 > r0, "r1 > r0");

	domain.clean();

	// ------------------------------------------
	// o.k. we build it that way: 
	// ------------------------------------------
	// num_k_values - 2 is the number of normal radial intervals
	// interval_length = (r1 - r0) / (num k values - 1)	
	// intervals:
	// first  [r0, r0+dr/2]
	// then   [r0 + dr/2 + ii * dr, +dr]
	// last   [r1 - dr/2, r1]
	// ------------------------------------------- 
			
	const double dr  = (r1 - r0) / static_cast<double>(num - 1);
	
	// 2pi ri dr
	double weight = constants::pi * (((r0 + dr / 2.0) * (r0 + dr / 2.0)) - (r0 * r0));
	
	// add first node at r0 as singular point node
	domain.add_node(new DomainNodeSingularPoint(true, weight, r0, 0.0));
	
	// add subsequent nodes
	for(unsigned int ii = 0; ii < num - 2 && num > 2; ii++) {
		domain.add_node(new DomainNodeRadialPlane(1.0, 0.0, r0 + dr/2.0 + dr * ii, r0 + dr/2.0 + dr * ii + dr)); 		
	}
	// add last node
	if(num > 1) {
		weight = constants::pi * ((r1 * r1) - ((r1 - dr / 2.0) * (r1 - dr / 2.0)));
		domain.add_node(new DomainNodeSingularPoint(true, weight, r1, 0));
	}		 
	domain.update();

}

void create_3D_domain_radial(DomainMaster& domain, const Vector3D& direction, const double& r0, const double& r1, const unsigned int num) {
		
	TDKP_ASSERT(r1 >= r0, "r1 >= r0"); 

	Vector3D mydir = direction;
	mydir.normalize();

	domain.clean();
	
	if(num == 1) { 
		domain.add_node(new DomainNodeSingularPoint(true, 1.0e-100, r0 * mydir(0), r0 * mydir(1), r0 * mydir(2)));
	} else {	
		const double dr  = (r1 - r0) / static_cast<double>(num - 1);
		
		double rr = r0 + dr / 2.0;
		double rl = r0;
		double weight = 4.0 / 3.0 * constants::pi * ((rr * rr * rr) - (rl * rl * rl));
																				
		// add first node at r0 as singular point node				
		domain.add_node(new DomainNodeSingularPoint(true, weight,  r0 * mydir(0), r0 * mydir(1), r0 * mydir(2)));
		
		// add subsequent nodes		
		for(unsigned int ii = 0; ii < num - 2; ii++) {
			rl = r0 + dr / 2 + dr * ii;
			rr = r0 + dr / 2 + dr * (ii + 1);
			weight = 4.0 / 3.0 * constants::pi * ((rr * rr * rr) - (rl * rl * rl));	
			double rii = (rr + rl) / 2.0;
			domain.add_node(new DomainNodeSingularPoint(true, weight, rii * mydir(0), rii * mydir(1), rii * mydir(2)));					
		}
		// add last node
		rl = r1 - dr / 2;
		rr = r1;
		weight = 4.0 / 3.0 * constants::pi * ((rr * rr * rr) - (rl * rl * rl));
		domain.add_node(new DomainNodeSingularPoint(true, weight, rr * mydir(0), rr * mydir(1), rr * mydir(2)));
	}
	domain.update();
}

void create_nonuniform_point_distribution(
	vector<double>& target, const double& r0, const double& r1, const unsigned int num, const double& factor
) {
	TDKP_ASSERT(num > 1, "");
	TDKP_ASSERT(r0 < r1, "");
	TDKP_ASSERT(r0 >= 0.0, "");
	TDKP_ASSERT(factor > 0, "");
	target.resize(num);
	
	// -----------------------------------------
	// o.k. here is the theory
	// i have a weight function F(x) (which is the
	// integral of f(x)) 
	// so the volume of the domain is V = F(r1) - F(r0)
	// and the volume of a point is dx = V / (num - 1)
	// therefore, if i have a point x0, then i want
	// F(x1) - F(x0) = dx, so
	// F(x1) = dx + F(x0) and 
	// x1 = Finv(dx + F(x0))
	// -----------------------------------------
	// i used for F(x) = x^(1/n)		
	double volume = pow(r1, 1.0 / factor) - pow(r0, 1.0 / factor); 
	double dx     = volume / (num - 1);	
	target[0] = r0;
	for(unsigned int ii = 1; ii < num; ii++) {
		target[ii] = pow(dx + pow(target[ii - 1], 1.0 / factor), factor);
	}
			
}

} // end of namespace
