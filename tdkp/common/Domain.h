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

#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <list>
#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"

namespace tdkp {

class DomainPoint {
public:
	DomainPoint(const double& x);	
	DomainPoint(const double& x, const double& y);
	DomainPoint(const double& x, const double& y, const double& z);	
	// standard copy constructor is o.k.!
	const double& get_coord(unsigned short idx) const { TDKP_BOUNDS_ASSERT(3 > idx, "idx < 3"); return coords[idx]; }
	const double* get_coords() const { return coords; }
	const double& get_weight() const { return weight; }
	double get_coord_abs() const { return sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]); }
	void set_weight(const double& weight_) { weight = weight_; }
	int  get_index() const { return index; }
	void set_index(int index_) { index = index_; }
	
private:
	double coords[3];
	double weight;
	int    index;

};

ostream& operator<<(ostream& out, const DomainPoint& point);


class DomainNode {
	
	friend class DomainNodeFreeze;
	
public:	
	
	virtual ~DomainNode();
	
	// ---------------------------------------
	// general function
	// ---------------------------------------
	virtual void refine() = 0;
	unsigned short get_dimension() const { return dimension; }
	virtual DomainNode* clone() const = 0;
	/** returns true if domain is approximated radially */
	virtual bool radial() const = 0;
	
	// ---------------------------------------
	// tree functions
	// ---------------------------------------
	unsigned int      get_number_of_children() const { return children.size(); }
	DomainNode&       get_child(unsigned int idx) { return const_cast<DomainNode&>(static_cast<const DomainNode&>(*this).get_child(idx)); }
	/** return child */
	const DomainNode& get_child(unsigned int idx) const { TDKP_BOUNDS_ASSERT(idx < children.size(), "idx < children.size()"); return *children[idx]; }
	/** pop child! returns child with largest index and pops pointer from vector. so you have to delete it later! */
	DomainNode* pop_child();
	 	
	// ----------------------------------------
	// leaf functions
	// ----------------------------------------
	/** return true if we are a leaf */
	bool leaf() const      { TDKP_BOUNDS_ASSERT(point != 0 || children.size() > 0, "point != 0 || children.size() > 0"); return point != 0; }
	DomainPoint&       get_point()       { TDKP_ASSERT(leaf(), "stupid man, only leaf has a point!"); return *point; }
	/** get leafs point */
	const DomainPoint& get_point() const { TDKP_ASSERT(leaf(), "stupid man, only leaf has a point!"); return *point; }
	/** return weight of domain node (== sum off all child nodes) */
	const double       get_weight() const { return weight; }
	
	const DomainNode& operator=(const DomainNode& rhs) { TDKP_GENERAL_EXCEPTION("sorry, this is not allowed!"); }
		
	
protected:
	explicit DomainNode(unsigned short dimension);
	explicit DomainNode(const DomainNode& copy);	
	vector<DomainNode*>  children;
	DomainPoint*         point;
	double               weight;
	unsigned short       dimension;
	
};

class DomainMaster {
public:

	typedef vector<DomainPoint*>::iterator point_iterator;
	typedef vector<DomainPoint*>::const_iterator point_const_iterator;
	
	DomainMaster();
	explicit DomainMaster(DomainNode* node);
	explicit DomainMaster(istream& in); // deserialize
	DomainMaster(const DomainMaster& copy);
	~DomainMaster();
	
	void           clean(); 
	void           refine();
	void           update();
	void           reindex_points();
	void           freeze();
	bool           frozen() const { return is_frozen; }	
	void           collapse();
	void           dump() const;
	bool           compare_points(const DomainMaster& domain) const;
	bool           operator==(const DomainMaster& rhs) const;
		
	double         get_total_weight() const;
	unsigned short get_dimension() const { TDKP_ASSERT(this->root_nodes.size() > 0, "empty domain master"); return root_nodes.front()->get_dimension(); }
	bool           radial() const { TDKP_ASSERT(this->root_nodes.size() > 0, "empty domain master"); return root_nodes.front()->radial(); }
	unsigned int   get_number_of_root_nodes() const { return root_nodes.size(); }
	const DomainNode& get_root_node(unsigned int idx) const;
		
	point_iterator begin();
	point_const_iterator begin() const;
	inline point_const_iterator end() const { return domain_points.end(); }
	unsigned int get_number_of_points() const { return domain_points.size(); }
	const DomainPoint& get_point(unsigned int idx) const { TDKP_BOUNDS_ASSERT(idx < domain_points.size(), "idx < domain_points.size()"); return *domain_points[idx]; }
	const DomainPoint& get_first_point() const { return *domain_points.front(); }
	const DomainPoint& get_last_point()  const { return *domain_points.back(); } 
		
	void add_node(DomainNode* node);
	
	void read_binary(istream& in); 		
	void write_binary(ostream& out) const;
	
	const DomainMaster& operator=(const DomainMaster& rhs);	
	
private:
	vector<DomainNode*>   root_nodes;
	vector<DomainPoint*>  domain_points;
	int  node_counter;	
	bool is_frozen;
	
	void update_domain_points(DomainNode& node);
	void collapse_node_tree(DomainNode& node);	
	void collect_all_indexed_points(DomainNode& node, vector<DomainPoint*>& points);
	
};


void create_nonuniform_point_distribution(vector<double>& target, const double& r0, const double& r1, const unsigned int num, const double& factor);
/** create trapezoidal domain */
void create_1D_domain_trapezoidal(DomainMaster& domain, const double& x0, const double& x1, const unsigned int num);
/** create simpson integration domain */
void create_1D_domain_simpson(DomainMaster& domain, const double& x0, const double& x1, unsigned int num);
/** 1d bandstructure domain for quantum wires (weight factor multiplied by 2, x0 & x1 > 0) */
void create_1D_domain_wire_bands(DomainMaster& domain, const double& x0, const double& x1, const unsigned int num);
/** create radial along [1 0] */
void create_2D_domain_radial(DomainMaster& domain, const double& r0, const double& r1, const unsigned int num);
/** create radial along a direction */
void create_3D_domain_radial(DomainMaster& domain, const Vector3D& direction, const double& r0, const double& r1, const unsigned int num);

/** rectangular 2D domain */
class DomainNodeRectangle : public DomainNode {
public:	
	DomainNodeRectangle(const double& llx, const double& lly, const double& urx, const double& ury);
	virtual ~DomainNodeRectangle();	
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return false; }
private:
	DomainNodeRectangle(DomainPoint* point, const double coords[4]);
	DomainNodeRectangle(const double coords[4]);
	explicit DomainNodeRectangle(const DomainNodeRectangle& copy);
	void init(const double coords[4]);
	/** lower left x, lower left y, upper right x, upper right y */
	double coords[4];
};

/** radial approximation of a 2D domain 
 * 
 * an integral over an area
 * int dx dy 
 * written in the radial approximation is
 * 2pi * int r dr
 * where the integral goes from r0 to r1  
 */
class DomainNodeRadialPlane : public DomainNode {
public:
	DomainNodeRadialPlane(const double& dir_x, const double& dir_y, const double& r0, const double& r1);
	virtual ~DomainNodeRadialPlane();	
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return true; }
private:
	DomainNodeRadialPlane(DomainPoint* point, const double& dir_x, const double& dir_y, const double& r0, const double& r1);
	explicit DomainNodeRadialPlane(const DomainNodeRadialPlane& copy);
	void init();
	
	/** direction of radial plane */
	double dir_x, dir_y;
	/** interval on radial vector */
	double r0, r1;	
};

/** singular point with arbitrary weight ... huston, your control! */
class DomainNodeSingularPoint : public DomainNode {
public:
	DomainNodeSingularPoint(const double& weight_); // 0D point (used for Quantum Dot Bandstructure)
	DomainNodeSingularPoint(bool radial, const double& weight_, const double& x);
	DomainNodeSingularPoint(bool radial, const double& weight_, const double& x, const double& y);
	DomainNodeSingularPoint(bool radial, const double& weight_, const double& x, const double& y, const double& z);
	virtual ~DomainNodeSingularPoint();	
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return is_radial; }
	
private:
	bool is_radial;
	explicit DomainNodeSingularPoint(const DomainNodeSingularPoint& copy);				
};

/** cubic 3D domain */
class DomainNodeCuboid : public DomainNode {
public:	
	DomainNodeCuboid(const double& llx, const double& lly, const double& llz, const double& urx, const double& ury, const double& urz);
	virtual ~DomainNodeCuboid();	
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return false; }
		
private:
	DomainNodeCuboid(DomainPoint* point, const double coords[6]);
	DomainNodeCuboid(const double coords[6]);
	explicit DomainNodeCuboid(const DomainNodeCuboid& copy);
	void init(const double coords[6]);
	double coords[6];	
};

/** line 1D domain (just to be consistent ;-)) */
class DomainNodeLine : public DomainNode {
public:
	DomainNodeLine(const double& lx, const double& rx);
	virtual ~DomainNodeLine();	
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return false; }
	
private:
	DomainNodeLine(DomainPoint* point, const double coords[2]);
	DomainNodeLine(const double coords[2]);
	explicit DomainNodeLine(const DomainNodeLine& copy);
	void init(const double coords[2]);
	/** left x, right x */
	double coords[2];	
};




/** frozen arbitrary domain node 
 *
 * when a domain master freezes (so, no refinements or updates are possible),
 * then it collapes replaces all of its domain nodes by this one
 *  
 * the idea is that if i don't refine anymore, i don't need to 
 * know what i am representing. so i can get rid of the details
 * 
 * the advantage: i just have to think about one way of storing
 * bandstructures and don't have to think how to store the different
 * types of domains
 */
class DomainNodeFreeze : public DomainNode {
public:
	explicit DomainNodeFreeze(DomainNode& take_over);
	explicit DomainNodeFreeze(istream& in); 
	virtual ~DomainNodeFreeze();
	
	// ---------------------------------------
	// general function
	// ---------------------------------------
	virtual void refine();
	virtual DomainNode* clone() const;
	virtual bool radial() const { return is_radial; }
	// -----------------------------------------
	// serialization
	// -----------------------------------------
	void write_binary(ostream& out) const;
	
protected:
	explicit DomainNodeFreeze(const DomainNodeFreeze& copy);
	
private:
	void read_binary(istream& in);
	bool is_radial;	
};

class DomainPointContribution {
public:	
	DomainPointContribution() : point_idx(0), contribution(0.0) {}
	DomainPointContribution(unsigned int point_idx_, const double& contribution_)
	  : point_idx(point_idx_), contribution(contribution_) {}
	
	/** point idx in from Domain Master */
	unsigned int point_idx;
	/* usually inside ]0,1] (except for extrapolations) */
	double contribution;		
};

class DomainMap {
public:
	typedef list<DomainPointContribution>::const_iterator PointMapIterator;
	
	DomainMap(const DomainMaster& source_, const DomainMaster& target_);
	
	PointMapIterator begin(unsigned int target_point_idx) {
		TDKP_ASSERT(point_map.size() > target_point_idx, "point_map.size() > target_point_idx");
		return point_map[target_point_idx].begin();	
	}

	PointMapIterator end(unsigned int target_point_idx) {
		TDKP_ASSERT(point_map.size() > target_point_idx, "point_map.size() > target_point_idx");
		return point_map[target_point_idx].end();
	}

	const DomainMaster& get_source_domain() { return source; }
	const DomainMaster& get_target_domain() { return target; }	
			
private:
	const DomainMaster& source;
	const DomainMaster& target;
	vector<list<DomainPointContribution> > point_map;
};


} // end of namespace

#endif /*DOMAIN_H_*/
