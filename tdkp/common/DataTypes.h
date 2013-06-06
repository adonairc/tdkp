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

#ifndef DATATYPES_H_
#define DATATYPES_H_

#include "tdkp/common/all.h"
#include "tdkp/geometry/Geometry.h"

namespace tdkp {

/** base class to hold node based data (based on global id) */
template<class T> 
class NodeData {
		
public:
	NodeData() : num_data_per_node(-1) {};	
	virtual ~NodeData() {};
	virtual const T& get_node_value(unsigned int vidx, unsigned int eqidx = 0) const = 0;
	virtual T&       get_node_value(unsigned int vidx, unsigned int eqidx = 0) = 0;	
	virtual void     set_node_value(unsigned int vidx, unsigned int eqidx, T value) = 0; // not const ref as tcl can not handle const ref ... bahhh
	virtual int      get_num_data_per_node() const { return this->num_data_per_node; }
	virtual int      get_length() const = 0;
	virtual void     set_length(int length, int num_data_per_node) = 0;
	virtual string   get_identifier(unsigned int eqidx) const { return string("anydata_" + eqidx); }
	virtual T*       get_element_data(const Geometry& geom) const throw(Exception*);
	virtual const NodeData<T>& operator=(const NodeData<T>& copy) { TDKP_GENERAL_EXCEPTION("operator= must be implemented in derived class"); }

protected:
	int            num_data_per_node;
	
};

/** base class to hold element based data (based on global id) */
template<class T>
class ElementData {
public:
	ElementData() : num_data_per_element(-1) {} ;
	virtual ~ElementData() {}
	virtual const T& get_element_value(unsigned int eidx, unsigned int eqidx = 0) const = 0;
	virtual T&       get_element_value(unsigned int eidx, unsigned int eqidx = 0) = 0;
	virtual void     set_element_value(unsigned int eidx, unsigned int eqidx, T value) = 0; // not const ref as tcl can not handle const ref ... bahhh
	virtual int      get_length() const = 0; 			
	virtual void     set_length(int length, int num_data_per_element) = 0;
	virtual int      get_num_data_per_element() const { return this->num_data_per_element; }
	virtual string   get_identifier(unsigned int eqidx) const { return string("anydata_" + eqidx); } 
	virtual T*       get_node_data(const Geometry& geom) const throw(Exception*);
	
protected:
	int num_data_per_element;		

};

/** interface class for XY data sets */
template<class T> 
class XYData {
public:
	XYData() {};
	virtual ~XYData() {};	
	virtual int    get_x_length()                const = 0;
	virtual int    get_num_y_sets()              const = 0;
	virtual int    get_num_x_sets()              const = 0; 
	virtual void   get_x(int xidx, vector<T> &x) const = 0;
	virtual void   get_y(int yidx, vector<T>& y) const = 0; 
	virtual string get_x_identifier(int xidx)    const = 0;
	virtual string get_y_identifier(int yidx) 	 const = 0;
};

/** real wrapper */
template<class T> 
class XYDataReal : public XYData<double> {
public:	
	XYDataReal(const XYData<T>& dat);
	~XYDataReal() {};	
	int    get_x_length()                const;
	int    get_num_y_sets()              const; 
	int    get_num_x_sets()              const; 
	void   get_x(int xidx, vector<double> &x) const;
	void   get_y(int yidx, vector<double>& y) const; 
	string get_x_identifier(int xidx)    const;
	string get_y_identifier(int yidx) 	 const;		
private:
	const XYData<T>& data;
};




template<class T> 
class StdElementData : public ElementData<T> {
public:	
	StdElementData() {} ;	
	StdElementData(int num_data_per_element_, int length_);
	StdElementData(const Geometry& geometry, const NodeData<T>& data);
	virtual ~StdElementData() {};
	virtual void     set_from_node_data(const Geometry& geometry, const NodeData<T>& data);
	virtual const T& get_element_value(unsigned int eidx, unsigned int eqidx = 0) const;
	virtual T&       get_element_value(unsigned int eidx, unsigned int eqidx = 0);
	virtual void     set_element_value(unsigned int eidx, unsigned int eqidx, T value);
	virtual int      get_num_data_per_element() const { return this->num_data_per_element; }
	virtual int      get_length() const;
	virtual void     set_length(int length, int num_data_per_node);
	virtual void     set_identifier(unsigned int eqidx, const char* ident);
	virtual string   get_identifier(unsigned int eqidx) const;
	virtual vector<T>& get_data_vector() { return this->data; }
	 
protected:
	void    	   initialize_identifier();
	vector<T>      data;
	vector<string> identifier;		
};

template<class T> 
class StdNodeData : public NodeData<T> {
public:
	StdNodeData() {} ;
	explicit StdNodeData(const NodeData<T>& copy);
	StdNodeData(int num_data_per_node_, int length_);
	StdNodeData(int num_data_per_node_, int length_, T* values);
	StdNodeData(const Geometry& geometry, const ElementData<T>& data);
	virtual ~StdNodeData() {};
	virtual void     set_from_element_data(const Geometry& geometry, const ElementData<T>& data);
	virtual const T& get_node_value(unsigned int vidx, unsigned int eqidx = 0) const;
	virtual T&       get_node_value(unsigned int vidx, unsigned int eqidx = 0);	
	virtual void     set_node_value(unsigned int vidx, unsigned int eqidx, T value);
	virtual int      get_num_data_per_node() const { return this->num_data_per_node; }
	virtual int      get_length() const;
	virtual void     set_length(int length, int num_data_per_node);

	virtual void     set_identifier(unsigned int eqidx, const char* ident);
	virtual string   get_identifier(unsigned int eqidx) const;
	
	virtual vector<T>& get_data_vector() { return this->data; }
	virtual const StdNodeData<T>& operator=(const NodeData<T>& copy);
protected:
	void    	   initialize_identifier();
	vector<T>      data;
	vector<string> identifier;
		
};


template<class T>
const StdNodeData<T>& StdNodeData<T>::operator=(const NodeData<T>& copy) {
	if(this == &copy) {
		return *this;	
	}
	this->set_length(copy.get_length(), copy.get_num_data_per_node());
	unsigned int tlength = copy.get_length();
	for(unsigned int ii = 0; ii < tlength; ii++) {
		for(int nn = 0; nn < this->num_data_per_node; nn++) {
			this->set_node_value(ii, nn, copy.get_node_value(ii,nn));	
		}	
	}
	return *this;	
}

template<class T> 
XYDataReal<T>::XYDataReal(const XYData<T>& dat) : data(dat) {

}

template<class T> 
int XYDataReal<T>::get_x_length() const {
	return data.get_x_length(); 	
}

template<class T> 
int XYDataReal<T>::get_num_y_sets() const {
	return data.get_num_y_sets(); 		
}

template<class T> 
void XYDataReal<T>::get_x(int xidx, vector<double> &x) const {
	vector<T> tmp; 
	x.clear();
	data.get_x(xidx, tmp); 
	for(typename vector<T>::iterator it = tmp.begin(); it != tmp.end(); it++) {
		x.push_back(tdkp_math::only_real(*it)); 
	}	
}

template<class T> 
void XYDataReal<T>::get_y(int yidx, vector<double> &y) const {
	vector<T> tmp;	
	data.get_y(yidx,tmp);
	y.clear();
	for(typename vector<T>::iterator it = tmp.begin(); it != tmp.end(); it++) {
		y.push_back(tdkp_math::only_real(*it)); 
	}
}

template<class T> 
int XYDataReal<T>::get_num_x_sets() const {
	return data.get_num_x_sets();	
}

template<class T> 
string XYDataReal<T>::get_x_identifier(int xidx) const {
	return data.get_x_identifier(xidx); 	
}
template<class T> 
string XYDataReal<T>::get_y_identifier(int yidx) const {
	return data.get_y_identifier(yidx); 	
}



// --------------------------------------------------------
// StdNodeData
// --------------------------------------------------------
template<class T> 
StdNodeData<T>::StdNodeData(const NodeData<T>& copy)  
{
	this->set_length(copy.get_length(),copy.get_num_data_per_node());
	unsigned int tlength = copy.get_length();
	for(unsigned int ii = 0; ii < tlength; ii++) {
		for(int nn = 0; nn < this->num_data_per_node; nn++) {
			data[ii * this->num_data_per_node + nn] = copy.get_node_value(ii,nn);	
		}	
	}		
}
template<class T> 
StdNodeData<T>::StdNodeData(int num_data_per_node_, int length_) {
	this->set_length(length_, num_data_per_node_);	
}

template<class T> 
StdNodeData<T>::StdNodeData(int num_data_per_node_, int length_, T* values) {
	this->set_length(length_, num_data_per_node_);
	const int total_length = length_ * num_data_per_node_;
	for(int vv = 0; vv < total_length; vv++) {
		this->data[vv] = values[vv];
	}
}

/** write element data to node data container */
template<class T>
StdNodeData<T>::StdNodeData(const Geometry& geometry, const ElementData<T>& source) {

	this->set_from_element_data(geometry, source);
	
}


template<class T>
void StdNodeData<T>::set_from_element_data(const Geometry& geometry, const ElementData<T>& source) {
	
	// get data and init
	T* node_data = source.get_node_data(geometry);
	this->set_length(geometry.get_num_nodes(), source.get_num_data_per_element());
	// set identifier	
	for(int vv = 0; vv < this->num_data_per_node; vv++) {
		this->set_identifier(vv, source.get_identifier(vv).c_str());
	}
	// set data
	for(unsigned int ii = 0; ii < data.size(); ii++) {
		this->data[ii] = node_data[ii];
	}
	delete[] node_data; 
	
}

template<class T> 
const T& StdNodeData<T>::get_node_value(unsigned int vidx, unsigned int eqidx) const {	
	TDKP_BOUNDS_ASSERT(vidx >= 0 && vidx < (this->data.size() / this->num_data_per_node), "node idx out of range");
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_node, "equation idx out of range");
	return this->data[vidx * this->num_data_per_node + eqidx];
}

template<class T>
T& StdNodeData<T>::get_node_value(unsigned int vidx, unsigned int eqidx) {
	return const_cast<T&>(static_cast<const StdNodeData&>(*this).get_node_value(vidx,eqidx)); // use const member function and recast	
}

template<class T> 
void StdNodeData<T>::set_node_value(unsigned int vidx, unsigned int eqidx, T value) {	
	TDKP_BOUNDS_ASSERT(vidx >= 0 && vidx < (this->data.size() / this->num_data_per_node), "node idx out of range");
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_node, "equation idx out of range");
	this->data[vidx * this->num_data_per_node + eqidx] = value;
}

template<class T> 
int StdNodeData<T>::get_length() const {
	TDKP_BOUNDS_ASSERT(this->num_data_per_node > 0, "node data is uninitialized");
	return this->data.size() / this->num_data_per_node;	
}

template<class T>
void StdNodeData<T>::set_length(int length, int num_data_per_node_) {
	this->num_data_per_node = num_data_per_node_;
	TDKP_ASSERT(this->num_data_per_node > 0, "num data per node is <= 0");	
	this->data.assign(length * this->num_data_per_node, 0);	
	this->initialize_identifier(); 	
}

template<class T>
void StdNodeData<T>::set_identifier(unsigned int eqidx, const char* ident) {
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_node, "equation idx out of range");		
	this->identifier[eqidx] = string(ident);
}

template<class T>
string StdNodeData<T>::get_identifier(unsigned int eqidx) const {
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_node, "equation idx out of range");	
	return this->identifier[eqidx];	
}

template<class T>
void StdNodeData<T>::initialize_identifier() {
	this->identifier.resize(this->num_data_per_node);
	ostringstream sout;
	for(int ii = 0; ii < this->num_data_per_node; ii++) {
		sout.str(""); sout << "anydata" << ii;
		this->identifier[ii] = sout.str();
	}
}


// -------------------------------------------------------
// StdElementData
// -------------------------------------------------------
template<class T> 
StdElementData<T>::StdElementData(int num_data_per_element_, int length_) {
	this->set_length(length_,num_data_per_element_);
}

template<class T>
StdElementData<T>::StdElementData(const Geometry& geometry_, const NodeData<T>& source) {
	this->set_from_node_data(geometry_, source);
}

template<class T>
void StdElementData<T>::set_from_node_data(const Geometry& geometry, const NodeData<T>& source) {
	// get data and init
	T* element_data = source.get_element_data(geometry);
	this->set_length(geometry.get_num_elements(), source.get_num_data_per_node());
	// set identifier	
	for(int vv = 0; vv < this->num_data_per_element; vv++) {
		this->set_identifier(vv, source.get_identifier(vv).c_str());
	}
	// set data
	for(unsigned int ii = 0; ii < data.size(); ii++) {
		this->data[ii] = element_data[ii];
	}
	delete[] element_data; 	
}

template<class T> 
const T& StdElementData<T>::get_element_value(unsigned int eidx, unsigned int equation_idx) const {	
	TDKP_BOUNDS_ASSERT(eidx >= 0 && eidx < (this->data.size() / this->num_data_per_element), "element idx out of range");
	TDKP_BOUNDS_ASSERT(equation_idx >= 0 && equation_idx < (unsigned)this->num_data_per_element, "equation idx out of range");
	return this->data[eidx * this->num_data_per_element + equation_idx];
}

template<class T>
T& StdElementData<T>::get_element_value(unsigned int eidx, unsigned int equation_idx) {
	return const_cast<T&>(static_cast<const StdElementData&>(*this).get_element_value(eidx,equation_idx)); // use const member function and recast	
}

template<class T> 
void StdElementData<T>::set_element_value(unsigned int eidx, unsigned int equation_idx, T value) {	
	TDKP_BOUNDS_ASSERT(eidx >= 0 && eidx < (this->data.size() / this->num_data_per_element), "element idx out of range");
	TDKP_BOUNDS_ASSERT(equation_idx >= 0 && equation_idx < (unsigned)this->num_data_per_element, "equation idx out of range");	
	this->data[eidx * this->num_data_per_element + equation_idx] = value;
}

template<class T> 
int StdElementData<T>::get_length() const {
	TDKP_BOUNDS_ASSERT(this->num_data_per_element > 0, "element data is uninitialized");
	return this->data.size() / this->num_data_per_element;	
}

template<class T>
void StdElementData<T>::set_length(int length, int num_data_per_element_) {
	this->num_data_per_element = num_data_per_element_;
	TDKP_BOUNDS_ASSERT(this->num_data_per_element > 0, "num data per element is < 0");	
	this->data.assign(length * this->num_data_per_element, 0.0);	
	this->initialize_identifier();
}

template<class T>
void StdElementData<T>::set_identifier(unsigned int eqidx, const char* ident) {
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_element, "equation idx out of range");		
	this->identifier[eqidx] = string(ident);
}

template<class T>
string StdElementData<T>::get_identifier(unsigned int eqidx) const {
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_element, "equation idx out of range");	
	return this->identifier[eqidx];	
}

template<class T>
void StdElementData<T>::initialize_identifier() {
	this->identifier.resize(this->num_data_per_element);
	ostringstream sout;
	for(int ii = 0; ii < this->num_data_per_element; ii++) {
		sout.str(""); sout << "anydata" << ii;
		this->identifier[ii] = sout.str();
	}
}


// TODO: replace T* by vector<T>
template<class T> 
T* NodeData<T>::get_element_data(const Geometry& geom) const throw(Exception*) {
	
	TDKP_ASSERT(num_data_per_node > 0, "NodeData::num_data_per_node must be set");
	TDKP_ASSERT((signed)geom.get_num_nodes()  == this->get_length(), "is this really the corresponding geometry object?");
	
	const int nelem = geom.get_num_elements();		
	T* data         = new T[nelem * this->num_data_per_node];	
	TDKP_POINTER_ASSERT(data);
	for(int ii = 0; ii < nelem * this->num_data_per_node; ii++) {
		data[ii] = 0;	
	}
	for(Geometry::element_const_iterator eit = geom.elements_begin(); eit != geom.elements_end(); eit++) {				
		int offset = (*eit)->get_index_global() * num_data_per_node;
		double divide = static_cast<double>((*eit)->get_num_nodes());			
		// average over node values		
		for(unsigned int vv = 0; vv < (*eit)->get_num_nodes(); vv++) {
			int node_index = (*eit)->get_node(vv).get_index_global();
			for(int ii = 0; ii < this->num_data_per_node; ii++) {
				data[offset + ii] += this->get_node_value(node_index, ii) / divide;
			}
		}
	}
	return data;	
}

// TODO: replace this by vector<T> 
template<class T> 
T* ElementData<T>::get_node_data(const Geometry& geom) const throw(Exception*){
		
	TDKP_ASSERT(get_num_data_per_element() > 0, "ElementData::num_data_per_element must be set");
	TDKP_ASSERT((signed)geom.get_num_elements() == this->get_length(), "is this really the corresponding geometry object?");
	
	// build new data structures
	const int nnode = geom.get_num_nodes();
	T* data           = new T[nnode * get_num_data_per_element()];
	int* track        = new int[nnode];
	TDKP_POINTER_ASSERT(data);
	for(int ii = 0; ii < nnode; ii++) {
		track[ii] = 0;
		for(int jj = 0; jj <  get_num_data_per_element(); jj++) {
			data[ii * get_num_data_per_element() + jj] = 0;	
		}
	} 	
	T* value_cache = new T[get_num_data_per_element()];
	// loop over all elements
	for(Geometry::element_const_iterator eit = geom.elements_begin(); eit != geom.elements_end(); eit++) {				
		for(int ii = 0; ii < get_num_data_per_element(); ii++) {
			value_cache[ii] = this->get_element_value((*eit)->get_index_global(), ii);
		}
		// add values to nodes values and track how many values were added
		for(unsigned int vv = 0; vv < (*eit)->get_num_nodes(); vv++) {
			int node_index = (*eit)->get_node(vv).get_index_global();			
			track[node_index]++;
			int offset = node_index * get_num_data_per_element();			
			for(int ii = 0; ii < get_num_data_per_element(); ii++) {
				data[offset + ii] += value_cache[ii];
			}
		}
	}	
	// average nodes
	for(int ii = 0; ii < nnode; ii++) {
		int offset = ii * get_num_data_per_element();
		for(int jj = 0; jj < get_num_data_per_element(); jj++) {
			data[offset + jj] /= double(track[ii]);	
		}
	}
	delete[] value_cache;
	delete[] track;
	return data;	
}

/** XYData scalator (for unit change) */
template<class T>
class XYDataScale : public XYData<T> {
public:	
	XYDataScale(const XYData<T>& base_, const T& x_scale_, const T& y_scale_)
	: base(base_), x_scale(x_scale_), y_scale(y_scale_) {};
	virtual ~XYDataScale() {};	
	virtual int    get_x_length()                const;
	virtual int    get_num_y_sets()              const;
	virtual int    get_num_x_sets()              const; 
	virtual void   get_x(int xidx, vector<T> &x) const;
	virtual void   get_y(int yidx, vector<T>& y) const; 
	virtual string get_x_identifier(int xidx)    const;
	virtual string get_y_identifier(int yidx) 	 const;
private:
	const XYData<T>& base;
	T x_scale;
	T y_scale;
}; 

template<class T>
int XYDataScale<T>::get_x_length() const {
	return base.get_x_length();	
}
template<class T>
int XYDataScale<T>::get_num_y_sets() const {
	return base.get_num_y_sets();	
}
template<class T>
int XYDataScale<T>::get_num_x_sets() const {
	return base.get_num_x_sets();	
}
template<class T>
void XYDataScale<T>::get_x(int xidx, vector<T> &x) const {
	base.get_x(xidx, x);
	for(typename vector<T>::iterator it = x.begin(); it != x.end(); it++) {
		(*it) *= x_scale;	
	}	
}
template<class T>
void XYDataScale<T>::get_y(int yidx, vector<T>& y) const {
	base.get_y(yidx, y);
	for(typename vector<T>::iterator it = y.begin(); it != y.end(); it++) {
		(*it) *= y_scale;	
	}		
}
template<class T>
string XYDataScale<T>::get_x_identifier(int xidx) const {
	return base.get_x_identifier(xidx);
}
template<class T>
string XYDataScale<T>::get_y_identifier(int yidx) const {
	return base.get_y_identifier(yidx);	
}

/** merge 2 or more element data sets into one (for plotting purposes) */
template<class T>
class MergedElementData : ElementData<T> {
public:
	MergedElementData() {} ;
	virtual ~MergedElementData() {}
	virtual const T& get_element_value(unsigned int eidx, unsigned int eqidx = 0) const;
	virtual T&       get_element_value(unsigned int eidx, unsigned int eqidx = 0);
	virtual void     set_element_value(unsigned int eidx, unsigned int eqidx, T value);
	virtual int      get_length() const; 			
	virtual void     set_length(int length, int num_data_per_element);
	virtual int      get_num_data_per_element() const;
	virtual string   get_identifier(unsigned int eqidx) const; 
	void add_data(const ElementData<T>* set, const string& identifier);
	void add_data(const ElementData<T>* set, const char* identifier);
private:
	struct IndexSet {
		int idx_dataset;
		int idx_eqidx;
	};
	vector<const ElementData<T>*> datasets;
	vector<string>                labels;
	vector<IndexSet>              indexes;
};  

template<class T>
void MergedElementData<T>::add_data(const ElementData<T>* set, const char*identifier) {
	this->add_data(set, string(identifier));	
}

/** add element data */
template<class T>
void MergedElementData<T>::add_data(const ElementData<T>* set, const string& identifier) {

	// ----------------------------------------------
	// test if this one is compatible
	// ----------------------------------------------
	if(datasets.size() > 0) {
		TDKP_ASSERT(set->get_length() == datasets.front()->get_length(), "passed set has not same length");
	}
	
	// ----------------------------------------------
	// add it
	// ----------------------------------------------
	datasets.push_back(set);
	labels.push_back(identifier);
	for(int ii = 0; ii < set->get_num_data_per_element(); ii++) {
		IndexSet a;
		a.idx_dataset = datasets.size() - 1;	
		a.idx_eqidx   = ii;
		indexes.push_back(a);
	}
		
}

template<class T>
int MergedElementData<T>::get_num_data_per_element() const {
	return indexes.size(); 	
}

template<class T>
const T& MergedElementData<T>::get_element_value(unsigned int eidx, unsigned int eqidx) const {
	TDKP_BOUNDS_ASSERT(eqidx < indexes.size(), "eqidx < indexes.size()");
	const IndexSet& a = indexes[eqidx];
	return datasets[a.idx_dataset]->get_element_value(eidx, a.idx_eqidx); 	
}
template<class T>
T& MergedElementData<T>::get_element_value(unsigned int eidx, unsigned int eqidx) {
	return const_cast<T&>(static_cast<const MergedElementData<T>&>(*this).get_element_value(eidx,eqidx));
}
template<class T>
void     MergedElementData<T>::set_element_value(unsigned int eidx, unsigned int eqidx, T value) {
	TDKP_GENERAL_EXCEPTION("read only instance");	
}
template<class T>
int MergedElementData<T>::get_length() const {
 	TDKP_BOUNDS_ASSERT(datasets.size() > 0, "datasets.size()"); 
 	return datasets.front()->get_length();
}
template<class T> 			
void MergedElementData<T>::set_length(int length, int num_data_per_element) {
	TDKP_GENERAL_EXCEPTION("read only instance");
}

template<class T>
string MergedElementData<T>::get_identifier(unsigned int eqidx) const { 
	TDKP_BOUNDS_ASSERT(eqidx < indexes.size(), "eqidx < indexes.size()");
	const IndexSet& a = indexes[eqidx];
	ostringstream sout;
	sout << labels[a.idx_dataset] << datasets[a.idx_dataset]->get_identifier(a.idx_eqidx);
	return sout.str(); 		
}


} //e nd of namespace tdkp

#endif /*DATATYPES_H_*/
