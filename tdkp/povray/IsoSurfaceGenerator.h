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

#ifndef ISOSURFACEGENERATOR_H_
#define ISOSURFACEGENERATOR_H_

#include <string>
#include <list>
#include <vector>
#include <map>

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/geometry/Node.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/common/DataTypes.h"

using namespace std;

namespace tdkp
{

class IsoSurfaceGenerator
{
public:
	IsoSurfaceGenerator();
	void set_geometry(Geometry* geom);
	void set_data(double* res) throw(Exception*);
	void set_data(const NodeData<double>* res) throw(Exception*);	
	void create_lists() throw(Exception*);
	void create_surface_plot() throw(Exception*);
	void create_surface_plot(const char* filename,   double value) throw(Exception*);
	void create_surface_plot(const string& filename, double value, bool neglect_structure = true) throw(Exception*);
	void clear();
	
	virtual ~IsoSurfaceGenerator();
	
private:
	struct DataPoint {
		Vector3D coords;
		Vector3D normal;
		int nn;
		vector<Vector3D> all_normals;
		int id;			
	};
	struct Edge {
		Node *v1,*v2;								
		bool has_value;
		DataPoint dp;
	};
	struct Face {
		int id[3];	
	};
			
	bool lists_initialized;
	// data plot
	vector<Edge*> edges;
	vector<list<Edge*> > node2edges;
	vector<list<Edge*> > elements2edges;
	vector<bool>         involved_node;
	vector<Face>         faces;
	
	// structure interface plot
	vector<Face>         structure_triangles;
	map<int,int>         node_povid2node_gid;
										
	typedef vector<Edge*>::iterator       edge_iterator;								
	typedef vector<Edge*>::const_iterator edge_const_iterator;	
	Geometry* geometry;
	double*   result;
	
	double get_interpolation_factor(double a, double b, double value) const;
	int    create_datapoints(double value);
	void   create_faces(double value) throw(Exception*);
    void   create_structure_triangles() throw(Exception*);
	void   add_normal_to_dp(Face &face, list<Edge*> &edge_list, double value);
	void   write_header_to_stream(ostream &fout, double value) const;
	void   write_datapoints_to_stream(ostream& fout, int num_dp) const;
	void   write_normals_to_stream(ostream& fout);
	void   write_faces_to_stream(ostream &fout) const;	
	void   write_region_interfaces_to_stream(ostream &fout) const;
		
};





}

#endif /*ISOSURFACEGENERATOR_H_*/
