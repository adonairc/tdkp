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

#ifndef DF3GENERATOR_H_
#define DF3GENERATOR_H_

#include "tdkp/common/all.h"
#include "tdkp/geometry/Geometry.h"
#include "tdkp/main/EigenSolution.h"

namespace tdkp
{

class DF3Generator
{
public:
	DF3Generator();
	virtual ~DF3Generator();
	void set_geometry(Geometry* geom) { this->geometry = geom; TDKP_ASSERT(geom->get_dimension() == 3, "geometry dimension == 3"); }
	void set_data(NodeData<double>* res) throw(Exception*);
	void set_df3_grid(short    nx, short ny, short nz, 
	                  double xmin, double ymin, double zmin, 
	                  double xmax, double ymax, double zmax);	
	                  
	void write_map(const char* filename); 
	void read_map(const char* filename);	                  	
	void write_df3_file(const char* filename) throw(Exception*);
	
private:

	class NodeMap {
	public:
		NodeMap() { 
			for(short ii = 0; ii < 4; ii++) {
				nodes[ii]     = 0;
				contribution[ii] = 0.0;	
			}	
			value = 0;
		}
		int    nodes[4];
		double contribution[4];
		double value;
	};

	class GridMap {
	public:
		GridMap(short nx_, short ny_, short nz_) : nx(nx_), ny(ny_), nz(nz_) {
			map.resize(int(nx) * int(ny) * int(nz));
		}
		NodeMap& get(int x, int y, int z) {
			return map[get_idx(x,y,z)];
		}
		NodeMap& operator()(int x, int y, int z) {
			return map[get_idx(x,y,z)];				
		}
	private:
		short nx, ny, nz;
		vector<NodeMap> map;		
		int get_idx(int x, int y, int z) {
			return z + y * nz + x * nz * ny;	
		}
		
	};
	
	double xmin, ymin, zmin, xmax, ymax, zmax;
	double dx, dy, dz;
	int    nx, ny, nz;	
	
	GridMap*    map;
	Geometry*   geometry;
	NodeData<double>* data;	
						
};

}

#endif /*DF3GENERATOR_H_*/
