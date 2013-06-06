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

#include "tdkp/povray/DF3Generator.h"
#include "tdkp/geometry/Element3DTR.h"

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))

namespace tdkp
{


DF3Generator::DF3Generator() :
 xmin(0), ymin(0), zmin(0), xmax(0), ymax(0), zmax(0), nx(0), ny(0), nz(0) {
 	this->map      = 0;
	this->geometry = 0;
	this->data     = 0;	 		
}

/** deletes map and data, but does not touch geometry object */
DF3Generator::~DF3Generator() {
	if(this->map != 0) {
		delete map;	
	}
	if(this->data != 0) {
		delete data;	
	}
}

/** set node based data 
 * 
 * a grid needs to be defined frist
 */
void DF3Generator::set_data(NodeData<double>* res) throw(Exception*) {
	if(this->map == 0) {
		TDKP_GENERAL_EXCEPTION("define grid before setting data");	
	}	
	if(this->data != 0) {
		delete this->data;	
	}
	//df
		
	this->data = res;
	// calculate values
	for(int xx = 0; xx < nx; xx++) {
		for(int yy = 0; yy < ny; yy++) {
			for(int zz = 0; zz < nz; zz++) {
				NodeMap& vmap = this->map->get(xx,yy,zz);
				vmap.value = 0;
				for(int ii = 0; ii < 4; ii++) {
					vmap.value += vmap.contribution[ii] * this->data->get_node_value(vmap.nodes[ii]);	
				}	
			}	
		}	
	}		
}

void DF3Generator::set_df3_grid(short    nx_, short ny_, short nz_, 
                  				double xmin_, double ymin_, double zmin_, 
			                    double xmax_, double ymax_, double zmax_) {

	TDKP_ASSERT(nx_ > 0 && ny_ > 0 && nz_ > 0, "nx > 0 && ny > 0 && nz > 0");
	TDKP_ASSERT(xmin_ < xmax_ && ymin_ < ymax_ && zmin_ < zmax_, "xmin < xmax && ymin < ymax && zmin < zmax");
	TDKP_ASSERT(this->geometry != 0, "geometry object needs to be set first");
	
	if(this->map != 0) {
		delete map;	
	}		
	nx = nx_; ny = ny_; nz = nz_; 
	xmin = xmin_; ymin = ymin_; zmin = zmin_;
	xmax = xmax_; ymax = ymax_; zmax = zmax_;
	dx = (xmax - xmin) / double(nx - 1);
	dy = (ymax - ymin) / double(ny - 1);	
	dz = (zmax - zmin) / double(nz - 1);		
	this->map = new GridMap(nx, ny, nz);

	int num  = nx * ny * nz;
	int next = 0;
	int count = 0;
	Logger::get_instance()->init_progress_bar("building grid map", num);
	// for all nodes	
	for(short xx = 0; xx < nx; xx++) {
		for(short yy = 0; yy < ny; yy++) {
			for(short zz = 0; zz < nz; zz++) {
				
				bool notfound = true;
				double contrib[4];
				double coords[] = {double(xx) * dx + xmin, double(yy) * dy + ymin, double(zz) * dz + zmin}; 
				// find element where node is in
				for(Geometry::element_const_iterator it = this->geometry->elements_begin(); it != this->geometry->elements_end(); it++) {
					// store contributions
					if(dynamic_cast<Element3DTR*>((*it))->get_contribution(coords, contrib)) {
						NodeMap& vmap = this->map->get(xx,yy,zz); 
						double check = 0.0;
						for(int aa = 0; aa < 4; aa++) {
							vmap.contribution[aa] = contrib[aa];
							check += contrib[aa];
							vmap.nodes[aa]     = (*it)->get_node(aa).get_index_global();
						}
						if(fabs(check - 1.0) > 1.0e-9) {						
							cout << "p: (" << coords[0] << ", " << coords[1]  << ", " << coords[2] << ") "
							     << " matches element " << (*it)->get_index_global() << " but has contrib > 1: " << check << "\n";
						} 
						notfound = false;
						break;
					}
				}	
				if(notfound) {
					ostringstream sout;
					sout << "there is a problem: df3 node	(" << coords[0] << ", " << coords[1] 
					     << ", " << coords[2] << ") is not inside an element\n";
					TDKP_GENERAL_EXCEPTION(sout.str());					     
				}
				if(count == next) {
					next = Logger::get_instance()->set_progress_bar(count, num);
				}	
				count++;			
			}
		}	
	}
	Logger::get_instance()->end_progress_bar();
	
	
	
}
void DF3Generator::write_df3_file(const char* filename) throw(Exception*) {

	TDKP_ASSERT(sizeof(short) == 2, "sizeof(short) == 2");
	
	double themin, themax;
	double val;
	// determine bounds
	themin = themax = this->map->get(0,0,0).value;
	for(int xx = 0; xx < nx; xx++) {
		for(int yy = 0; yy < ny; yy++) {
			for(int zz = 0; zz < nz; zz++) {
				val = this->map->get(xx,yy,zz).value;
				themin = min(themin, val);
				themax = max(themax, val); 		
			}
		}
	}	
	cout << "values bounded by: " << themin << " and " << themax << endl;

			     
   	// open file
   	fstream fout(filename, ios::out | ios::binary);
   	if(!fout) {
   		TDKP_GENERAL_EXCEPTION("could not open file for writing");
   	}
   	// write in big endian
   	// 3 shorts with num x, num y, num z
   	fout.put(nx >> 8);
   	fout.put(nx & 0xff);
   	fout.put(ny >> 8);
   	fout.put(ny & 0xff);
   	fout.put(nz >> 8);
   	fout.put(nz & 0xff);
   	// and then data
   	unsigned short sdata;
   	for(int zz = 0; zz < nz; zz++) {
 		for(int yy = 0; yy < ny; yy++) {
 			for(int xx = 0; xx < nx; xx++) {
 				if(!((this->map->get(xx,yy,zz).value - themin) / (themax - themin) >= 0.0 && (this->map->get(xx,yy,zz).value - themin) / (themax - themin) <= 1.0)) {
					cout << "value out of bounds for " << xx << ", " << yy << ", " << zz << ": " 
					     << (this->map->get(xx,yy,zz).value - themin) / (themax - themin) << endl;
 				}
	            sdata = ushort(256 * 256 * (this->map->get(xx,yy,zz).value - themin) / (themax - themin));
	            FIX_SHORT(sdata); // transfers from little endian to big endian
	            fout.write((char*)&sdata, 2);	            
 			}	            
     	}     
	}
	fout.close();
}

void DF3Generator::write_map(const char* filename) {

	fstream fout(filename, ios::out | ios::binary);
	if(!fout) {
		TDKP_GENERAL_EXCEPTION("could not open file for writing");	
	}	
	
	fout.write((char*)(&nx), sizeof(int));
	fout.write((char*)(&ny), sizeof(int));
	fout.write((char*)(&nz), sizeof(int));		
	fout.write((char*)(&xmin), sizeof(double));
	fout.write((char*)(&ymin), sizeof(double));
	fout.write((char*)(&zmin), sizeof(double));	
	fout.write((char*)(&xmax), sizeof(double));		
	fout.write((char*)(&ymax), sizeof(double));	
	fout.write((char*)(&zmax), sizeof(double));					
	fout.write((char*)(&dx), sizeof(double));		
	fout.write((char*)(&dy), sizeof(double));	
	fout.write((char*)(&dz), sizeof(double));						

	for(int xx = 0; xx < nx; xx++) {
		for(int yy = 0; yy < ny; yy++) {
			for(int zz = 0; zz < nz; zz++) {
				NodeMap& vmap = this->map->get(xx,yy,zz);
				fout.write((char*)(vmap.contribution), sizeof(double) * 4);
				fout.write((char*)(vmap.nodes),     sizeof(int) * 4);				
			}
		}
	}

	fout.close(); 
	
}
void DF3Generator::read_map(const char* filename) {

	fstream fout(filename, ios::in | ios::binary);
	if(!fout) {
		TDKP_GENERAL_EXCEPTION("could not open file for writing");		
	}	
	if(this->map != 0) {
		TDKP_GENERAL_EXCEPTION("can only read maps in freshly created df3 generators");	
	}
		
	fout.read((char*)(&nx), sizeof(int));
	fout.read((char*)(&ny), sizeof(int));
	fout.read((char*)(&nz), sizeof(int));		
	fout.read((char*)(&xmin), sizeof(double));
	fout.read((char*)(&ymin), sizeof(double));
	fout.read((char*)(&zmin), sizeof(double));	
	fout.read((char*)(&xmax), sizeof(double));		
	fout.read((char*)(&ymax), sizeof(double));	
	fout.read((char*)(&zmax), sizeof(double));					
	fout.read((char*)(&dx), sizeof(double));		
	fout.read((char*)(&dy), sizeof(double));	
	fout.read((char*)(&dz), sizeof(double));						

	this->map = new GridMap(nx,ny,nz);

	for(int xx = 0; xx < nx; xx++) {
		for(int yy = 0; yy < ny; yy++) {
			for(int zz = 0; zz < nz; zz++) {
				NodeMap& vmap = this->map->get(xx,yy,zz);
				fout.read((char*)(vmap.contribution), sizeof(double) * 4);
				fout.read((char*)(vmap.nodes),     sizeof(int) * 4);				
			}
		}
	}

	fout.close(); 
		
}


}
