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

#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "tdkp/povray/IsoSurfaceGenerator.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"

using namespace std;

namespace tdkp
{

IsoSurfaceGenerator::IsoSurfaceGenerator() {
	this->lists_initialized = false;
	this->geometry          = 0;
	this->result            = 0;
}

IsoSurfaceGenerator::~IsoSurfaceGenerator() {
	if(this->lists_initialized) {
		for(edge_iterator it = edges.begin(); it != edges.end(); it++) {
			delete (*it);	
		}				
	}
	if(this->result != 0) {
		delete[] this->result; 
	} 		
}

void IsoSurfaceGenerator::set_geometry(Geometry* geom) {
	this->geometry = geom;
}

void IsoSurfaceGenerator::set_data(double* res) throw(Exception*) {
	if(this->geometry == 0) {
		TDKP_GENERAL_EXCEPTION("geometry must be set before data");	
	}
	this->result = res;	
}

void IsoSurfaceGenerator::set_data(const NodeData<double>* res) throw(Exception*) {
	if(this->geometry == 0) {
		TDKP_GENERAL_EXCEPTION("geometry must be set before data");	
	}

	TDKP_ASSERT(res->get_num_data_per_node() == 1, "res->get_num_data_per_node() != 1");
	TDKP_ASSERT(res->get_length() == (signed)geometry->get_num_nodes(), "res->get_length() == geometry->get_num_nodes()");
	this->result = new double[res->get_length()];
	for(int ii = 0; ii < res->get_length(); ii++) {
		this->result[ii] = res->get_node_value(ii); 		
	}

}

/** reset for new isosurface */
void IsoSurfaceGenerator::clear() {
	
	for(unsigned int ii = 0; ii < edges.size(); ii++) {
		edges[ii]->has_value = false;
		edges[ii]->dp.nn = 0;
		edges[ii]->dp.all_normals.clear();
		edges[ii]->dp.normal = Vector3D(0.0,0.0,0.0);
		edges[ii]->dp.coords = Vector3D(0.0,0.0,0.0);
		edges[ii]->dp.id = 0;
	}
	structure_triangles.clear();
	faces.clear();
	node_povid2node_gid.clear();
}
	
void IsoSurfaceGenerator::create_lists() throw(Exception*) {
	if(this->geometry == 0 || this->result == 0) {
		TDKP_GENERAL_EXCEPTION("geometry and data must be set before creating lists");	
	}	
	if(this->lists_initialized) {
		TDKP_GENERAL_EXCEPTION("create_lists must not be called twice");	
	}
	// ----------------------------------------------------------
	// init node2edges / elements2edges
	// ----------------------------------------------------------
	this->node2edges.resize(this->geometry->get_num_nodes());
	this->elements2edges.resize(this->geometry->get_num_elements());
	
	Geometry::element_iterator eit;
	list<Edge*>::iterator edge_it;
	bool swapped;
	int  tmp;
	// --------------------------------------------------------------------
	// loop over elements and build edge list
	// edge (ii,jj) equals (jj,ii)
	// --------------------------------------------------------------------	
	int next  = 0;
	Logger::get_instance()->init_progress_bar("creating edge list", this->geometry->get_num_elements());
	int count = 0;
	int global_idx[4];
	for(eit = this->geometry->elements_begin(); eit != this->geometry->elements_end(); eit++) {
		if((*eit)->get_num_nodes() != 4) {
			TDKP_GENERAL_EXCEPTION("iso surface generator works only for tetraeders ...");				
		}
		// get global indices
		for(unsigned int ii = 0; ii < 4; ii++) {
			global_idx[ii] = (*eit)->get_node(ii).get_index_global();
		}
		// loop over nodes
		for(int ii = 0; ii < 4; ii++) {
			for(int jj = ii + 1; jj < 4; jj++) {
				// swap if necessary
				if(global_idx[ii] > global_idx[jj]) {
					swapped = true;	
					tmp     = ii;
					ii      = jj;
					jj      = tmp;
				} else {
					swapped = false;	
				}
				// look up if edge (ii,jj) is already defined
				bool defined = false;
				for(edge_it = (this->node2edges[global_idx[ii]]).begin(); edge_it != (this->node2edges[global_idx[ii]]).end(); edge_it++) {
					if((*edge_it)->v2->get_index_global() == global_idx[jj]) {
						defined = true;	
						break;
					}	
				}
				// edge not defined -> create new				
				if(!defined) {
					Edge* new_edge      = new Edge(); 
					new_edge->v1        = &(*eit)->get_node(ii);
					new_edge->v2        = &(*eit)->get_node(jj);
					new_edge->has_value = false;
					new_edge->dp.nn     = 0;
					this->edges.push_back(new_edge);
					this->node2edges[global_idx[ii]].push_back(new_edge);
					this->node2edges[global_idx[jj]].push_back(new_edge);
					this->elements2edges[count].push_back(new_edge);
				} else {
				// edge already defined, add to element2edges map
					this->elements2edges[count].push_back((*edge_it));					
				}
				if(swapped) {
					tmp     = ii;
					ii      = jj;
					jj      = tmp;						
				}
			}			
		}
				
		if(next == count++) {
			next = Logger::get_instance()->set_progress_bar(count,	this->geometry->get_num_elements());
		}
	}
	Logger::get_instance()->end_progress_bar();
	this->lists_initialized = true;
}

void IsoSurfaceGenerator::create_surface_plot() throw(Exception*) {
	if(!this->lists_initialized) {
		TDKP_GENERAL_EXCEPTION("edge lists must be initialized before call");		
	}
	double min = this->result[0]; 
	double max = this->result[0];		
	int    num = this->geometry->get_num_nodes();
	for(int ii = 0; ii < num; ii++) {
		if(this->result[ii] < min) {
			min = this->result[ii];	
		} else if(this->result[ii] > max) {
			max = this->result[ii];	
		}			
	}	
	double in = 0.0;
	while(in == 0) {
		std::cout << "data for surface plot is between " << min << " and " << max << ". which value do you want me to plot? ";
		cin >> in;
		std::cout << std::endl;
	}
	string filename;
	while(filename.size() == 0) {
		std::cout << "enter name of file to store povray file: ";
		cin >> filename;
		std::cout << std::endl;				 	
	}
	this->create_surface_plot(filename,in);
}

void IsoSurfaceGenerator::create_surface_plot(const char* filename, double value) throw(Exception*) {
	string tmp(filename);
	this->create_surface_plot(tmp, value);		
}

void IsoSurfaceGenerator::create_surface_plot(const string& filename, double value, bool neglect_structure) throw(Exception*) {
	if(!this->lists_initialized) {
		TDKP_GENERAL_EXCEPTION("edge lists must be initialized before call");		
	}
	if(filename.size() == 0) {
		TDKP_GENERAL_EXCEPTION("missing filename");	
	}	
	fstream fout(filename.c_str(), std::ios::out);
	if(fout.bad()) {
		TDKP_GENERAL_EXCEPTION("can not open file for writing");	
	}
	
	this->clear();
			
	// -----------------------------------------------------------
	// create datapoints on edges
	// -----------------------------------------------------------	
	int num_dp;	
	num_dp = this->create_datapoints(value); 	
	if(num_dp == 0) {
		TDKP_GENERAL_EXCEPTION("there is no iosurface for that value");	
	}
	ostringstream sout;
	sout << "there are " << this->edges.size() << " edges and " << num_dp << " of them cut the isosurface.";
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	this->create_faces(value);
	this->create_structure_triangles(); 
	
	// -----------------------------------------------------------
	// write datapoints to stream
	// -----------------------------------------------------------		
	this->write_header_to_stream(fout, value);
	fout << " mesh2 {\n"
	     << "   node_vectors {\n"
	     << "     " << num_dp << ",\n";		     
	this->write_datapoints_to_stream(fout, num_dp);
	fout << "   }\n"
		 << "	normal_vectors {\n"
		 << "     " << num_dp << ",\n";		     
	this->write_normals_to_stream(fout);	 
	fout << "   }\n"
		 << "   face_indices {\n"
         << "     " << this->faces.size() << ",\n";
	this->write_faces_to_stream(fout);
	fout << "   }\n" 
		 << "   pigment {rgb <0.9,0.8, 0.1>}\n"
		 << "}\n";
	if(!neglect_structure) {
		fout << "//  ------- structure triangles ------------\n";
		this->write_region_interfaces_to_stream(fout);
	}		 
	fout.close();
			
}

/** calculates the x in value = a * x + (1-x) * b
 * 
 * @return x in [0,1]. if v not between a and b, -1 is returned
 */
double IsoSurfaceGenerator::get_interpolation_factor(double a, double b, double value) const {
	if(((a < b) && (a < value && value < b)) || ((a > b) && (b < value && value < a))) {		
		return (value - b) / (a - b);
	} else if(a == value) {
		return 1.0;
	} else if(b == value) {
		return 0.0;	
	} else {
		return -1.0;
	}
}

/** loop over edges and create datapoints where edge bisects isosurface
 * 
 * stores also some information if a node is involved in a cut (for
 * faster lookup if an element has to be checked)
 * 
 * @param out output stream where we write the datapoints to
 * @param value isosurface value
 * @return number of datapoints created
 */
int IsoSurfaceGenerator::create_datapoints(double value) {
	
	edge_iterator edge_it;	

	double x;
	int num_dp_created = 0;
	
	// init vector 
	bool myfalse = false;	
	this->involved_node.assign(this->geometry->get_num_nodes(),myfalse);
	
	// loop over edges and create datapoints							
	int next = 0;
	int count = 0;
	Logger::get_instance()->init_progress_bar("analyzing edges and creating datapoints", this->edges.size());	
	for(edge_it = this->edges.begin(); edge_it != this->edges.end(); edge_it++) {
		x = this->get_interpolation_factor(
				this->result[(*edge_it)->v1->get_index_global()], 
				this->result[(*edge_it)->v2->get_index_global()],
				value
			);
		// -----------------------------------------------------------------------	
	 	// we do not handle the case where the value is exact a nodal value ... 
	 	// if we handle it, we also need to check create_faces() to scan for
	 	// such cases. 
	 	// -----------------------------------------------------------------------	
		TDKP_ASSERT((x >= 0 && x <= 1.0) || (x == -1.0), "(x >= 0 && x <= 1.0) || (x == -1.0) (x = " << x << ", a = " << this->result[(*edge_it)->v1->get_index_global()] << ", b = " << this->result[(*edge_it)->v2->get_index_global()] << ")");
		if(x == 0.0 || x == 1.0) {
			TDKP_GENERAL_EXCEPTION("the plot value is equal to a value at a given node ...");	
		}
		if(x != -1.0) {
/*			std::cout << "created datapoint: x:" << x 
					  << " v1: " << (*edge_it)->v1->get_index_global()
					  << " (val: " << this->result[(*edge_it)->v1->get_index_global()] << ") "
					  << " v2: " << (*edge_it)->v2->get_index_global()
					  << " (val: " << this->result[(*edge_it)->v2->get_index_global()] << ")\n";					  
*/						
			// create datapoints
			(*edge_it)->has_value = true;
			(*edge_it)->dp.coords = Vector3D(
				(*edge_it)->v1->get_coords(),
				(*edge_it)->v2->get_coords(),
				x
			);				
			(*edge_it)->dp.id = num_dp_created++;
			// mark node as involved ... 
			this->involved_node[(*edge_it)->v1->get_index_global()] = true;
			this->involved_node[(*edge_it)->v2->get_index_global()] = true;									
		}					
		if(next == count++) {
			next = Logger::get_instance()->set_progress_bar(count,	this->edges.size());
		}
	}
	Logger::get_instance()->end_progress_bar();
	return num_dp_created;
}

void IsoSurfaceGenerator::create_faces(double value) throw(Exception*){
	
	Geometry::element_iterator eit;
	list<Edge*>::iterator      edge_it;
	int  elem_idx = 0;
	int  next     = 0;	
	int  face_indices[6];
	int  check_edges = 0; 
	Face tmp_face;
	int  global_idx[4];

	Logger::get_instance()->init_progress_bar("scanning elements and creating intersecting triangles", this->geometry->get_num_elements());				
	for(eit = this->geometry->elements_begin(); eit != this->geometry->elements_end(); eit++) {				
		// get global indices
		bool involved = false;		
		for(unsigned int ii = 0; ii < 4; ii++) {
			global_idx[ii] = (*eit)->get_node(ii).get_index_global();
			if(this->involved_node[global_idx[ii]]) {
				involved = true;	
				break;
			}
		}
		// only proceed if at least a node sits on an edge cutting the isosurface
		if(involved) {
			// check all edges and count number edges cutting the isosurface
			int cuts = 0; 
			check_edges = 0;
			for(edge_it = this->elements2edges[elem_idx].begin(); edge_it != this->elements2edges[elem_idx].end(); edge_it++) {
				if((*edge_it)->has_value) {
					face_indices[cuts] = (*edge_it)->dp.id;
					cuts++;	
				}	
				check_edges++;
			}
			TDKP_ASSERT(check_edges == 6, "check_edges == 6");
			// if three edges are cut, we have a triangle
			if(cuts == 3) {
				for(int ii = 0; ii < 3; ii++) {
					tmp_face.id[ii] = face_indices[ii];	
				}							
				// add to faces
				add_normal_to_dp(tmp_face, this->elements2edges[elem_idx], value);
				this->faces.push_back(tmp_face);		
			} else if (cuts == 4) {
			// cut is a rectangle -> we split it into two triangles
				for(int ii = 0; ii < 3; ii++) {
					tmp_face.id[ii] = face_indices[ii];	
				}
				// add to faces
				add_normal_to_dp(tmp_face, this->elements2edges[elem_idx], value);
				this->faces.push_back(tmp_face);		
				for(int ii = 1; ii < 4; ii++) {
					tmp_face.id[ii - 1] = face_indices[ii];	
				}
				// add to faces
				add_normal_to_dp(tmp_face, this->elements2edges[elem_idx], value);
				this->faces.push_back(tmp_face);						
			} else if (cuts > 4) {
				// very unlikely case of isosurface beeing exactly on a node ... so we don't
				// handle this case now, but just quit
				std::cout << "================== E R R O R =====================\n";
				for(edge_it = this->elements2edges[elem_idx].begin(); edge_it != this->elements2edges[elem_idx].end(); edge_it++) {
					if((*edge_it)->has_value) {
						std::cout << "edge (" << (*edge_it)->v1->get_index_global() << ","	<< (*edge_it)->v2->get_index_global() << ") "
								  << " res: " << this->result[(*edge_it)->v1->get_index_global()] << ", " 
								  << this->result[(*edge_it)->v2->get_index_global()] << "\n";
					}	
				}												
				TDKP_GENERAL_EXCEPTION("value is equal result at node in the element. case not handled yet. sorry.");
			} // else <- cut does not go through tetraeder
			else if(cuts > 0){
				std::cout << "cuts only " << cuts << "\n";	
			}
		}
		if(next == elem_idx++) {
			next = Logger::get_instance()->set_progress_bar(elem_idx,	this->geometry->get_num_elements());
		}
	}
	Logger::get_instance()->end_progress_bar();
	
}

void IsoSurfaceGenerator::add_normal_to_dp(Face &face, list<Edge*> &edge_list, double value) {
	
    // 1. create plane vectors
    Vector3D inplane_1, inplane_2;
    Vector3D point_coordinates[3];      
    Vector3D to_smaller[3];
    double   area;
    
    // loop over edges to get data of our data points
    for(list<Edge*>::const_iterator it = edge_list.begin(); it != edge_list.end(); it++) {
        for(int ii = 0; ii < 3; ii++) {
            // if face data point is equal to edge data pint
            if((*it)->has_value && (*it)->dp.id == face.id[ii]) {
                if(this->result[(*it)->v1->get_index_global()] < value) {
                    to_smaller[ii] = Vector3D((*it)->v1->get_coords()) - (*it)->dp.coords;      
                    to_smaller[ii].normalize();                 
                } else if(this->result[(*it)->v2->get_index_global()] < value) {
                    to_smaller[ii] = Vector3D((*it)->v2->get_coords()) - (*it)->dp.coords;      
                    to_smaller[ii].normalize();                 
                } else {
                    std::cout << "problem at edge here:\n";
                    std::cout << "dp: " << (*it)->dp.id << "\n";
                    for(int jj = 0; jj < 3; jj++) {
                        std::cout << face.id[jj] << "  ";
                    }
                    cout << "\n";
                                    for(list<Edge*>::const_iterator lit = edge_list.begin(); lit != edge_list.end(); lit++) {
                                            std::cout << "dpid: " << (*lit)->dp.id << "  ("
                                    << (*lit)->v1->get_index_global() << ","
                            << (*lit)->v2->get_index_global() << ") "
                            << "values: " << this->result[(*lit)->v1->get_index_global()] << "  "
                            << this->result[(*lit)->v2->get_index_global()] << ".\n";
                    }                                  
                    TDKP_GENERAL_EXCEPTION("should not be here");   
                }   
                TDKP_ASSERT(point_coordinates[ii].norm() == 0.0, "point_coordinates[ii].norm() == 0.0");
                // store for later 
                point_coordinates[ii] = (*it)->dp.coords;                       
                break;
            }
        }
    } 
    

    // ----------------------------------------------------------------             
    // create inplane vector 
    // ----------------------------------------------------------------
    inplane_1 = point_coordinates[1] - point_coordinates[0];
    inplane_2 = point_coordinates[2] - point_coordinates[0];   

	// 2. take cross prod (for normal)
	Vector3D normal = Vector3D::cross_product(inplane_1, inplane_2);
	area = 0.5 * normal.norm();
	normal.normalize();
	inplane_1.normalize();
	inplane_2.normalize();
	
	double leabs;
	if( (leabs = tdkp_math::abs(Vector3D::dot_product(inplane_1,normal))) > 1.0e-9) {
		std::cout << "tdkp_math::abs(Vector3D::dot_product(inplane_1,normal)) > 1.0e-9\n";
	} 
	if( (leabs = tdkp_math::abs(Vector3D::dot_product(inplane_2,normal))) > 1.0e-9) {
		std::cout << "tdkp_math::abs(Vector3D::dot_product(inplane_2,normal)) > 1.0e-9\n";
	}
		
	// 3. swap sign of normal if necessary
//	if(Vector3D::dot_product(normal, to_smaller) < 0.0) {		
//		normal = normal * -1.0;	
//	}	
	
	
	// 4. add normal to dot normals (but do not normalize it (read expl. in write_normals_to_stream)
	// loop over edges to get data of our data points
	for(list<Edge*>::const_iterator it = edge_list.begin(); it != edge_list.end(); it++) {
		for(int ii = 0; ii < 3; ii++) {
			// if face data point is equal to edge data pint
			if((*it)->dp.id == face.id[ii]) {
				if(Vector3D::dot_product(normal, to_smaller[ii]) < 0.0) {		
					normal = normal * -1.0;	
				}	/*										
				if(Vector3D::dot_product((*it)->dp.normal,normal) < 0.0) {
					std::cout << "WARNING: normal at " << (*it)->dp.id << " " << (*it)->dp.normal << " gets contribution in " << normal << "\n";	
				}		*/		
				(*it)->dp.normal = (*it)->dp.normal + normal * area;
				(*it)->dp.nn++;
				(*it)->dp.all_normals.push_back(normal);
			}
		}
	}		
}

void IsoSurfaceGenerator::write_header_to_stream(ostream &fout, double value) const {
	if(this->geometry != 0) {
		// ----------------------------------------------------------------
		// calculate bounding box
		// ----------------------------------------------------------------
		edge_const_iterator it;
		Vector3D min_coords;
		Vector3D max_coords;
		Vector3D min_structure, max_structure;
		Vector3D middle;
		Vector3D camera;
		Vector3D light;
		Vector3D coords;
		
		bool first = true;
		for(it = this->edges.begin(); it != this->edges.end(); it++) {
			if((*it)->has_value) {
				coords = (*it)->dp.coords;				
				if(first) {		
					min_coords = max_coords = coords;	
					first = false;
				} else {
					for(int ii = 0; ii < 3; ii++) {				
						if(coords(ii) < min_coords(ii)) {
							min_coords(ii) = coords(ii);	
						} else if(coords(ii) > max_coords(ii)) {
							max_coords(ii) = coords(ii);	
						}
					}
				}
			}
		}
		map<int,int>::const_iterator sit = this->node_povid2node_gid.begin();
		min_structure = max_structure = Vector3D(this->geometry->get_node((*sit).second).get_coords());			
		for(; sit != this->node_povid2node_gid.end(); sit++) {
			coords = Vector3D(this->geometry->get_node((*sit).second).get_coords());									
			for(int ii = 0; ii < 3; ii++) {				
				if(coords(ii) < min_structure(ii)) {
					min_structure(ii) = coords(ii);	
				} else if(coords(ii) > max_structure(ii)) {
					max_structure(ii) = coords(ii);	
				}
			}
		}
		
		// calculate middle of bounding box (take (min + max) / 2)
		// take light on line min_coord - max_coord
		middle = (min_coords + max_coords) * 0.5;	
	    light  = 3.5 *  max_structure - 2.5 * min_structure;
		
		// take camera at x axis 2.5 x distance between max x and middle
		camera = middle;
		camera(0) = 4.5 * (max_coords(0) - middle(0));		
						
		// ----------------------------------------------------------------
		// write header 
		// ----------------------------------------------------------------
		fout.precision(12);
		fout << "// -------------------------------------------------------\n"
		     << "//               povray file created by tdkp\n\n"
		     << "// plot at value: " << value << "\n"
		     <<	"// \n"
		     << "// grid: \n"
		     << "//   elements: " << this->geometry->get_num_elements() << "\n"
		     << "//   nodes:    " << this->geometry->get_num_nodes() << "\n"
		     << "// bounding box: \n"
		     << "//   min: (" << min_coords(0) << ", " << min_coords(1) << ", " << min_coords(2) << ")\n"
		     << "//   max: (" << max_coords(0) << ", " << max_coords(1) << ", " << max_coords(2) << ")\n"
		     << "// structure box:\n"
		     << "//   min: " << min_structure << "\n"
		     << "//   max: " << max_structure << "\n"
		     << "// -------------------------------------------------------\n"
			 << "\n" 
			 << "// background { color rgb<1, 1, 1> } \n"
			 << "//camera {\n"
			 << "//   location " << camera << "\n"
			 << "//   look_at <" << middle(0) << ", " << middle(0) << ", " << middle(2) << ">\n"
			 << "//   sky <0,0,1>\n"   
			 << "//}\n"
			 << "//light_source { <" << light(0) << ", " << light(1) << ", " << light(2) << "> color rgb<1, 1, 1> }\n"
			 << "//light_source { <" << -light(0) << ", " << light(1) << ", " << light(2) << "> color rgb<1, 1, 1> }\n"
			 << " //light_source { <" << light(0) << ", " << -light(1) << ", " << light(2) << "> color rgb<1, 1, 1> }\n"
			 << "// light_source { <" << light(0) << ", " << light(1) << ", " << -light(2) << "> color rgb<1, 1, 1> }\n"
			 << "// light_source { <" << -light(0) << ", " << -light(1) << ", " << light(2) << "> color rgb<1, 1, 1> }\n"			 			 			 			 			 
			 << "\n";		     		     
	} 
}

void IsoSurfaceGenerator::write_datapoints_to_stream(ostream &fout, int num_dp) const {
	string comma = "";
	int count = 0;
	for(edge_const_iterator edge_it = this->edges.begin(); edge_it != this->edges.end(); edge_it++) {		
		if((*edge_it)->has_value) {
			fout << comma << (*edge_it)->dp.coords;  
			comma = ",\n    ";		         
			TDKP_ASSERT(count == (*edge_it)->dp.id, "count == (*edge_it)->dp.id");
			count++;
		}		
	}
	fout << "\n";
	TDKP_ASSERT(count == num_dp, "count == num_dp");		
}
	
void IsoSurfaceGenerator::write_faces_to_stream(ostream &fout) const {
	string comma = "";
	for(vector<Face>::const_iterator it = this->faces.begin(); it != this->faces.end(); it++) {
		fout << comma << "<" << (*it).id[0] << "," << (*it).id[1] << "," << (*it).id[2] << ">";
		comma = ",\n    ";		         	
	}		
}

void IsoSurfaceGenerator::write_normals_to_stream(ostream& fout) {
	string comma = "";
	for(edge_const_iterator edge_it = this->edges.begin(); edge_it != this->edges.end(); edge_it++) {		
		if((*edge_it)->has_value) {
			(*edge_it)->dp.normal.normalize();
			fout << comma << (*edge_it)->dp.normal << " // (" << (*edge_it)->dp.id << ")  "; 
/*			for(int ii = 0; ii < (*edge_it)->dp.all_normals.size(); ii++) {
				fout << (*edge_it)->dp.all_normals[ii] << "  ";	
			}*/
			comma = ",\n    ";		         
		}		
	}
	fout << "\n";	
}

/** write region interfaces to stream 
 */
void IsoSurfaceGenerator::create_structure_triangles() throw(Exception*) {
	
	Geometry::boundary_const_iterator it;
	int pov_idx = 0;
	vector<int> nodes;
	map<int,int> node_gid2node_povid;
	// loop over all elements
	for(it = this->geometry->element_boundaries_begin(); 
	    it != this->geometry->element_boundaries_end(); it++) {
		if((*it)->get_location() == 'f') {
			nodes.resize((*it)->get_num_vertices());
			for(unsigned int ii = 0; ii < (*it)->get_num_vertices(); ii++) {
				nodes[ii] = (*it)->get_vertex(ii).get_index_global();	
			}			
			// add missing nodes to global map
			for(int ii = 0; ii < 3; ii++) {
				map<int,int>::iterator gid_it = node_gid2node_povid.find(nodes[ii]);
				if(gid_it == node_gid2node_povid.end()) {
					this->node_povid2node_gid[pov_idx] = nodes[ii];
					node_gid2node_povid[nodes[ii]] = pov_idx++;	
				}
			}		
			Face tmp_face;	
			for(int ii = 0; ii < 3; ii++) {
				tmp_face.id[ii] = node_gid2node_povid[nodes[ii]];
			}
			this->structure_triangles.push_back(tmp_face);	
		}
	}	
}

void IsoSurfaceGenerator::write_region_interfaces_to_stream(ostream &fout) const {
	string comma = "";	
	fout << " mesh2 {\n"
	     << "   node_vectors {\n"
	     << "     " << this->node_povid2node_gid.size() << ",\n";		     
	for(map<int,int>::const_iterator it = this->node_povid2node_gid.begin(); it != this->node_povid2node_gid.end(); it++) {
		fout << comma << Vector3D(this->geometry->get_node((*it).second).get_coords());	
		comma = ",\n    ";		         	
	} 	
	fout << "   }\n"
		 << "   face_indices {\n"
	     << "     " << this->structure_triangles.size() << ",\n";		     
	comma = " ";
	for(vector<Face>::const_iterator it = this->structure_triangles.begin(); it != this->structure_triangles.end(); it++) {
		fout << comma << "<" << (*it).id[0] << "," << (*it).id[1] << "," << (*it).id[2] << ">";
		comma = ",\n    ";		         	
	}				
	fout << "   }\n" 
		 << "   pigment {rgbf <0.9,0.8, 0.1, 0.9>}\n"
		 << "   finish { phong 0.8 ambient 0.05 }\n"
   		 << "}\n";		
}
   											
} // end namespace tdkp
