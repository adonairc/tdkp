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

#include "tdkp/probdefs/IntrinsicStrain.h"
#include "tdkp/io/InputParser.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/Configuration.h"

namespace tdkp {

void IntrinsicStrain::prepare() {	

	if(this->sparsity_copy == 0 || this->reference_lattice_constant <= 0.0) {
		TDKP_GENERAL_EXCEPTION("class is not ready to be prepared(?) - misses reference lattice constant or sparsity pattern");	
	} 
	
	ostringstream sout;
	sout << "IntrinsicStrain: the system is strained (or unstrained) into the following directions:\n";
	for(int ii = 0; ii < 3; ii++) {
		sout << "direction " << ii << ": "
		     << (strained_axes[ii] ? "strained" : "free")
		     << "\n"; 		 	
	}
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());
	
	this->properties.resize(0);
	// ----------------------------------------------
	// prepare initial intrinsic strains
	// ----------------------------------------------
	this->prepare_initial_intrinsic_strains();
	
	sout.str("");
	sout << "IntrinsicStrain: used elastic coefficents are:\n";
	sout << setw(15) << "Material" 
	     << setw(10) << "C11"
	     << setw(10) << "C12"
	     << setw(10) << "C44";


	region_properties prop;		
	double C11, C12, C44;
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {		
		// init coefficient matrices to zero
		for(int cc = 0; cc < 3; cc++) {			
			for(int dd = 0; dd < 3; dd++) {
				for(int ff = 0; ff < 3; ff++) {
					for(int gg = 0; gg < 3; gg++) {						
						prop.coefficent_matrices[cc][dd][ff][gg] = 0.0;
					}
				}
			}			
		}			
		
		C11 = this->material_db.get_material(ii)->get("elastic_constant_C11");
		C12 = this->material_db.get_material(ii)->get("elastic_constant_C12");
		C44 = this->material_db.get_material(ii)->get("elastic_constant_C44");		

		// -----------------------------------
		// write constants to screen
		// ----------------------------------- 
		sout << "\n" 
		     << setw(15) << this->material_db.get_material_name(ii)
		     << setw(10) << C11
		     << setw(10) << C12
		     << setw(10) << C44; 
 		  			
		// D_DX D_DX
		prop.coefficent_matrices[D_DX][D_DX][0][0] = C11;			
		prop.coefficent_matrices[D_DX][D_DX][1][1] = C44;
		prop.coefficent_matrices[D_DX][D_DX][2][2] = C44;
		// D_DY D_DY
		prop.coefficent_matrices[D_DY][D_DY][0][0] = C44;			
		prop.coefficent_matrices[D_DY][D_DY][1][1] = C11;
		prop.coefficent_matrices[D_DY][D_DY][2][2] = C44;			 
		// D_DZ D_DZ
		prop.coefficent_matrices[D_DZ][D_DZ][0][0] = C44;			
		prop.coefficent_matrices[D_DZ][D_DZ][1][1] = C44;
		prop.coefficent_matrices[D_DZ][D_DZ][2][2] = C11;			 			
		// D_DX D_DY
		prop.coefficent_matrices[D_DX][D_DY][0][1] = C12;			
		prop.coefficent_matrices[D_DX][D_DY][1][0] = C44;						
		// D_DY D_DX
		prop.coefficent_matrices[D_DY][D_DX][0][1] = C44;
		prop.coefficent_matrices[D_DY][D_DX][1][0] = C12;
		// D_DX D_DZ
		prop.coefficent_matrices[D_DX][D_DZ][0][2] = C12;
		prop.coefficent_matrices[D_DX][D_DZ][2][0] = C44;
		// D_DZ D_DX
		prop.coefficent_matrices[D_DZ][D_DX][0][2] = C44;
		prop.coefficent_matrices[D_DZ][D_DX][2][0] = C12;			
		// D_DY D_DZ
		prop.coefficent_matrices[D_DY][D_DZ][1][2] = C12;
		prop.coefficent_matrices[D_DY][D_DZ][2][1] = C44;
		// D_DZ D_DY
		prop.coefficent_matrices[D_DZ][D_DY][1][2] = C44;
		prop.coefficent_matrices[D_DZ][D_DY][2][1] = C12;			
				
		if(Configuration::get_instance()->get("output_strain_coefficent_matrices") == 1.0) {
			ostringstream mout;
			mout << "IntrinsicStrain: elastic coefficent matrices in material " 
			     << this->material_db.get_material_name(ii) << ":\n"; 								
			for(int uu = 0; uu < 3; uu++) {
				for(int jj = 0; jj < 3; jj++) {
					for(int kk = 0; kk < 3; kk++) {
						for(int ll = 0; ll < 3; ll++) {
							mout << "[" << uu << "][" << jj << "][" << kk 
								 << "][" << ll << "] = "
								 <<	prop.coefficent_matrices[uu][jj][kk][ll]
								 << "\n";				
						}	
					}	
				}
			}
			Logger::get_instance()->emit(LOG_INFO, mout.str());		
		}					
											
		// build strains
		prop.intrinsic_strain = this->initial_intrinsic_strains[ii];
		// build stresses loads
		for(int jj = 0; jj < 3; jj++) {
			int next = (jj + 1) % 3;
			int prev = (jj + 2) % 3;		
			prop.intrinsic_stress_load(jj, jj) = 
				prop.coefficent_matrices[jj][jj][jj][jj]     
			  * prop.intrinsic_strain(jj,jj) 
			  + prop.coefficent_matrices[jj][prev][jj][prev] 
			  * prop.intrinsic_strain(prev,prev) 
			  + prop.coefficent_matrices[jj][next][jj][next]
			  * prop.intrinsic_strain(next,next) 
			  + prop.coefficent_matrices[prev][prev][jj][jj]
			  * prop.intrinsic_strain(prev,jj)   
			  + prop.coefficent_matrices[next][next][jj][jj]
			  * prop.intrinsic_strain(next,jj);
			// shear / offdiagonal				       
			prop.intrinsic_stress_load(jj,prev) = 
				prop.coefficent_matrices[prev][prev][jj][jj] 
			  * prop.intrinsic_strain(jj,prev);
			prop.intrinsic_stress_load(jj,next) = 
				prop.coefficent_matrices[next][next][jj][jj]
			  * prop.intrinsic_strain(jj,next);			

		} 
		// store for later						
		this->properties.push_back(prop);
	}	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
}

/** calculate initial intrinsic strains
 * 
 * here, as we have with std. zincblende isotropic
 * material, we just create a diagonal strain tensor
 */ 
void IntrinsicStrain::prepare_initial_intrinsic_strains() {
	this->initial_intrinsic_strains.resize(0);	
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		// disregard materials without lattice constant (gas)		
		if(this->material_db.get_material(ii)->valid_key("lattice_constant")) {									
			// see comment on set_reference_lattice_constant			
			const double lat_a = this->material_db.get_material(ii)->get("lattice_constant");												
			// build strains
			RMatrix<double> intrinsic_strain(3,3);
			for(int jj = 0; jj < 3; jj++) {	
				if(strained_axes[jj]) {		 
					intrinsic_strain(jj,jj) = (this->reference_lattice_constant - lat_a) / lat_a;
				} else {
					intrinsic_strain(jj,jj) = 0.0;	
				}
			}
			TDKP_ASSERT(intrinsic_strain(0,1) == intrinsic_strain(1,0), "intrinsic_strain(0,1) == intrinsic_strain(1,0)");
			TDKP_ASSERT(intrinsic_strain(0,2) == intrinsic_strain(2,0), "intrinsic_strain(0,2) == intrinsic_strain(2,0)");	
			TDKP_ASSERT(intrinsic_strain(2,1) == intrinsic_strain(1,2), "intrinsic_strain(2,1) == intrinsic_strain(1,2)");		
			ostringstream sout;
			sout << "IntrinsicStrain: the intrinsic strain in material " << this->material_db.get_material_name(ii) << " is:\n"
			     << intrinsic_strain;
			Logger::get_instance()->emit(LOG_INFO, sout.str());
			this->initial_intrinsic_strains.push_back(intrinsic_strain);
		} else {
			ostringstream sout;
			sout << "IntrinsicStrain: the intrinsic strain in material " << this->material_db.get_material_name(ii) << " is 0 as material has no lattice constant!";
			TDKP_LOGMSG(LOG_WARN, sout.str());			
			RMatrix<double> intrinsic_strain(3,3);
			this->initial_intrinsic_strains.push_back(intrinsic_strain);	
		}
	} 			
}

IntrinsicStrain::IntrinsicStrain(const Geometry& geometry_, MaterialDatabase& material_database_)
: LinearProblem<double>(geometry_, material_database_) {
	this->num_equations_per_node   = -1;
	this->sparsity_copy			     = 0; 
	this->reference_lattice_constant = -1.0;
	this->strained_axes[0] 			 = true;
	this->strained_axes[1]           = true;
	this->strained_axes[2]           = true;
		
	this->num_equations_per_node = this->geometry.get_dimension();
	this->set_sparsity_pattern();	
	
}

IntrinsicStrain::~IntrinsicStrain() {
	if(this->sparsity_copy != 0) {
		delete[] this->sparsity_copy; this->sparsity_copy = 0;
	}
}


/** set reference lattice constant
 * 
 * the intrinsic stress is defined via the difference between the lattice constant and this reference 
 * lattice constant. the intrinsic stress is given as
 * \f$ \sigma_{ii}^{in} = (C_{11} + 2\,C_{12}) \frac{a_{0} - a_{d}}{a_{0}} \f$
 * where \f$a_{d}\f$ is the material lattice constant and \f$a_{0}\f$ is the reference lattice constant
 * 
 */
void IntrinsicStrain::set_reference_lattice_constant(double base_lattice) {
	TDKP_ASSERT(base_lattice > 0.0, "base_lattice > 0.0");
	this->reference_lattice_constant = base_lattice; 	
}

/** build sparsity pattern depending on problem dimensionality */
void IntrinsicStrain::set_sparsity_pattern() throw(Exception*) {
	TDKP_ASSERT(this->geometry.get_dimension() > 0 && this->geometry.get_dimension() <= 3, "this->geometry.get_dimension() > 0 && this->geometry.get_dimension() <= 3");
	// allocate 
	if(this->sparsity_copy == 0 || (unsigned int)this->sparsity_length != this->geometry.get_dimension() * this->geometry.get_dimension()) {
		this->sparsity_length = this->geometry.get_dimension() * this->geometry.get_dimension();			
		if(this->sparsity_copy == 0) {
			delete[] this->sparsity_copy;	
		}
		this->sparsity_copy = new int[this->sparsity_length * 2];
		TDKP_POINTER_ASSERT(this->sparsity_copy);
	}
	for(unsigned int ii = 0; ii < this->geometry.get_dimension(); ii++) {
		for(unsigned int jj = 0; jj < this->geometry.get_dimension(); jj++) {
			this->sparsity_copy[(ii * this->geometry.get_dimension() + jj) * 2] = ii;	
			this->sparsity_copy[(ii * this->geometry.get_dimension() + jj) * 2 + 1] = jj;			
		}		
	} 			
}

const int* IntrinsicStrain::get_node_sparsity_pattern(int &num) const {
	TDKP_BOUNDS_ASSERT(this->sparsity_length > 0, "this->sparsity_length > 0");
	num = this->sparsity_length;
	return this->sparsity_copy;
}


/** calculate the local element matrices (stiff) and load 
 * 
 * stiff and mass is stored in lhs, overlap in rhs
 * lhs and rhs must be min. of length (elem->get_num_nodes()^2 * nonzeros in interaction matrix)
 * in node_internal_indices, the nonzero node indices of the ELEMENT are stored (not globals)
 * means: if there are 4 nodes ([0 ..3]) and node 2 is on the boundary (internal index == -1)
 * then node_internal_indicies will be [0 1 3].
 *
 * @param elem the element where the matrices should be calculated
 * @param lhs array where the element matrix will be stored: lhs[(ii * kpn + jj) * num_sparse + sparse_idx]
 * @param rhs array where the element overlap matrix will be stored
 * @param internal_idx the nonzero indices -> lhs[ ii * kpn + jj 
 */ 
void IntrinsicStrain::calculate_element_matrices(const Element* elem, double* lhs, double *rhs, int* internal_idx, int &n) const {
	
	TDKP_ASSERT(this->properties.size() <= (unsigned)this->material_db.get_num_materials(),"this->properties.size() <= (unsigned)this->material_db.get_num_materials()");
	TDKP_ASSERT(this->geometry.get_dimension() > 0 && this->geometry.get_dimension() < 4, "this->geometry.get_dimension() > 0 && this->geometry.get_dimension() < 4");
	
	double lstiff[3][3][Element::max_num_nodes][Element::max_num_nodes]; /* stiff matrix */
	double lmass[Element::max_num_nodes][Element::max_num_nodes]; 	   /* mass matrix */	
	int    nnode;              /* number of nodes */
	int    neq;                /* number of kp equations */
	int    lsize; 			   /* size of lhs, rhs */
	int    nsparse;            /* number of sparse interaction matrix elements */
	                            	
	const region_properties& prop = this->properties[elem->get_region().get_material().get_id()];				
	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	neq   = this->get_num_equations_per_node();
	lsize = nnode * nnode * neq * neq;
	
	const int* sparse_pat = this->get_node_sparsity_pattern(nsparse);
	
	TDKP_ASSERT(nnode <= Element::max_num_nodes, "nnode <= Element::max_num_nodes (working arrays to small ...)");
		
	// ------------------------------------------------------------------------------	
	// build node_internal_indices (nonzero nodes, which we really calculate)
	// ------------------------------------------------------------------------------	
	for(int ii = 0; ii < nnode; ii++) {
		if(elem->get_node(ii).get_index_internal() != -1) {
			internal_idx[n++] = ii;	
		}	
	}
	// ------------------------------------------------------------------------------	
	// set lhs and rhs to zero		
	// ------------------------------------------------------------------------------	
	for(int ii = 0; ii < lsize; ii++) {
		lhs[ii] = 0.0;
		rhs[ii] = 0.0;	
	}
	// ------------------------------------------------------------------------------			
	// calculate local stiffness and mass matrix for nonzero indices
	// ------------------------------------------------------------------------------	
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {			
			// build stiffness matrix for all derivatives
			for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
				for(unsigned int bb = 0; bb < this->geometry.get_dimension(); bb++) {
					lstiff[aa][bb][ii][jj] = elem->get_element_integral_2nd_order(aa, internal_idx[ii], bb, internal_idx[jj]);								
				}	
			}
			lmass[ii][jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
		}	
	}	

	TDKP_ASSERT((unsigned int)nsparse == this->geometry.get_dimension() * this->geometry.get_dimension(), "nsparse test");
	
			
	// ------------------------------------------------------------------------------	
	// assemble ....
	// ------------------------------------------------------------------------------
	int offset = 0;
	int test;
	// for all nonzero nodes in element
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {			
			offset = (ii * n + jj) * nsparse;	
			// assemble "stiffness" contributions
			// for all differential operators (aa,bb)
			for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
				for(unsigned int bb = 0; bb < this->geometry.get_dimension(); bb++) {
					// for the 2 or 3 equations per node										
					for(unsigned int gg = 0; gg < this->geometry.get_dimension(); gg++) {
						for(unsigned int hh = 0; hh < this->geometry.get_dimension(); hh++) {
							test = gg * this->geometry.get_dimension() + hh;
							TDKP_BOUNDS_ASSERT(sparse_pat[test * 2] == (signed)gg && sparse_pat[test * 2 + 1] == (signed)hh, "sparsepat test");	
							lhs[offset + gg * this->geometry.get_dimension() + hh] += prop.coefficent_matrices[aa][bb][gg][hh] *  lstiff[aa][bb][ii][jj];
						}
					}							
				}
			}
		}	
		// ------------------------------------------------
		// load vector
		// ------------------------------------------------
		for(int gg = 0; gg < neq; gg++) {
			rhs[ii * neq + gg] = 0.0;	
			for(unsigned int ff = 0; ff < this->geometry.get_dimension(); ff++) {
				rhs[ii * neq + gg] += - prop.intrinsic_stress_load(gg,ff) * elem->get_single_integral_1st_order(ff,ii);
			}
		} 
	}					
}

void IntrinsicStrain::build_faky_solution() {
	
	solution.assign(this->geometry.get_num_nonzero_nodes() * geometry.get_dimension(), 0.0);
	// ux = (x - x0), uy = (y - y0) (x0,y0) = 0
	unsigned int ndim = this->geometry.get_dimension();
	for(unsigned int ii = 0; ii < this->geometry.get_num_nodes(); ii++) {
		if(this->geometry.get_node(ii).get_index_internal() != -1) {
			const Node& node = this->geometry.get_node(ii);
			double x = node.get_coord(0);
			double y = node.get_coord(1);
			double vx = (x*x*x / 3.0 - x * 64);
			double vy = (y*y*y / 3.0 - y * 64);
			double f  = 1.0e-5; 						
			solution[ndim * node.get_index_internal() + 0] = f * vx * (y + 8) * (y - 8); 			                                                
			solution[ndim * node.get_index_internal() + 1] = f * vy * (x + 8) * (x - 8);
		}
	}
		
}
/** evaluate strains 
 *
 * derivatives has derivatives[ii * 3 + jj], where ii is the solution variable index and jj is the derviative index 
 */
void IntrinsicStrain::evaluate_strain(
	const double derivatives[3][3], 
	const RMatrix<double>& intrinsic_strain, 
	StrainTensor& strain_tensor
) const {
	
	int derivative_index[3][3] = {
		{iexx, iexy, iexz},
		{ieyx, ieyy, ieyz},
		{iezx, iezy, iezz}
	};
		
	// build strain tensor
	for(unsigned int aa = 0; aa < this->geometry.get_dimension(); aa++) {
		strain_tensor.set(
			derivative_index[aa][aa],
			(derivatives[aa][aa] + intrinsic_strain(aa,aa))
		);
		for(unsigned int bb = aa + 1; bb < this->geometry.get_dimension(); bb++) {
			strain_tensor.set(
				derivative_index[aa][bb],
				0.5 * (derivatives[aa][bb] + intrinsic_strain(aa,bb) + derivatives[bb][aa] + intrinsic_strain(bb,aa))
			);			
		}	
	}
	// set remaining strain components
	for(unsigned int aa = this->geometry.get_dimension(); aa < 3; aa++) {
		strain_tensor.set(
			derivative_index[aa][aa],
			intrinsic_strain(aa,aa)
		);		
		for(unsigned int bb = aa + 1; bb < this->geometry.get_dimension(); bb++) {
			strain_tensor.set(
				derivative_index[aa][bb],
				0.5 * (intrinsic_strain(aa,bb) + intrinsic_strain(bb,aa))
			);				
		}
	}
}

/** return calculated strain field as newly created object */
StrainField* IntrinsicStrain::get_strain_field() const {
					
	TDKP_ASSERT(solution.size() > 0, "solution was not yet calculated ...");	
	TDKP_ASSERT(solution.size() == this->geometry.get_num_nonzero_nodes() * geometry.get_dimension(), "solution has wrong length");
		
	Geometry::element_const_iterator it; 
	StrainField* field = new StrainField(this->geometry);
	// loop over all elements
	double derivatives[3][3];
	double values[3][Element::max_num_nodes];
	int    neq = this->get_num_equations_per_node();
	TDKP_ASSERT(this->get_num_equations_per_node() == (signed)geometry.get_dimension(), "this->get_num_equations_per_node() == geometry.get_dimension()"); 
	
	// ------------------------------------------------------
	// for all elements in structure
	// ------------------------------------------------------
	for(it = this->geometry.elements_begin(); it != this->geometry.elements_end(); it++) {

		if((*it)->enabled()) {
		
			const RMatrix<double>& intrinsic_strain = this->properties[(*it)->get_region().get_material().get_id()].intrinsic_strain;
						
			// -------------------------------------------------						
			// get nodal displacement values
			// -------------------------------------------------
			for(unsigned int ii = 0; ii < (*it)->get_num_nodes(); ii++) {
				if((*it)->get_node(ii).get_index_internal() != -1) {
					for(int nn = 0; nn < neq; nn++) {
						values[nn][ii] = this->solution[(*it)->get_node(ii).get_index_internal() * neq + nn];
					}
				} else {
					for(int nn = 0; nn < neq; nn++) {
						values[nn][ii] = 0.0;	
					}
				}	
			}					
			// -------------------------------------------------
			// build derivatives at element mid 
			// -------------------------------------------------
			for(unsigned int ee = 0; ee < this->geometry.get_dimension(); ee++) {
				for(unsigned int ff = 0; ff < this->geometry.get_dimension(); ff++) {
					derivatives[ee][ff] = (*it)->differentiate_solution(ff, values[ee]);		
				}	
			}
			
			// -------------------------------------------------
			// evaluate strain
			// ------------------------------------------------- 
			this->evaluate_strain(
				derivatives, 
				intrinsic_strain,
				field->get((*it)->get_index_global())
			);		
			
			// -------------------------------------------------
			// now do this at every node for higher order elements
			// -------------------------------------------------
			if((*it)->get_element_order() > 1) {
				
				// -------------------------------------------------
				// for every node
				// -------------------------------------------------
				vector<double> local_coords;
				for(unsigned int nn = 0; nn < (*it)->get_num_nodes(); nn++) {
					
					// -------------------------------------------------
					// get nodes local coordinates
					// -------------------------------------------------
					(*it)->get_node_local_coords(nn, local_coords);
				
					// -------------------------------------------------
					// build derivatives at nodes local coordinates 
					// -------------------------------------------------
					for(unsigned int ee = 0; ee < this->geometry.get_dimension(); ee++) {
						for(unsigned int ff = 0; ff < this->geometry.get_dimension(); ff++) {
							derivatives[ee][ff] = (*it)->differentiate_solution(ff, values[ee], &local_coords[0]);		
						}	
					}
					
					// -------------------------------------------------
					// evaluate strain
					// ------------------------------------------------- 
					this->evaluate_strain(
						derivatives, 
						intrinsic_strain,
						field->get_nodal_strain((*it)->get_index_global(), nn)						
					);
				}		
			}			
		}
	}	
	
	return field;
		
}

void IntrinsicStrain::dump() {
	StrainField* field  = get_strain_field();
	double* node_data = field->get_node_data(this->geometry);
	InputParser parser;	
	if(this->geometry.get_dimension() == 1) {
		Logger::get_instance()->emit(LOG_INFO, "IntrinsicStrain: writing xy data file");
		ofstream fout("strain_node.dat");
		for(int vv = 0; vv < (signed)this->geometry.get_num_nodes(); vv++) {
			int offset = vv * 6;
			fout << this->geometry.get_node(vv).get_coord(0) << " \t";
			double tr_E = node_data[offset + iexx] 
			            + node_data[offset + ieyy]
			            + node_data[offset + iezz];
			fout << tr_E << " \t";
			for(int ee = 0; ee < 6; ee++) {
				fout << node_data[offset + ee] << " \t";	
			}
			fout << "\n";
				
		} 							
		fout.close();
	} else {
		Logger::get_instance()->emit(LOG_INFO, "IntrinsicStrain: writing writing strain field as node field");
		double* node_data = field->get_node_data(this->geometry);
		StdNodeData<double> node_container(
			field->get_num_data_per_element(),
			this->geometry.get_num_nodes(),
			node_data
		);		
		parser.write_ascii(this->geometry, node_container, "strain_node.dat");
	}	
	delete[] node_data;
	delete field;
}


// -----------------------------------------------
// intrinsic strain wurzite implementation
// -----------------------------------------------
IntrinsicStrainWurzite::IntrinsicStrainWurzite(const Geometry& geometry_, MaterialDatabase& material_database_) 
: IntrinsicStrain(geometry_, material_database_),
  reference_lattice_a(-1.0),
  reference_lattice_c(-1.0),
  DIR_XX(D_DX),
  DIR_YY(D_DY),
  DIR_ZZ(D_DZ) 
{  	  	
}

/** throw error, as we need a and c for wurzite */
void IntrinsicStrainWurzite::set_reference_lattice_constant(double base_lattice) {
	TDKP_GENERAL_EXCEPTION("wurzite crystal needs a and c reference lattice constants!");	
}

/** set reference lattice constants a and c */
void IntrinsicStrainWurzite::set_reference_lattice_constant(const double& lattice_a, const double& lattice_c) {
	reference_lattice_a = lattice_a;
	reference_lattice_c = lattice_c;	
}

/** set axis permutation
 * 
 * initially, the x axis is a (axis 0), y axis is a(axis 1) and z axis is c(axis 2)
 * as we assume that the quantization is along axis 0 one can set the permutation
 * of these axes. so if e.q. a c-plane qw is desired, then
 * set_axes(1,2,0) (a0-axis will be mapped to 1, a1 will be mapped to 2, c2 will be mapped to 0)
 */ 
void IntrinsicStrainWurzite::set_axes(short XX, short YY, short ZZ) {
	TDKP_ASSERT(XX != YY, "XX != YY");
	TDKP_ASSERT(XX != ZZ, "XX != ZZ");
	TDKP_ASSERT(ZZ != YY, "ZZ != YY");
	TDKP_ASSERT(XX >= 0 && XX < 3, "XX between 0 and 2");
	TDKP_ASSERT(YY >= 0 && YY < 3, "YY between 0 and 2");
	TDKP_ASSERT(ZZ >= 0 && ZZ < 3, "ZZ between 0 and 2");
	DIR_XX = XX;
	DIR_YY = YY;
	DIR_ZZ = ZZ;		
}


void IntrinsicStrainWurzite::prepare() {	

	if(this->sparsity_copy == 0 || this->reference_lattice_a <= 0.0 || this->reference_lattice_c <= 0.0) {
		TDKP_GENERAL_EXCEPTION("class is not ready to be prepared(?) - misses reference lattice constant or sparsity pattern");	
	} 
	
	ostringstream sout;
	sout << "IntrinsicStrainWurzite: the system is strained (or unstrained) into the following directions:\n";
	for(int ii = 0; ii < 3; ii++) {
		sout << "direction " << ii << ": "
		     << (strained_axes[ii] ? "strained" : "free")
		     << "\n"; 		 	
	}
	
	short permuted_axes[3] = {DIR_XX, DIR_YY, DIR_ZZ};
	
	sout << "IntrinsicStrainWurzite: the axes are mapped to the following directions:\n";
	for(int ii = 0; ii < 3; ii++) {	
		sout << " axis " << ii << " <- ";
		short target_axis = -1;
		for(int jj = 0; jj < 3; jj++) {
			if(permuted_axes[jj] == ii) {
				target_axis = jj;
			}
		}			
		switch(target_axis) {
			case 0:
				sout << " x-axis (a/m plane)\n";
				break;
			case 1:	
				sout << " y-axis (a/m plane)\n";
				break;
			case 2:
				sout << " z-axis (c-axis)\n"; 
				break;
			default:
				TDKP_GENERAL_EXCEPTION("invalid axis");				
		}		
	}
	Logger::get_instance()->emit(LOG_INFO_DEVEL1, sout.str());
	sout.str("");
	
	this->properties.resize(0);
	// ----------------------------------------------
	// prepare initial intrinsic strains
	// ----------------------------------------------
	this->prepare_initial_intrinsic_strains();
	

	region_properties prop;		
	sout.str("");
	sout << "IntrinsicStrainWurtzite: used elastic coefficents are:\n";
	sout.setf(ios::left);
	sout << setw(15) << "Material" 
	     << setw(10) << "C11"
	     << setw(10) << "C12"
	     << setw(10) << "C13"
	     << setw(10) << "C33"
	     << setw(10) << "CXYXY"
	     << setw(10) << "CXZXZ(C44)";
	     
	double C11,   C22,  C33,  C12,  C13, C23,   C44,  C55, C66;
	C11 = C22 = C33 = C12 = C13 = C23 = C44 = C55 = C66 = -1.0;
	double iC11, iC22, iC33, iC12, iC13, iC23, iC44, iC55, iC66;	
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {		
		// init coefficient matrices to zero
		for(int cc = 0; cc < 3; cc++) {			
			for(int dd = 0; dd < 3; dd++) {
				for(int ff = 0; ff < 3; ff++) {
					for(int gg = 0; gg < 3; gg++) {						
						prop.coefficent_matrices[cc][dd][ff][gg] = 0.0;
					}
				}
			}			
		}			
		// ----------------------------------
		// initial elastic constants
		// ----------------------------------
		iC11 = this->material_db.get_material(ii)->get("elastic_constant_C11");
		iC12 = this->material_db.get_material(ii)->get("elastic_constant_C12");
		iC13 = this->material_db.get_material(ii)->get("elastic_constant_C13");
		iC33 = this->material_db.get_material(ii)->get("elastic_constant_C33");
		// --------------------------------------
		// ATTENTION: C44 in literature IS CXZXZ and CYZYZ
		// BECAUSE CXYXY IS IN HEXAGONAL SYMMETRY GIVEN BY (C11 - C12) / 2
		// TO ENFORCE RADIAL SYMMETRY IN XY PLANE (see landau and lifshitz)
		// MANY BOOKS ARE WRONG ON THIS ...
		// --------------------------------------
		iC44 = (iC11 - iC12) / 2; 
		iC66 = this->material_db.get_material(ii)->get("elastic_constant_C44");		
		iC55 = iC66;
		iC22 = iC11; // our assumptions (x/y are equivalent)
		iC23 = iC13; // only have that for permutations
		
		// -----------------------------------
		// write constants to screen
		// ----------------------------------- 
		sout << "\n" 
		     << setw(15) << this->material_db.get_material_name(ii)
		     << setw(10) << iC11
		     << setw(10) << iC12
		     << setw(10) << iC13 
		     << setw(10) << iC33 
		     << setw(10) << iC44 
		     << setw(10) << iC66; 		     
		
		// ---------------------------------
		// remap directions
		// ---------------------------------
		if(DIR_XX == D_DX) {			
			if(DIR_YY == D_DY) {
				// XYZ -> do permutations
				C11 = iC11;
				C12 = iC12;
				C13 = iC13;
				C23 = iC23;
				C22 = iC22;
				C33 = iC33;
				C44 = iC44;
				C55 = iC55;
				C66 = iC66;				
			} else {
				// XZY 
				C11 = iC11;
				C12 = iC13;
				C13 = iC12;
				C22 = iC33;
				C23 = iC23;
				C33 = iC22;
				C44 = iC55;
				C55 = iC44;
				C66 = iC66;	
			}
		} else if(DIR_XX == D_DY) {
			if(DIR_YY == D_DX) {
				// YXZ
				C11 = iC22;
				C12 = iC12;
				C13 = iC23;
				C23 = iC13;
				C22 = iC11;
				C33 = iC33;
				C44 = iC44;
				C55 = iC66;
				C66 = iC55;	
			} else {
				// ZXY
				C11 = iC33;
				C12 = iC13;
				C13 = iC23;
				C22 = iC11;
				C23 = iC12;
				C33 = iC22;
				C44 = iC55;
				C55 = iC66;
				C66 = iC44;
			}
		} else if(DIR_XX == D_DZ) {
			if(DIR_YY == D_DX) {
				// YZX
				C11 = iC22;
				C12 = iC23;
				C13 = iC12;
				C22 = iC33;
				C23 = iC13;
				C33 = iC11;
				C44 = iC66;
				C55 = iC44;
				C66 = iC55;				
			} else {
				// ZYX	
				C11 = iC33;
				C12 = iC23;
				C13 = iC13;
				C22 = iC22;
				C23 = iC12;
				C33 = iC11;
				C44 = iC66;
				C55 = iC55;
				C66 = iC44;
			}
		} else {
			TDKP_GENERAL_EXCEPTION("very strange direction");	
		}
			
		TDKP_ASSERT(C11 > 0.0, "C11");
		TDKP_ASSERT(C22 > 0.0, "C22");		
		TDKP_ASSERT(C33 > 0.0, "C33");		
		TDKP_ASSERT(C12 > 0.0, "C12");		
		TDKP_ASSERT(C13 > 0.0, "C13");		
		TDKP_ASSERT(C23 > 0.0, "C23");		
		TDKP_ASSERT(C44 > 0.0, "C44");		
		TDKP_ASSERT(C55 > 0.0, "C55");		
		TDKP_ASSERT(C66 > 0.0, "C66");		
														
		// D_DX D_DX
		prop.coefficent_matrices[D_DX][D_DX][0][0] = C11;			
		prop.coefficent_matrices[D_DX][D_DX][1][1] = C44;
		prop.coefficent_matrices[D_DX][D_DX][2][2] = C55;
		// D_DY D_DY
		prop.coefficent_matrices[D_DY][D_DY][0][0] = C44;			
		prop.coefficent_matrices[D_DY][D_DY][1][1] = C22;
		prop.coefficent_matrices[D_DY][D_DY][2][2] = C66;			 
		// D_DZ D_DZ
		prop.coefficent_matrices[D_DZ][D_DZ][0][0] = C55;			
		prop.coefficent_matrices[D_DZ][D_DZ][1][1] = C66;
		prop.coefficent_matrices[D_DZ][D_DZ][2][2] = C33;			 			
		// D_DX D_DY
		prop.coefficent_matrices[D_DX][D_DY][0][1] = C12;			
		prop.coefficent_matrices[D_DX][D_DY][1][0] = C44;						
		// D_DY D_DX
		prop.coefficent_matrices[D_DY][D_DX][0][1] = C44;
		prop.coefficent_matrices[D_DY][D_DX][1][0] = C12;
		// D_DX D_DZ
		prop.coefficent_matrices[D_DX][D_DZ][0][2] = C13;
		prop.coefficent_matrices[D_DX][D_DZ][2][0] = C55;
		// D_DZ D_DX
		prop.coefficent_matrices[D_DZ][D_DX][0][2] = C55;
		prop.coefficent_matrices[D_DZ][D_DX][2][0] = C13;			
		// D_DY D_DZ
		prop.coefficent_matrices[D_DY][D_DZ][1][2] = C23;
		prop.coefficent_matrices[D_DY][D_DZ][2][1] = C66;
		// D_DZ D_DY
		prop.coefficent_matrices[D_DZ][D_DY][1][2] = C66;
		prop.coefficent_matrices[D_DZ][D_DY][2][1] = C23;			
						
		if(Configuration::get_instance()->get("output_strain_coefficent_matrices") == 1.0) {
			ostringstream mout;
			mout << "IntrinsicStrainWurzite: elastic coefficent matrices in material " 
			     << this->material_db.get_material_name(ii) << ":\n"; 								
			for(int uu = 0; uu < 3; uu++) {
				for(int jj = 0; jj < 3; jj++) {
					for(int kk = 0; kk < 3; kk++) {
						for(int ll = 0; ll < 3; ll++) {
							mout << "[" << uu << "][" << jj << "][" << kk 
								 << "][" << ll << "] = "
								 <<	prop.coefficent_matrices[uu][jj][kk][ll]
								 << "\n";				
						}	
					}	
				}
			}
			Logger::get_instance()->emit(LOG_INFO, mout.str());		
		}								
											
		// build strains
		prop.intrinsic_strain = this->initial_intrinsic_strains[ii];
				
		// build stresses loads
		for(int jj = 0; jj < 3; jj++) {
			int next = (jj + 1) % 3;
			int prev = (jj + 2) % 3;			
			// diagonal contribution
			prop.intrinsic_stress_load(jj, jj) = 
				prop.coefficent_matrices[jj][jj][jj][jj]     
			  * prop.intrinsic_strain(jj,jj) 
			  + prop.coefficent_matrices[jj][prev][jj][prev] 
			  * prop.intrinsic_strain(prev,prev) 
			  + prop.coefficent_matrices[jj][next][jj][next]
			  * prop.intrinsic_strain(next,next) 
			  + prop.coefficent_matrices[prev][prev][jj][jj]
			  * prop.intrinsic_strain(prev,jj)   
			  + prop.coefficent_matrices[next][next][jj][jj]
			  * prop.intrinsic_strain(next,jj);
			  /*
			cout << "s" << jj << " = " 
			     << prop.coefficent_matrices[jj][jj][jj][jj] << "C" << "* " 
			     << prop.intrinsic_strain(jj,jj) << " e" << jj << " + "
			     << prop.coefficent_matrices[jj][prev][jj][prev] << "C" << "* " 
			     << prop.intrinsic_strain(prev,prev)  << " e" << prev << " + "			        
			     << prop.coefficent_matrices[jj][next][jj][next] << "C" << "* " 
			     << prop.intrinsic_strain(next,next)  << " e" << next << " + "			        
			     << prop.coefficent_matrices[prev][prev][jj][jj] << "C" << "* " 
			     << prop.intrinsic_strain(prev,jj)  << " e" << prev << "_" << jj << " + "			        
			     << prop.coefficent_matrices[next][next][jj][jj] << "C" << "* " 
			     << prop.intrinsic_strain(next,jj)  << " e" << next << "_" << jj << "\n";			        
			*/			  			  
			// shear / offdiagonal				       
			prop.intrinsic_stress_load(jj,prev) = 
				prop.coefficent_matrices[prev][prev][jj][jj] 
			  * prop.intrinsic_strain(jj,prev);
			prop.intrinsic_stress_load(jj,next) = 
				prop.coefficent_matrices[next][next][jj][jj]
			  * prop.intrinsic_strain(jj,next);
		} 
		// store for later						
		this->properties.push_back(prop);
	}	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	
}
void IntrinsicStrainWurzite::prepare_initial_intrinsic_strains() {
	this->initial_intrinsic_strains.resize(0);	
	
	double reference_lattice[3];
	double material_lattice[3];
	
	reference_lattice[0] = DIR_ZZ == D_DX ? reference_lattice_c : reference_lattice_a;
	reference_lattice[1] = DIR_ZZ == D_DY ? reference_lattice_c : reference_lattice_a;
	reference_lattice[2] = DIR_ZZ == D_DZ ? reference_lattice_c : reference_lattice_a;
	
	ostringstream sout;
	sout << "IntrinsicStrainWurzite: reference lattice constants are:\n"
	     << " x = " << reference_lattice[0] << ", y = "  << reference_lattice[1] 
	     << ", z = " << reference_lattice[2]; 
	TDKP_LOGMSG(LOG_INFO_DEVEL1, sout.str());
		
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		
		// skip materials without valid lattice constants
		if(this->material_db.get_material(ii)->valid_key("lattice_constant_a") && this->material_db.get_material(ii)->valid_key("lattice_constant_c")) {		
									
			// see comment on set_reference_lattice_constant			
			const double lat_a = this->material_db.get_material(ii)->get("lattice_constant_a");
			const double lat_c = this->material_db.get_material(ii)->get("lattice_constant_c");
	
			material_lattice[0] = D_DX == DIR_ZZ ? lat_c : lat_a;
			material_lattice[1] = D_DY == DIR_ZZ ? lat_c : lat_a;
			material_lattice[2] = D_DZ == DIR_ZZ ? lat_c : lat_a;
	
			sout.str("");
																					
			// build strains		
			RMatrix<double> intrinsic_strain(3,3);
			for(int jj = 0; jj < 3; jj++) {	
				if(strained_axes[jj]) {		 
					intrinsic_strain(jj,jj) = (reference_lattice[jj] - material_lattice[jj]) / reference_lattice[jj];
				} else {
					intrinsic_strain(jj,jj) = 0.0;	
				}
			}
			TDKP_ASSERT(intrinsic_strain(0,1) == intrinsic_strain(1,0), "intrinsic_strain(0,1) == intrinsic_strain(1,0)");
			TDKP_ASSERT(intrinsic_strain(0,2) == intrinsic_strain(2,0), "intrinsic_strain(0,2) == intrinsic_strain(2,0)");	
			TDKP_ASSERT(intrinsic_strain(2,1) == intrinsic_strain(1,2), "intrinsic_strain(2,1) == intrinsic_strain(1,2)");		
			sout.str("");
			sout << "IntrinsicStrainWurzite: the intrinsic strain in material " << this->material_db.get_material_name(ii) << " is:\n"
			     << intrinsic_strain;
			Logger::get_instance()->emit(LOG_INFO, sout.str());
			this->initial_intrinsic_strains.push_back(intrinsic_strain);

		} else {
			ostringstream sout;
			sout << "IntrinsicStrainWurzite: the intrinsic strain in material " << this->material_db.get_material_name(ii) << " is 0 as material has no lattice constant!";
			TDKP_LOGMSG(LOG_WARN, sout.str());			
			RMatrix<double> intrinsic_strain(3,3);
			this->initial_intrinsic_strains.push_back(intrinsic_strain);	
		}		
		
	} 			
}

} // end of namespace
