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

#include "CoulombMatrixElement.h"
#include <omp.h>

namespace tdkp {
	
CoulombMatrixElement::CoulombMatrixElement(unsigned int wavefunction_length_)
: wavefunction_length(wavefunction_length_) 
{
}
CoulombMatrixElement::~CoulombMatrixElement() {}
	

/** remove requests and results */
void CoulombMatrixElement::clear_requests_and_results() {
	TDKP_TRACE("clearing requests");
 	container.requests.clear();	
 	container.coulomb_matrix_elements.clear();
}
/** only remove wavefunctions */
void CoulombMatrixElement::clear_wavefunctions() {
	wavefunctions.clear();
	container.wavefunction_types.clear();
}
/** remove results, requests and wavefunctions */
void CoulombMatrixElement::clear_all() {
	clear_wavefunctions();
	clear_requests_and_results();
}

/** set new q space discretization
 *
 * removes results 
 */
void CoulombMatrixElement::set_q_range(double q_min, double q_max, unsigned int num_q_values, double progression) {
	TDKP_ASSERT(num_q_values > 0, "num_q_values > 0");	 
	double dq = 1.0;
	if(progression == 1.0) {
		if(num_q_values > 0) {
			dq = (q_max - q_min) / static_cast<double>(num_q_values - 1);
		}
		container.q_values.assign(num_q_values, 0.0);	
		for(unsigned int ii = 0; ii < num_q_values; ii++) {
			container.q_values[ii] = q_min + dq * static_cast<double>(ii);			
		}
	} else {		
		create_nonuniform_point_distribution(
			container.q_values, q_min, q_max, num_q_values, progression
		);
	}	
}
/** set new q space discretization
 *
 * any result is removed 
 */
void CoulombMatrixElement::set_q_range(vector<double>& q_values_) {	
	container.q_values = q_values_;
	sort(container.q_values.begin(), container.q_values.end());
	container.coulomb_matrix_elements.clear();		
}
/** add wavefunctions 
 *
 * @return index of the passed wavefunction 
 */	
unsigned int CoulombMatrixElement::add_wavefunction(const EigenSolution<cplx>* wave, const char* name) {
	
	TDKP_ASSERT(wave->get_length() == (signed)wavefunction_length, "wave->get_length() == wavefunction_length"); 
	wavefunctions.push_back(wave);
	container.wavefunction_types.push_back(string(name));
	return wavefunctions.size() - 1;
		
}

void CoulombMatrixElement::request_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4) {

	// some checks
	TDKP_ASSERT(container.coulomb_matrix_elements.size() == 0, "you have to clear the object before you can reuse it!");
	TDKP_ASSERT(n1 < wavefunctions.size(), "index n1 is out of range: " << n1 << " max, " << wavefunctions.size());
	TDKP_ASSERT(n2 < wavefunctions.size(), "index n2 is out of range: " << n2 << " max, " << wavefunctions.size());
	TDKP_ASSERT(n3 < wavefunctions.size(), "index n3 is out of range: " << n3 << " max, " << wavefunctions.size());
	TDKP_ASSERT(n4 < wavefunctions.size(), "index n4 is out of range: " << n4 << " max, " << wavefunctions.size());
	// one action ;-) 	
	container.requests.push_back(CMEIdx(n1,n2,n3,n4));	
}


const vector<cplx>& CoulombMatrixElement::get_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4) {
	return container.get_matrix_element(n1,n2,n3,n4);	
}

/** return true of this is the same idx */
bool CMEIdx::operator==(const CMEIdx& rhs) const {
	return (n1 == rhs.n1 && n2 == rhs.n2 && n3 == rhs.n3 && n4 == rhs.n4);	
}

/** return true of this is the same idx */
bool CoulombMatrixElement::TwoInt::operator==(const TwoInt& rhs) const {
	return (n1 == rhs.n1 && n4 == rhs.n4);	
}

const EigenSolution<cplx>* CoulombMatrixElement::get_wavefunction(unsigned int idx) const {
	TDKP_BOUNDS_ASSERT(idx < wavefunctions.size(), "idx < wavefunctions.size()");
	return wavefunctions[idx]; 	
}

// -------------------------------------------------------
// implementation of coulomb matrix element quantized 
// -------------------------------------------------------

CoulombMatrixElementQuantized::CoulombMatrixElementQuantized(
	const Geometry& geometry_, 
	MaterialDatabase& material_database_,
	unsigned int kp_basis_size)
: CoulombMatrixElement(geometry_.get_num_nodes()),
  geometry(geometry_), 
  material_database(material_database_),
  bloch_strategy(0)
{
	switch(geometry.get_dimension()) {
		case 1:
			bloch_strategy = new CoulombBlochStrategyWell(kp_basis_size);
			break;
		case 2:
			bloch_strategy = new CoulombBlochStrategyWire(kp_basis_size);
			break;
		case 3:		
			bloch_strategy = new CoulombBlochStrategyDot(kp_basis_size);
			break;
		default:
			TDKP_GENERAL_EXCEPTION("invalid dimensionality");			
	}
}

CoulombMatrixElementQuantized::~CoulombMatrixElementQuantized() {
	delete bloch_strategy; bloch_strategy = 0;
}


/** calculate coulomb matrix elements 
 * 
 */
void CoulombMatrixElementQuantized::calculate() {
		
	TimeMeasurements::get_instance().start("coulomb matrix element");
		
	// --------------------------------------------------
	// calculate only once
	// note: to reduce computational time, we can restrict the integration
	// to the relevant parts ...
	// --------------------------------------------------
	const int num_nonzero_nodes = geometry.get_num_nonzero_nodes();

	// --------------------------------------------------
	// some tests
	// --------------------------------------------------
	TDKP_ASSERT(container.requests.size() > 0, "requests.size() > 0");
	TDKP_ASSERT(container.q_values.size() > 0, "q_values.size() > 0");
	TDKP_ASSERT(container.q_values[0] > 0.0, "");
	for(unsigned int ii = 0; ii < wavefunctions.size(); ii++) {
		TDKP_ASSERT(wavefunctions[ii]->get_length() == num_nonzero_nodes, "wavefunctions[ii]->get_length() == num_nonzero_nodes"); 	
	}
			 
	// --------------------------------------------------
	// calculate element mid points
	// --------------------------------------------------
	ElementMiddlePoints midpoints(geometry);
	
	// --------------------------------------------------
	// first find out which integrations over z we need (n1/n4 pair)  
	// --------------------------------------------------
	vector<TwoInt> z_integrations;
	for(unsigned int ii = 0; ii < container.requests.size(); ii++) {
		TwoInt n1n4(container.requests[ii].n1, container.requests[ii].n4);
		// -----------------------------------------------------
		// add if pair doesn't exist
		// -----------------------------------------------------
		if(find(z_integrations.begin(), z_integrations.end(), n1n4) == z_integrations.end()) {
			z_integrations.push_back(n1n4);
		}		
	}

	// --------------------------------------------------------
	// allocate space we need to store the n1n4 z integrations
	// --------------------------------------------------------
	vector<vector<cplx> > n1n4_z_integration_results(z_integrations.size());
	for(unsigned int ii = 0; ii < z_integrations.size(); ii++) {
		n1n4_z_integration_results[ii].assign(num_nonzero_nodes, 0.0);	
	}
	container.coulomb_matrix_elements.resize(container.requests.size());
	for(unsigned int ii = 0; ii < container.requests.size(); ii++) {
		container.coulomb_matrix_elements[ii].assign(container.q_values.size(), 0.0e0);	
	}
		
	ostringstream sout;
	sout << "i have to calculate " << container.requests.size() << " coulomb integrals "
	     << "for " << container.q_values.size() << " q-points. "
	     << "this leads to the construction of " 
	     << ((num_nonzero_nodes + container.requests.size()) * container.q_values.size())
	     << " matrices of size " << num_nonzero_nodes;
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
						
		
	#pragma omp parallel default(shared)
	{
		// --------------------------------------------------------
		// build coulomb integrator
		// --------------------------------------------------------		
		CoulombIntegrator<cplx>   integrator_cplx(geometry, material_database, midpoints);
		CoulombIntegrator<double> integrator_double(geometry, material_database, midpoints);
		CoulombBlochStrategy*     bloch_strategy_clone = bloch_strategy->clone_obj();

		// get processors nodes (for first stage)
		int start_node = 0; int end_node = 0;	
		get_omp_range(start_node, end_node, geometry.get_num_nodes());
				
		// get processors requests (integration over z_prime)				
		int start_request = 0; int end_request = 0;	
		get_omp_range(start_request, end_request, container.requests.size());
				
		// progress bar
		int next = 0;
		int cc   = 0;
		const unsigned int total = container.q_values.size() * (end_node - start_node + end_request - start_request);
		if(omp_get_thread_num() == 0 && Logger::get_instance()->get_level() > LOG_INFO) {	
			Logger::get_instance()->init_progress_bar("calculating coulomb matrix element", total);
		}		
		
		// --------------------------------------------------------
		// for every q value
		// --------------------------------------------------------		
		for(unsigned int qq = 0; qq < container.q_values.size(); qq++) {			
			// --------------------------------------------------------
			// so, for everz z_prime point
			// --------------------------------------------------------		
			for(int zp = start_node; zp < end_node; zp++) {				
				const Node& node_z_prime = geometry.get_node(zp);
				if(node_z_prime.get_index_internal() != -1) {
					// prepare coulomb function
					bloch_strategy_clone->prepare_coulomb_function(container.q_values[qq], geometry.get_node(zp).get_coords());
					// build matrix
					integrator_double.build_matrix(bloch_strategy_clone->get_coulomb_function());
					// --------------------------------------------------------
					//  do the calculations of n1/n4 pairs (so we integrate over z for every z')
					// --------------------------------------------------------	
					for(unsigned int ii = 0; ii < z_integrations.size(); ii++) {
						const TwoInt& tt = z_integrations[ii]; 					
						// calculate int dz n1*(z) n4(nz) v(z,z') 					
						n1n4_z_integration_results[ii][node_z_prime.get_index_internal()]
						  = bloch_strategy_clone->evaluate_M_product(
						  		this->get_wavefunction(tt.n1),
						  		this->get_wavefunction(tt.n4),
						  		integrator_double
						   	);
						if(isnan(n1n4_z_integration_results[ii][node_z_prime.get_index_internal()].real()) || isnan(n1n4_z_integration_results[ii][node_z_prime.get_index_internal()].imag())) {
							TDKP_TRACE("q: " << qq << " level 1 nan detected set " << ii << " zp" << zp << " node: " << node_z_prime.get_index_internal());	
						}						   	
						   	
					}								
				}
				if(omp_get_thread_num() == 0) {
					if(next == cc && Logger::get_instance()->get_level() > LOG_INFO) {
						next = Logger::get_instance()->set_progress_bar(cc, total);
						TimeMeasurements::get_instance().track_memory_usage();	
					}
					cc++;					
				}								
			}
	
			#pragma omp barrier
	
			// -------------------------------------------------------
			// and now, for every n2/n3 combination
			// -------------------------------------------------------
			for(int rr = start_request; rr < end_request; rr++) {
				// --------------------------------------------------------
				// find appropriate n1n4 product
				// --------------------------------------------------------
				int idx = -1;
				for(unsigned int ii = 0; ii < z_integrations.size(); ii++) {
					TwoInt n1n4(container.requests[rr].n1,container.requests[rr].n4);
					if(z_integrations[ii] == n1n4) {
						idx = ii;	
					}	
				} 
				TDKP_ASSERT(idx != -1, "idx != -1");
				// --------------------------------------------------------
				// build matrix and finally evaluate coulomb integral
				// --------------------------------------------------------
				integrator_cplx.build_matrix(n1n4_z_integration_results[idx]);								
				container.coulomb_matrix_elements[rr][qq] = bloch_strategy_clone->evaluate_M_product(
					this->get_wavefunction(container.requests[rr].n2),
					this->get_wavefunction(container.requests[rr].n3),
					integrator_cplx
				);
				if(isnan(container.coulomb_matrix_elements[rr][qq].real()) || isnan(container.coulomb_matrix_elements[rr][qq].imag())) {
					TDKP_TRACE("q: " << qq << " level 2 nan detected " << rr);	
				}					
				if(omp_get_thread_num() == 0) {		
					if(next == cc && Logger::get_instance()->get_level() > LOG_INFO) {
						next = Logger::get_instance()->set_progress_bar(cc, total);
						TimeMeasurements::get_instance().track_memory_usage();	
					}
					cc++;
				}				
			}
				
			#pragma omp barrier
			
		}	
		delete bloch_strategy_clone;
	} // end parallel block
	if(Logger::get_instance()->get_level() > LOG_INFO) {
		Logger::get_instance()->end_progress_bar();
	}
	
	TimeMeasurements::get_instance().stop("coulomb matrix element");
	
} 


// --------------------------------------------------------
// implementation of bulk coulomb matrix elements
// --------------------------------------------------------
CoulombMatrixElementBulk::CoulombMatrixElementBulk()
: CoulombMatrixElement(1)
{

}
CoulombMatrixElementBulk::~CoulombMatrixElementBulk() {}

/** calculate coulomb matrix element for bulk crystals */	
void CoulombMatrixElementBulk::calculate() {
			
	TDKP_ASSERT(container.q_values.size() > 0, "");
	TDKP_ASSERT(container.q_values[0] > 0.0, "");
	TDKP_ASSERT(container.requests.size() > 0, "requests.size() > 0");
	
	container.coulomb_matrix_elements.resize(container.requests.size());	
	// -------------------------------------------------
	// for every request
	// -------------------------------------------------
	for(unsigned int rr = 0; rr < container.requests.size(); rr++) {
		
		// ---------------------------------------
		// init result space
		// ---------------------------------------
		container.coulomb_matrix_elements[rr].assign(container.q_values.size(), 0.0);		
		// ---------------------------------------
		// get wavefunctions
		// ---------------------------------------
		const EigenSolution<cplx>* n1 = get_wavefunction(container.requests[rr].n1);
		const EigenSolution<cplx>* n2 = get_wavefunction(container.requests[rr].n2);
		const EigenSolution<cplx>* n3 = get_wavefunction(container.requests[rr].n3);
		const EigenSolution<cplx>* n4 = get_wavefunction(container.requests[rr].n4);

		if(n1->get_basis_size() == n4->get_basis_size() && n2->get_basis_size() == n3->get_basis_size() ) {
			cplx n1n4 = 0.0;
			cplx n2n3 = 0.0;
			// n1n4 g-factor 
			for(unsigned int ii = 0; ii < n1->get_basis_size(); ii++) {
				n1n4 += conj(n1->get_node_value(0,ii)) * n4->get_node_value(0,ii);	
			}
			// n2n3 g-factor
			for(unsigned int ii = 0; ii < n2->get_basis_size(); ii++) {
				n2n3 += conj(n2->get_node_value(0,ii)) * n3->get_node_value(0,ii);	
			}
			for(unsigned int qq = 0; qq < container.q_values.size(); qq++) {		
				// so zeta is coulomb fourier transform times weight
				// but here, we multiplied times *q^2		
				container.coulomb_matrix_elements[rr][qq] = 4.0 * constants::pi 			                                         
				                                          * n1n4 * n2n3; 	
			} 			 			
		}
	} 	
}



// --------------------------------------------------------
// implementation of coulomb matrix element data
// --------------------------------------------------------
const int CoulombMatrixElementData::magic = 19012007;
/** constructor for reading coulomb matrix elements from binary file */
CoulombMatrixElementData::CoulombMatrixElementData(const char* filename) {
	ifstream fin(filename, ios::in | ios::binary);
	if(fin) {
		this->read_binary(fin);
		fin.close();		
	} else {
		TDKP_GENERAL_EXCEPTION("can not read from binary file " << filename);	
	}
	
}
/** write coulomb matrix elements to file */
void CoulombMatrixElementData::write_binary(const char* filename) const {
	ofstream fout(filename, ios::out | ios::binary);
	if(fout) {
		this->write_binary(fout);	
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to binary file " << filename);	
	}
}

/** write coulomb matrix elements to binary file */
void CoulombMatrixElementData::write_binary(ostream& stream) const {
			
	int num_items;
	int str_length;
	// -------------------------------------------
	// write header 
	// -------------------------------------------
	stream.write((char*)&magic, sizeof(int));
	// -------------------------------------------
	// write wavefunction types
	// -------------------------------------------
	num_items = wavefunction_types.size();
	stream.write((char*)&num_items, sizeof(int));
	for(unsigned int ii = 0; ii < wavefunction_types.size(); ii++) {
		// length of char string (including \0)
		str_length = wavefunction_types[ii].size() + 1; 
		stream.write((char*)&str_length, sizeof(int));		
		// write string including null termination
		stream.write((char*)wavefunction_types[ii].c_str(),sizeof(char) * str_length); 	
	}
	// -------------------------------------------
	// write q_values
	// -------------------------------------------
	num_items = q_values.size();
	stream.write((char*)&num_items, sizeof(int));
	stream.write((char*)&q_values[0], sizeof(double) * num_items);
	// -------------------------------------------
	// write requests
	// -------------------------------------------
	num_items = requests.size(); 
	stream.write((char*)&num_items, sizeof(int));
	struct { unsigned int n1,n2,n3,n4; } ntuple;
	for(unsigned int ii = 0; ii < requests.size(); ii++) {
		ntuple.n1 = requests[ii].n1;
		ntuple.n2 = requests[ii].n2;
		ntuple.n3 = requests[ii].n3;
		ntuple.n4 = requests[ii].n4;		
		stream.write((char*)&ntuple, sizeof(ntuple)); 	
	}	
	// -------------------------------------------
	// write coulomb matrix elements
	// -------------------------------------------
	num_items = coulomb_matrix_elements.size();
	stream.write((char*)&num_items, sizeof(int));
	for(unsigned int ii = 0; ii < coulomb_matrix_elements.size(); ii++) {
		num_items = coulomb_matrix_elements[ii].size();
		stream.write((char*)&num_items, sizeof(int));
		stream.write((char*)&coulomb_matrix_elements[ii][0], sizeof(cplx) * coulomb_matrix_elements[ii].size());		
	}
	stream.write((char*)&magic, sizeof(int));
}

/** read coulomb matrix elements from binary file */
void CoulombMatrixElementData::read_binary(istream& stream) { 
	int num_items;
	int str_length;
	int test_magic;
	// -------------------------------------------
	// read header 
	// -------------------------------------------
	stream.read((char*)&test_magic, sizeof(int));
	TDKP_ASSERT(test_magic == magic, "magic key does not match ...");
	// -------------------------------------------
	// read wavefunction types
	// -------------------------------------------	
	stream.read((char*)&num_items, sizeof(int));
	TDKP_ASSERT(num_items >= 0, "num wavefunction types has an invalid number");
	wavefunction_types.resize(num_items);
	for(unsigned int ii = 0; ii < wavefunction_types.size(); ii++) {
		// length of char string (including \0)
		stream.read((char*)&str_length, sizeof(int));
		TDKP_ASSERT(str_length > 0, "");
		char* buf = new char[str_length];
		stream.read((char*)buf, sizeof(char) * str_length);
		TDKP_ASSERT(buf[str_length - 1] == '\0', "");
		wavefunction_types[ii] = string(buf);
		delete[] buf; 	
	}
	// -------------------------------------------
	// read q_values
	// -------------------------------------------
	stream.read((char*)&num_items, sizeof(int));
	TDKP_ASSERT(num_items >= 0, "");
	q_values.resize(num_items);
	stream.read((char*)&q_values[0], sizeof(double) * num_items);
	// -------------------------------------------
	// read requests
	// -------------------------------------------	 
	stream.read((char*)&num_items, sizeof(int));
	TDKP_ASSERT(num_items >= 0, "");
	requests.resize(num_items, CMEIdx(0,0,0,0));
	struct { unsigned int n1,n2,n3,n4; } ntuple;
	for(unsigned int ii = 0; ii < requests.size(); ii++) {
		stream.read((char*)&ntuple, sizeof(ntuple));
		requests[ii] = CMEIdx(ntuple.n1,ntuple.n2,ntuple.n3,ntuple.n4);
	}	
	// -------------------------------------------
	// read coulomb matrix elements
	// -------------------------------------------	
	stream.read((char*)&num_items, sizeof(int));
	coulomb_matrix_elements.resize(num_items);
	for(unsigned int ii = 0; ii < coulomb_matrix_elements.size(); ii++) {
		stream.read((char*)&num_items, sizeof(int));
		coulomb_matrix_elements[ii].resize(num_items);		
		stream.read((char*)&coulomb_matrix_elements[ii][0], sizeof(cplx) * coulomb_matrix_elements[ii].size());		
	}
	stream.read((char*)&test_magic, sizeof(int));
	TDKP_ASSERT(test_magic == magic, "magic key at end does not match ...");
}
/** return matrix element, but throw error if it is not calculated */
const vector<cplx>& CoulombMatrixElementData::get_matrix_element(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4) const {
	
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no results calculated yet");
	TDKP_ASSERT(coulomb_matrix_elements.size() == requests.size(), "coulomb_matrix_elements.size() == requests.size()");
	// ----------------------------------------------
	// find the results
	// ----------------------------------------------
	int idx = -1;
	CMEIdx ureq(n1,n2,n3,n4);
	for(unsigned int ii = 0; ii < requests.size(); ii++) {
		if(ureq == requests[ii]) {
			idx = ii;
			break;	
		}		
	}
	TDKP_ASSERT(idx != -1, "the matrix element <" << n1 << " " << n2 << " | v | " << n3 << " " << n4 << "> has not been calculated");
	return coulomb_matrix_elements[idx];
}

int CoulombMatrixElementData::get_x_length() const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet"); 
	return q_values.size();		
}
int CoulombMatrixElementData::get_num_y_sets() const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	return coulomb_matrix_elements.size();	
}
int CoulombMatrixElementData::get_num_x_sets() const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	return 1;	
} 
void CoulombMatrixElementData::get_x(int xidx, vector<double> &x) const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	x = q_values;
}
void CoulombMatrixElementData::get_y(int yidx, vector<double>& y) const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	TDKP_ASSERT(yidx < (signed)coulomb_matrix_elements.size(), "yidx < coulomb_matrix_elements.size()");
	y.assign(q_values.size(), 0.0);
	for(unsigned int ii = 0; ii < q_values.size(); ii++) {
		y[ii] = tdkp_math::abs(coulomb_matrix_elements[yidx][ii]);
	}
} 
string CoulombMatrixElementData::get_x_identifier(int xidx) const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	return string("q");
}


string CoulombMatrixElementData::get_y_identifier(int yidx) const {
	TDKP_ASSERT(coulomb_matrix_elements.size() > 0, "no result calculated yet");
	TDKP_ASSERT(yidx < (signed)requests.size(), "yidx < request.size()");
	const CMEIdx& r = requests[yidx];
	ostringstream sout;
	sout << "<" 
	     << wavefunction_types[r.n1] << r.n1 << ","
	     << wavefunction_types[r.n2] << r.n2 << "|v|"
	     << wavefunction_types[r.n3] << r.n3 << ","
	     << wavefunction_types[r.n4] << r.n4 << ">";
	return sout.str();	
}


const string& CoulombMatrixElementData::get_wavefunction_type(unsigned int ni) const {
	TDKP_ASSERT(ni < wavefunction_types.size(), "");
	return wavefunction_types[ni];	
}
unsigned int  CoulombMatrixElementData::get_num_wavefunctions() const {
	return wavefunction_types.size();	
}

} // end of namespace
