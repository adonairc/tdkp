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

#include "tdkp/utilities/MatrixElements.h"
#include "tdkp/common/all.h"
#include "tdkp/probdefs/KP8x83D.h"
#include "tdkp/probdefs/KP4x43D.h"
#include "tdkp/probdefs/KP4x41D2D.h"
#include "tdkp/probdefs/KP6x61D2D.h"
#include "tdkp/probdefs/KP8x81D2D.h"
#include "tdkp/probdefs/EffectiveMass.h"

namespace tdkp {

MatrixElements::MatrixElements(MomentumOperator& momentum_operator_)
: momentum_operator(momentum_operator_),
  num_cb_bands(0),
  num_vb_bands(0)
{
}

MatrixElements::~MatrixElements()
{

}

unsigned int MatrixElements::get_num_cb_bands() const {
	return this->num_cb_bands;
}
unsigned int MatrixElements::get_num_vb_bands() const {
	return this->num_vb_bands;
}
unsigned int MatrixElements::get_num_k_values() const {
	return this->domain.get_number_of_points();
}

/** direct data access */
const double& MatrixElements::get_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk) const {
	return this->values_square[get_values_idx(ee,hh,dir,kk)];
}

/** return index where we stored that value */
unsigned int MatrixElements::get_values_idx(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk) const {
	const unsigned int num_k_values = domain.get_number_of_points();
	TDKP_ASSERT(num_k_values > 0, "num_k_values > 0");
	unsigned int offset = kk
	                    + num_k_values * dir
	                    + num_k_values * 3 * hh
	                    + num_k_values * 3 * this->num_vb_bands * ee;
	TDKP_BOUNDS_ASSERT(offset < this->values_square.size(), "index out of bounds: " << offset << ", max = " << values_square.size() << ", call: (" << ee << ", " << hh << ", " << dir << ", " << kk<< ")");
	return offset;
}

/** direct data access */
double& MatrixElements::get_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk) {
	return const_cast<double&>(static_cast<const MatrixElements&>(*this).get_abs_square(ee,hh,dir,kk)); // use const member function and recast
}

void MatrixElements::set_abs_square(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, const double& value) {
	TDKP_ASSERT(tdkp_math::abs(this->get_abs_square(ee,hh,dir,kk)) == 0.0, "tdkp_math::abs(this->get_abs_square(ee,hh,dir,kk)) == 0.0");
	this->get_abs_square(ee,hh,dir,kk) = value;
}

/** return the matrix element <c|p|v> 
 *
 * value_idx may be only 0 for kp8x8. for eff mass type cb, 
 * value_idx = 0 is momentum matrix element with cb spin up and
 * value_idx = 1 is cb spin down 
 * 
 * @param ee cb band index
 * @param hh vb band index
 * @param dir direction of momentum operator
 * @param kk  k index
 * @param value_idx index of requested value  
 */
const complex<double>& MatrixElements::get(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx) const {
	unsigned int idx = (momentum_operator.eval_results_length() - 1) * get_values_idx(ee,hh,dir,kk) + value_idx;
	TDKP_BOUNDS_ASSERT(idx < this->values_raw.size(), "idx < this->values_raw.size()");
	return this->values_raw[idx];	 		
}

/** set the matrix element <c|p|v> (idx value_idx)*/
void MatrixElements::set(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx, const complex<double>& value) {
	unsigned int idx = (momentum_operator.eval_results_length() - 1) * get_values_idx(ee,hh,dir,kk) + value_idx;
	TDKP_BOUNDS_ASSERT(idx < this->values_raw.size(), "idx ( " << idx << ") < this->values_raw.size() (" << this->values_raw.size() << ")");
	this->values_raw[idx] = value;	 	
}

/** return the matrix element <c|p|v> */
complex<double>& MatrixElements::get(unsigned int ee, unsigned int hh, unsigned int dir, unsigned int kk, int value_idx) {
	return const_cast<complex<double>& >(static_cast<const MatrixElements&>(*this).get(ee,hh,dir,kk,value_idx)); // use const member function and recast	
}


/*
MatrixElements::value_iterator
MatrixElements::begin(unsigned int ee, unsigned int hh, unsigned int dir) {
	// iterator plus offset
	return this->values.begin() + get_values_idx(ee,hh,dir,0);
}

MatrixElements::const_value_iterator
MatrixElements::begin(unsigned int ee, unsigned int hh, unsigned int dir) const {
	// iterator plus offset
	return this->values.begin() + get_values_idx(ee,hh,dir,0) + this->get_num_k_values();
}

MatrixElements::const_value_iterator
MatrixElements::end(unsigned int ee, unsigned int hh, unsigned int dir) const {
	// end is offset + num_k_values
	return this->values.begin() + get_values_idx(ee,hh,dir,0) + this->get_num_k_values();
}
*/
int MatrixElements::get_x_length() const {
	TDKP_ASSERT(this->values_square.size() > 0, "no matrix elements calculated!");
	return this->get_num_k_values();
}
int MatrixElements::get_num_y_sets() const {
	return this->get_num_cb_bands()
	     * this->get_num_vb_bands()
	     * 3;
}
void MatrixElements::get_x(int xidx, vector<double> &x) const {
	x.resize(this->get_num_k_values());
	for(unsigned int ii = 0; ii < domain.get_number_of_points(); ii++) {
		x[ii] = domain.get_point(ii).get_coord(xidx);
	}
}

void MatrixElements::get_y(int yidx, std::vector<double>& y) const {

	int dir_idx   = yidx % 3;
	yidx = ((yidx - dir_idx) / 3);
	int hole_band = yidx % this->get_num_cb_bands();
	yidx = (yidx - hole_band) / this->get_num_vb_bands();
	int elec_band = yidx;

	y.resize(this->get_num_k_values());
	for(unsigned int kk = 0; kk < this->get_num_k_values(); kk++) {
		y[kk] = this->get_abs_square(elec_band,hole_band,dir_idx,kk);
	}

}

string MatrixElements::get_x_identifier(int xidx) const {
	ostringstream sout;
	sout << "k" << xidx;
	return string(sout.str());
}
string MatrixElements::get_y_identifier(int yidx) const {
	int dir_idx   = yidx % 3;
	yidx = ((yidx - dir_idx) / 3);
	int hole_band = yidx % this->get_num_vb_bands();
	yidx = (yidx - hole_band) / this->get_num_vb_bands();
	int elec_band = yidx;
	ostringstream str;
	str << "cb" << elec_band << "vb" << hole_band << "p" << dir_idx;
	return str.str();
}

void MatrixElements::calculate(const BandstructureDomain<cplx>& cb_bands, const BandstructureDomain<cplx>& vb_bands) {
	this->calculate(cb_bands.get_number_of_bands(), vb_bands.get_number_of_bands(), cb_bands, vb_bands);
}

/** solve bulk matrix elements within the effective mass approximaton
 *
 * attention, the bulk cb matrix element includes the degeneracy!
 */
void MatrixElements::calculate_bulk_effmass(const BandstructureDomain<complex<double> >& vb_bands) {
	// -------------------------------------------
	// build faky single band bandstructure
	// -------------------------------------------
	DomainMaster single_gamma_point(
		new DomainNodeSingularPoint(
			vb_bands.get_domain().radial(), 0.0, 0.0, 0.0, 0.0
		)
	);
	BandstructureDomain<cplx> cb_bands(
		1, 1, 1, single_gamma_point
	);
	cplx one(1.0, 0.0);
	cb_bands.add_eigensolution(0, 0, new EigenSolution<cplx>(0.0, &one, 1, 1));
	// -------------------------------------------
	// and use standard function ...
	// -------------------------------------------
	this->calculate(cb_bands, vb_bands);
}

/** takes eff mass matrix elements at k = 0 and expands them to the whole domain
 */
MatrixElements MatrixElements::get_disp_matrix_elements(const DomainMaster& new_domain) {
	
	TDKP_ASSERT(this->values_square.size() > 0, "no matrix elements calculated!");
	TDKP_ASSERT(this->domain.get_dimension() == 0, "create_disp_matrix_elements does only make sense for effective mass matrix elements!");
	
	MatrixElements ret(momentum_operator);
	ret.num_cb_bands = num_cb_bands;
	ret.num_vb_bands = num_vb_bands;
	ret.domain       = new_domain;
	 
	ret.values_square.assign(ret.num_cb_bands * ret.num_vb_bands * 3 * ret.get_num_k_values(), 0.0);
	ret.values_raw.assign((ret.momentum_operator.eval_results_length() - 1) * ret.num_cb_bands * ret.num_vb_bands * 3 * ret.get_num_k_values(), 0.0);	 
	 
	for(unsigned int cc = 0; cc < num_cb_bands; cc++) {
		for(unsigned int vv = 0; vv < num_vb_bands; vv++) {
			for(unsigned int pp = 0; pp < 3; pp++) {
				for(unsigned int kk = 0; kk < new_domain.get_number_of_points(); kk++) {	 
					ret.set_abs_square(cc, vv, pp, kk, this->get_abs_square(cc,vv,pp,0));
					for(int vl_idx = 0; vl_idx < momentum_operator.eval_results_length() - 1; vl_idx++) {
						ret.set(cc, vv, pp, kk, vl_idx,  this->get(cc,vv,pp,0,vl_idx));
					}
				}
			}
		}
	}

	return ret;
	
}

/** calculates the matrix elements for the passed bands
 *
 * @max_cb_bands use cb bands [0, max_cb_bands[
 * @max_vb_bands use vb bands [0, max_cb_bands[
 */
void MatrixElements::calculate(unsigned int max_cb_bands, unsigned int max_vb_bands, const BandstructureDomain<cplx>& cb_bands, const BandstructureDomain<cplx>& vb_bands) {

	TimeMeasurements::get_instance().start("matrix element calculation");

	ostringstream sout;

	// --------------------------------------------------
	// test compatibility with momentum operator
	// --------------------------------------------------
	if(!momentum_operator.compatible(cb_bands.get_basis_size(), cb_bands.get_solution_length(), vb_bands.get_basis_size(), vb_bands.get_solution_length())) {
		sout << "your momentum operator "
		     << "(" << momentum_operator.get_description() << ") "
		     << "is not suited for your the matrix element calculation!";
		TDKP_GENERAL_EXCEPTION(sout.str());
	}

	// --------------------------------------------------
	// test compatibility of bandstructures ...
	// but if cb bands are effmass, we accept it anyway
	// -------------------------------------------------
	if(cb_bands.get_basis_size() != 1 && !cb_bands.get_domain().compare_points(vb_bands.get_domain())) {
		sout << "the bandstructure you supplied was not calculated on the same points of the Brillouin zone.";
		TDKP_GENERAL_EXCEPTION(sout.str());
	}

	// --------------------------------------------------
	// get storage space
	// --------------------------------------------------
	this->num_cb_bands = max_cb_bands;
	this->num_vb_bands = max_vb_bands;
	this->domain       = vb_bands.get_domain();
	this->values_square.assign(this->num_cb_bands * this->num_vb_bands * 3 * this->get_num_k_values(), 0.0);
	this->values_raw.assign((momentum_operator.eval_results_length() - 1) * this->num_cb_bands * this->num_vb_bands * 3 * this->get_num_k_values(), 0.0);

	if(!domain.frozen()) {
		TDKP_GENERAL_EXCEPTION("a domain that comes from a bandstructure object must be frozen in order to prevent modifications!");
	}

	// --------------------------------------------------
	// some local variables
	// --------------------------------------------------
	const int cb_basis_size	= cb_bands.get_basis_size();
	const int vb_basis_size	= vb_bands.get_basis_size();
//	vector<cplx> eval_results(momentum_operator.eval_results_length());

	// --------------------------------------------------
	// some user info
	// --------------------------------------------------
	sout << "calculating momentum matrix elements:\n"
	     << "num cb bands:		" << num_cb_bands  << "\n"
	     << "cb basis size:     " << cb_basis_size << "\n"
	     << "num vb bands:      " << num_vb_bands  << "\n"
		 << "vb basis size:     " << vb_basis_size << "\n"
		 << "momentum operator " << (momentum_operator.is_k_dependent() ? "does depend on k" : "is independent of k") << "\n"
		 << "num k points:      " << domain.get_number_of_points();
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
	if(Logger::get_instance()->get_level() == LOG_INFO) {
		Logger::get_instance()->init_progress_bar("calculating matrix elements", 3 * this->get_num_k_values());
	}
	int next, current; next = current = 0;

	// --------------------------------------------------
	// calculate the bastards!
	// --------------------------------------------------
	// for all directions
	for(unsigned short dd = 0; dd < 3; dd++) {
		TimeMeasurements::get_instance().track_memory_usage();
		momentum_operator.set_operator_direction(dd);
		// ----------------------------------------------------
		// build operator for non-k dependent stuff
		// ----------------------------------------------------
		if(!momentum_operator.is_k_dependent()) {
			momentum_operator.lock();
		}

		// ----------------------------------------------------
		// for all k values
		// ----------------------------------------------------
		for(unsigned int kk = 0; kk < this->get_num_k_values(); kk++) {
			if(Logger::get_instance()->get_level() == LOG_INFO) {
				if(current == next) {
					next = Logger::get_instance()->set_progress_bar(dd * this->get_num_k_values() + kk, 3 * this->get_num_k_values());
				}
				current++;
			} else {
				sout.str("");
				sout << "calculating matrix element p" << dd << " at kk = " << kk;
				Logger::get_instance()->emit(LOG_INFO_DEVEL2, sout.str());
			}
			// -----------------------------------------------------
			// tell the k dependent operator the current k point
			// -----------------------------------------------------
			if(momentum_operator.is_k_dependent()) {
				momentum_operator.set_k_value(domain.get_point(kk));
				momentum_operator.lock(); // prepares the object
			}

			// -----------------------------------------------------
			// prepare the k idx for the cb band
			// if cb band is effective mass, we restrict ourself to kk = 0 there
			// -----------------------------------------------------
			const int cb_kidx = cb_bands.get_basis_size() == 1 ? 0 : kk;

			// -----------------------------------------------------
			// loop over all solution pairs
			// -----------------------------------------------------
			#pragma omp parallel for
			for(int vb_idx = 0; vb_idx < (signed)num_vb_bands; vb_idx++) {
				for(unsigned int cb_idx = 0; cb_idx < num_cb_bands; cb_idx++) {
				
					const EigenSolution<cplx>& cb_wave = cb_bands.get_eigensolution(cb_kidx, cb_idx);
					const EigenSolution<cplx>& vb_wave = vb_bands.get_eigensolution(kk, vb_idx);
					vector<cplx> eval_results(momentum_operator.eval_results_length());

					// -----------------------------------------------------
					// evaluate matrix element
					// -----------------------------------------------------
					momentum_operator.eval(eval_results, cb_wave.get_data(), vb_wave.get_data());
					TDKP_ASSERT(eval_results[0].imag() == 0.0, "|<c|p|v>|^2 = eval_results[0].imag() == 0.0 failed (was " << eval_results[0].imag() << " for cb " << cb_idx << ", vb " << vb_idx << " at k = " << kk << " for p" << dd << ")");  
					this->set_abs_square(cb_idx,vb_idx,dd,kk, eval_results[0].real());
					for(unsigned int gg = 0; gg < eval_results.size() - 1; gg++) {
						this->set(cb_idx,vb_idx,dd,kk,gg,eval_results[gg + 1]);	
					}
				}

			}

			// -----------------------------------------------------
			// release momentum operator
			// -----------------------------------------------------
			if(momentum_operator.is_k_dependent()) {
				momentum_operator.release(); // prepares the object
			}

		}

		// ----------------------------------------------------
		// release k independent operator
		// ----------------------------------------------------
		if(!momentum_operator.is_k_dependent()) {
			momentum_operator.release();
		}
	}
	if(Logger::get_instance()->get_level() == LOG_INFO) {
		Logger::get_instance()->end_progress_bar();
	}
	
	TimeMeasurements::get_instance().stop("matrix element calculation");
	
}


/** read matrix elements from binary file */
void MatrixElements::read_binary(const char* filename) {
	ifstream fin(filename, ios::binary);
	if(fin) {
		this->read_binary(fin);
		fin.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not read matrix elements from file " << filename);	
	}		
}

/** write matrix elements to binary file */
void MatrixElements::write_binary(const char* filename) const {
	ofstream fout(filename, ios::binary);
	if(fout) {
		this->write_binary(fout);
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("can not write matrix elements to file " << filename);	
	}			
}

/** read matrix elements from binary stream */
void MatrixElements::read_binary(istream& in) {

	int magic;	
	in.read((char*)&magic, sizeof(int));
	TDKP_ASSERT(magic == 101010, "magic key at matrix element start is wrong");
	
	unsigned int vsqlength, vrwlength;
	
	// ------------------------------------------------
	// read the number of bands, the array length and the array 
	// ------------------------------------------------
	in.read((char*)&num_cb_bands, sizeof(unsigned int));
	in.read((char*)&num_vb_bands, sizeof(unsigned int));
	in.read((char*)&vsqlength,    sizeof(unsigned int));
	in.read((char*)&vrwlength,    sizeof(unsigned int));
	TDKP_ASSERT(vsqlength > 0, "vsqlength > 0");
	TDKP_ASSERT(vrwlength > 0, "vrwlength > 0");
	values_square.resize(vsqlength);
	values_raw.resize(vrwlength);
	in.read((char*)&values_square[0], sizeof(double) * vsqlength);
	in.read((char*)&values_raw[0],    sizeof(complex<double>&) * vrwlength);
	
	// ------------------------------------------------
	// store the domain
	// ------------------------------------------------
	this->domain.read_binary(in);		
	in.read((char*)&magic, sizeof(int));
	TDKP_ASSERT(magic == 101010, "magic key at matrix element end is wrong");
	
}

/** write matrix elements to binary stream */
void MatrixElements::write_binary(ostream& out) const {

	int magic = 101010;
	out.write((char*)&magic, sizeof(int));
	
	unsigned int vsqlength = values_square.size();
	unsigned int vrwlength = values_raw.size();
	
	// ------------------------------------------------
	// store the number of bands, the array length and the array 
	// ------------------------------------------------
	out.write((char*)&num_cb_bands, sizeof(unsigned int));
	out.write((char*)&num_vb_bands, sizeof(unsigned int));
	out.write((char*)&vsqlength,    sizeof(unsigned int));
	out.write((char*)&vrwlength,    sizeof(unsigned int));
	out.write((char*)&values_square[0], sizeof(double) * vsqlength);
	out.write((char*)&values_raw[0],    sizeof(complex<double>) * vrwlength);
	
	// ------------------------------------------------
	// store the domain
	// ------------------------------------------------
	this->domain.write_binary(out);		
	out.write((char*)&magic, sizeof(int));
	
}









// --------------------------------------------------------
// base momentum operator implementation
// --------------------------------------------------------

MomentumOperator::MomentumOperator()
: momentum_operator_direction(0)
{

}
MomentumOperator::~MomentumOperator() {

}



// --------------------------------------------------------
// base bulk kp momentum operator
// --------------------------------------------------------
MomentumOperatorBulk::MomentumOperatorBulk()
: locked(false),
  tensor(0)
{
	kvalues[0] = kvalues[1] = kvalues[2] = 0.0e0;
}
void  MomentumOperatorBulk::release() {
	TDKP_ASSERT(locked, "you can not release an unlocked momentum operator");
	locked = false;
}
void  MomentumOperatorBulk::set_k_value(const DomainPoint& point) {
	TDKP_ASSERT(!locked, "you may not change k values of a locked structure");
	kvalues[0] = point.get_coord(0);
	kvalues[1] = point.get_coord(1);
	kvalues[2] = point.get_coord(2);
}


// --------------------------------------------------------
// bulk kp 8x8 momentum operator
// --------------------------------------------------------
/** constructor for zinc-blende kp 8x8 matrix elements */
MomentumOperatorBulk8x8::MomentumOperatorBulk8x8(const KPMatrix8x8EndersForeman& matrix_)
: matrix(matrix_),
  k_dependent(false)
{
	tensor.resize(64);
	TDKP_ASSERT(matrix.ready(), "matrix.ready() failed - use .calculate() before passing to momentum operator");
}

/** constructor for wurtzite kp 8x8 matrix elements */
MomentumOperatorBulk8x8::MomentumOperatorBulk8x8(const KPMatrix8x8Wurtzite& matrix_)
: matrix(matrix_),
  k_dependent(false)
{
	tensor.resize(64);
	TDKP_ASSERT(matrix.ready(), "matrix.ready() failed - use .calculate() before passing to momentum operator");
}

/** locking means preparing!
 *
 *
 */
void  MomentumOperatorBulk8x8::lock() {
	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;

	cplx i(0.0, 1.0);

	// -------------------------------------------------
	// see enders prb 1995 eq. (7')
	// the momentum matrix element is not only the matrix
	// element between bulk bloch functions, but also contains
	// k dependent stuff coming from folding in of higher
	// bands with k	
	//
	// attention: the factor m / hbar is applied in the
	// eval routine! so don't try to find it here
	// -------------------------------------------------
	
	

	// -------------------------------------------------
	// starting with bulk block momentum stuff
	// -------------------------------------------------
	tensor.assign(64, 0.0e0);

	RMatrix<cplx>* left  = matrix.get_first_order_strainless_matrix(KPMatrixBase::op_left, get_operator_direction());
	RMatrix<cplx>* right = matrix.get_first_order_strainless_matrix(KPMatrixBase::op_right, get_operator_direction());
	for(short ii = 0; ii < 8; ii++) {
		for(short jj = 0; jj < 8; jj++) {
			tensor[ii * 8 + jj] = ((*left)(ii,jj) + (*right)(ii,jj)) * i;
		}
	}
	delete left;
	delete right;

	// ------------------------------------------
	// second order stuff
	// ------------------------------------------
	if(is_k_dependent()) {
		// --------------------------------------------------------
		// o.k. its simple:  i store the p.d.e. matrix
		// dx H dx (but actually the partially integrated, so -H dx dx)
		// in terms of kx -> dx H dx -> -kx H kx, but as H already -H,
		// i can omit the negative sign
		// --------------------------------------------------------
		for(unsigned int short ee = 0; ee < 3; ee++) {
			RMatrix<cplx>* left  = matrix.get_second_order_strainless_matrix(get_operator_direction(), ee);
			RMatrix<cplx>* right = matrix.get_second_order_strainless_matrix(ee, get_operator_direction());
			for(short ii = 0; ii < 8; ii++) {
				for(short jj = 0; jj < 8; jj++) {
					tensor[ii * 8 + jj] += 	((*left)(ii,jj) + (*right)(ii,jj)) * kvalues[ee];
				}
			}
			delete left;
			delete right;
		}
	}
		
}



bool  MomentumOperatorBulk8x8::compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const {
	if(cb_basis_size == 8 && cb_length == 1 && vb_basis_size == 8 && vb_length == 1) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << 8 << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << 8 << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << 1 << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << 1;
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;
}


/** returns |<cb|p|vb>|^2 and <cb|p|vb> */
void MomentumOperatorBulk8x8::eval(vector<cplx>& result_ret, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {
	TDKP_BOUNDS_ASSERT(cb_wave.size() == vb_wave.size() && cb_wave.size() == 8, "cb_wave.size() == vb_wave.size() && cb_wave.size() == 8");
	complex<double> result(0.0);
	complex<double> tmp;
		
	// matrix vector multiplication ... conj(transpose(cb)) * p * vb
	for(unsigned int jj = 0; jj < 8; jj++) {
		tmp = 0.0;
		for(unsigned int kk = 0; kk < 8; kk++) {
			tmp += tensor[jj * 8 + kk] * vb_wave[kk];
		}
		result += conj(cb_wave[jj]) * tmp;
	}	
	// <c|p|v>
	tmp = constants::m0 / constants::hbar * result;
	result_ret[1] = tmp;
	// |<c|p|v>|^2
	result_ret[0] = tmp.real() * tmp.real() + tmp.imag() * tmp.imag();	 
	
}

/** constructor for bulk 6x6 matrix */
MomentumOperatorBulk6x6::MomentumOperatorBulk6x6(const double& optical_matrix_parameter)
: optical_matrix_param(optical_matrix_parameter)
{
	tensor.resize(6);
}

MomentumOperatorBulk6x6::MomentumOperatorBulk6x6(const Material& material)
: optical_matrix_param(material.get("optical_matrix_element"))
{
	tensor.assign(6, 0.0);
}

/** locking momentom operator assembles it */
void  MomentumOperatorBulk6x6::lock() {

	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;

	// ----------------------------------------
	// see eval for details
	// P = sqrt(hbar2_2m0 * Ep)
	// -> |<c|p|v>|^2 = m^2 / hbar^2 * hbar^2 / 2m0 * Ep
	//                -> m0 / 2.0
	// ----------------------------------------
	tensor.assign(6, 0.0);
	tensor[get_operator_direction()] = tensor[get_operator_direction() + 3]
	                                 = constants::m0 / 2.0 * optical_matrix_param;
}

/** internal test of class matrix elements whether the supplied bands could be appropriate */
bool  MomentumOperatorBulk6x6::compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const {
	if(cb_basis_size == 1 && cb_length == 1 && vb_basis_size == 6 && vb_length == 1) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << 1 << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << 6 << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << 1 << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << 1;
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;
}

/** evaluate bulk momentum operator */
void MomentumOperatorBulk6x6::eval(vector<cplx>& result, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	TDKP_BOUNDS_ASSERT(cb_wave.size() == 1 && vb_wave.size() == 6, "cb_wave.size() == 1 && vb_wave.size() == 6");

	// ----------------------------------------------------
	// o.k. the cb basis function has no spin index, so
	// we assume |S down> and |S up>
	// means: we include the spin degeneracy in here
	// when returning the square values. for the non-squared
	// values, we return both, for  |S up> and |S down> 
	// ----------------------------------------------------
	// therefore, in the 6x6 matrix, we have
	// the basis ordered as |Xu> |Yu> |Zu> |Xd> |Yd> |Zd>
	// and |S> is assumed to be == 1
	// therefore, e.g. px:
	// |<cb|px|vb>|^2 = |<Su|px|Xu>|^2 + |<Sd|px|Xd>|^2
	//                = m/hbar |P<fSu|fXu>|^2 + |P<fSd|fXd>|^2
	//                = m/hbar P^2 * (|<fSu|fXu>|^2 + spin down)
	// ----------------------------------------------------
	double ret = 0.0;
	cplx tmp(0.0);
	cplx i(0.0,1.0);
	for(short ii = 0; ii < 3; ii++) {
		ret += tensor[ii].real() * (vb_wave[ii].real() * vb_wave[ii].real() + vb_wave[ii].imag() * vb_wave[ii].imag());
		tmp += i * sqrt(tensor[ii]) * vb_wave[ii];
	}
	result[1] = tmp;
	tmp = 0.0e0;
	for(short ii = 3; ii < 6; ii++) {
		ret += tensor[ii].real() * (vb_wave[ii].real() * vb_wave[ii].real() + vb_wave[ii].imag() * vb_wave[ii].imag());
		tmp += i * sqrt(tensor[ii]) * vb_wave[ii];
	}	
	result[0] = ret;
	result[2] = tmp;
	
	tmp = conj(result[1]) * result[1] + conj(result[2]) * result[2]; 
	TDKP_BOUNDS_ASSERT(tdkp_math::abs(result[0] - tmp) < 1.0e-10, "tdkp_math::abs(result[0] - tmp) < 1.0e-10");   

}
/** constructor for 6x6 bulk wurtzite matrix elements
 * 
 * no strain dependence available
 * P1 and P2 (optical momentum matrix elements) are calculated using the
 * designated routines of the 8x8 wurtzite matrix which use chuangs
 * equation (18)
 */
MomentumOperatorBulk6x6WZ::MomentumOperatorBulk6x6WZ(const Material& material)
: P1(0.0), P2(0.0)
{
	tensor.assign(6, 0.0);
	KPMatrix8x8Wurtzite kp_matrix(&material);
	kp_matrix.calculate();
	// calculate P1 and P2 
        if(material.is_set("optical_momentum_matrix_element_P1")) {
                P1 = material.get("optical_momentum_matrix_element_P1");
         } else {
                P1 = kp_matrix.get_optical_momentum_matrix_element_P1();
        }
        if(material.is_set("optical_momentum_matrix_element_P2")) {
                P2 = material.get("optical_momentum_matrix_element_P2");                
        } else {
                P2 = kp_matrix.get_optical_momentum_matrix_element_P2();
        }
	// -------------------------------------------------
	// i'm not sure if the WZ kp8x8 matrix is correctly derived
	// therefore i use here chuangs kane matrix 8x8 (eq. 5 in chuang)
	// ------------------------------------------------- 
}
	
/** initialize momentum matrix tensor */
void MomentumOperatorBulk6x6WZ::MomentumOperatorBulk6x6WZ::lock() {
	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;
	tensor.assign(6,0.0);	
	
	cplx i(0.0, 1.0);
	double P2sqrt2 = P2 / constants::sqrt2;
	switch(get_operator_direction()) {
		// ----------------------------------
		// px
		// ----------------------------------
		case 0:
			tensor[0] = - P2sqrt2;
			tensor[1] =   P2sqrt2;
			tensor[3] =   P2sqrt2;
			tensor[4] = - P2sqrt2;
			break;
		// ----------------------------------
		// py
		// ----------------------------------			
		case 1:
			tensor[0] = - i * P2sqrt2;
			tensor[1] = - i * P2sqrt2;
			tensor[3] = - i * P2sqrt2;
			tensor[4] = - i * P2sqrt2;
			break;
		// ----------------------------------
		// pz
		// ----------------------------------			
		case 2:
			tensor[2] = P1;
			tensor[5] = P1;			
			break;
		default:
			TDKP_GENERAL_EXCEPTION("must not reach that point");			
	}

	for(unsigned int ii = 0; ii < 6; ii++) {
		tensor[ii] *= constants::m0 / constants::hbar;
	}	
}

/** test if supplied bandstructure is compatible */
bool MomentumOperatorBulk6x6WZ::compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const {
	if(cb_basis_size == 1 && cb_length == 1 && vb_basis_size == 6 && vb_length == 1) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << 1 << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << 6 << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << 1 << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << 1;
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;	
}

/** evaluate matrix element */
void MomentumOperatorBulk6x6WZ::eval(vector<cplx>& result, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	cplx tmp(0);	
	// spin up
	for(short ii = 0; ii < 3; ii++) {
		tmp += tensor[ii] * vb_wave[ii];
	}
	result[1] = tmp;
	tmp = 0.0e0;
	// spin down
	for(short ii = 3; ii < 6; ii++) {
		tmp += tensor[ii] * vb_wave[ii];
	}
	result[2] = tmp;
	result[0] = conj(result[1]) * result[1] + conj(result[2]) * result[2];
		
}




/** constructor for bulk 4x4 matrix element */
MomentumOperatorBulk4x4::MomentumOperatorBulk4x4(const double& optical_matrix_parameter)
: optical_matrix_param(optical_matrix_parameter)
{
	tensor.assign(8,0.0);
}
MomentumOperatorBulk4x4::MomentumOperatorBulk4x4(const Material& material)
: optical_matrix_param(material.get("optical_matrix_element"))
{
	tensor.assign(8,0.0);
}

/** lock mechanism for operator initializes it */
void MomentumOperatorBulk4x4::lock() {

	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;

	// --------------------------------------------------
	// we build 2 times a 4 vector as the basis is more complicated
	// as for the 6x6 enders/foreman matrix
	// --------------------------------------------------
	const double P = sqrt(constants::hbar_square / (2.0 * constants::m0) * optical_matrix_param);
	tensor.assign(8,0.0);
	cplx i(0.0, 1.0);

	// ---------------------------------------------------------
	// build tensor for spin up and down cb wave function
	//
	// how that was obtained:
	//   - the kp 4x4 matrix uses the angular momentum basis
	//   - the sebi transformation transforms |X>,|Y>,|Z> (+ spin)
	//     into that basis (equal to chuangs transformation)
	//   - the angular momentum basis is given by (4.2.22 - 4.2.23)
	//     in chuang, p. 135
	//   - let the angular momentum basis with total angular moment
	//     of 3/2 be u1 - u4
	//   - the momentum matrix element for the effmass cb is given by
	//     <S up/down | pi | uj> where we have to sum over the uj's
    // ----------------------------------------------------------
	switch(get_operator_direction()) {
		// --------------------------------------------------------
		// px operator
		// --------------------------------------------------------
		case 0:
			// |S up> spin up
			tensor[0] = - 1.0 / constants::sqrt2 * i * P;
			tensor[2] =   1.0 / constants::sqrt6 * i * P;
			// |S down> spin down
			tensor[5] = - 1.0 / constants::sqrt6 * i * P;
			tensor[7] =   1.0 / constants::sqrt2 * i * P;
			break;
		// --------------------------------------------------------
		// py operator
		// --------------------------------------------------------
		case 1:
			// |S up>
			tensor[0] =   1.0 / constants::sqrt2 * P;
			tensor[2] =   1.0 / constants::sqrt6 * P;
			// |S  down>
			tensor[5] =   1.0 / constants::sqrt6 * P;
			tensor[7] =   1.0 / constants::sqrt2 * P;
			break;
		// --------------------------------------------------------
		// pz operator
		// --------------------------------------------------------
		case 2:
			// |S up>
			tensor[1] = constants::sqrt2 / constants::sqrt3 * i * P;
			// |S down>
			tensor[6] = constants::sqrt2 / constants::sqrt3 * i * P;
			break;
		default:
			TDKP_GENERAL_EXCEPTION("strange unhandled momentum direction");
	}


}

bool  MomentumOperatorBulk4x4::compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const {
	if(cb_basis_size == 1 && cb_length == 1 && vb_basis_size == 4 && vb_length == 1) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << 1 << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << 4 << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << 1 << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << 1;
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;
}
void MomentumOperatorBulk4x4::eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	TDKP_ASSERT(cb_wave.size() == 1 && vb_wave.size() == 4, "cb_wave.size() == 1 && vb_wave.size() == 4");

	cplx cb_spin_up(0.0);
	cplx cb_spin_down(0.0);

	for(short ii = 0; ii < 4; ii++) {
		cb_spin_up   += tensor[ii] * vb_wave[ii];
		cb_spin_down += tensor[ii + 4] * vb_wave[ii];
	}
	// <S up   | pi | V>
	cb_spin_up   = constants::m0 / constants::hbar * cb_spin_up;
	// <S down | pi | V>
	cb_spin_down = constants::m0 / constants::hbar * cb_spin_down;
	// square and sum
	results[0] = (cb_spin_up.real()   * cb_spin_up.real()   + cb_spin_up.imag() * cb_spin_up.imag())
	           + (cb_spin_down.real() * cb_spin_down.real() + cb_spin_down.imag() * cb_spin_down.imag());

	results[1] = cb_spin_up;
	results[2] = cb_spin_down;	           

}


/** contstructor for simple effective mass matrix element calculation
 * 
 * |<c|p|v>|^2 = (m/hbar)^2 P^2 = (m/hbar)^2 * hbar^2 / (2m) Ep
 *             = m/2 Ep
 * 
 * and
 * 
 * |<c|mu|v>|^2 = ec^2 / (m^2 omega^2) |<c|p|v>|^2
 *  
 */
MomentumOperatorBulkEffectiveMass::MomentumOperatorBulkEffectiveMass(
	const double& optical_matrix_parameter_
) 
: optical_matrix_param(optical_matrix_parameter_)
{
}
void MomentumOperatorBulkEffectiveMass::lock() {
	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;
}
bool MomentumOperatorBulkEffectiveMass::compatible(
	short cb_basis_size, int cb_length, short vb_basis_size, int vb_length
) const {
	if(cb_basis_size == 1 && cb_length == 1 && vb_basis_size == 1 && vb_length == 1) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << 1 << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << 1 << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << 1 << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << 1;
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;	
}

void MomentumOperatorBulkEffectiveMass::eval(
	vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave
) const {
	TDKP_ASSERT(cb_wave.size() == 1 && vb_wave.size() == 1, "");
	// the * 2 is used to include the spin summation!
	results[0] = constants::m0 / 2.0 * optical_matrix_param * 2.0; 
	results[1] = sqrt(results[0]);
}


// --------------------------------------------------------------------
// base quantized momentum operator
// --------------------------------------------------------------------
MomentumOperatorQuantized::MomentumOperatorQuantized(int num_equation_per_node_, short cb_basis_size, short vb_basis_size, Geometry& geometry_, MaterialDatabase& material_database_)
: LinearProblem<complex<double> >(geometry_, material_database_),  
  fem_assembler(
  	geometry_, *this,
  	new NoSolver<cplx>(num_equation_per_node_ * geometry_.get_num_nodes(), nonsymmetric_matrix)
  ), 
  locked(false),
  obj_cb_basis_size(cb_basis_size),
  obj_vb_basis_size(vb_basis_size)
{

	num_equations_per_node = num_equation_per_node_;

	// ------------------------------------------------
	// initialize geometry, material and problem class
	// ------------------------------------------------
	geometry_.set_boundary_conditions(new BCIncludeAll(geometry_));
	//TDKP_ASSERT(this->geometry.verify(), "this->geometry.verify()");
	
	// ------------------------------------------------
	// TODO: add strain dependence
	// ------------------------------------------------ 
	

	// -------------------------------------------------
	// remaining variables
	// -------------------------------------------------
	k_values[0] = k_values[1] = k_values[2] = 0.0;

}

MomentumOperatorQuantized::~MomentumOperatorQuantized() {

}

bool MomentumOperatorQuantized::compatible(short cb_basis_size, int cb_length, short vb_basis_size, int vb_length) const {
	if(cb_basis_size == obj_cb_basis_size
	   && cb_length == (signed)this->geometry.get_num_nodes()
	   && vb_basis_size == obj_vb_basis_size
	   && vb_length == (signed)this->geometry.get_num_nodes()) {
		return true;
	}
	ostringstream sout;
	sout << "cb_basis_size (you) = " << cb_basis_size << ", required = " << obj_cb_basis_size << "\n"
	 	 << "vb_basis_size (you) = " << vb_basis_size << ", required = " << obj_vb_basis_size << "\n"
	 	 << "cb_length     (you) = " << cb_length     << ", required = " << obj_cb_basis_size * this->geometry.get_num_nodes() << "\n"
	 	 << "vb_length     (you) = " << vb_length     << ", required = " << obj_vb_basis_size * this->geometry.get_num_nodes();
	Logger::get_instance()->emit(LOG_ERROR, sout.str());
	return false;
}

void MomentumOperatorQuantized::lock() {
	TDKP_ASSERT(!locked, "momentum operator already locked");
	locked = true;
	// ---------------------------------------------------------
	// assemble system
	// ---------------------------------------------------------
	fem_assembler.reset_matrix_to_zero();
	fem_assembler.assemble_system();
}

void MomentumOperatorQuantized::release() {
	TDKP_ASSERT(locked, "you can't release a non-locked operator");
	locked = false;
}

void MomentumOperatorQuantized::set_k_value(const DomainPoint& point) {
	TDKP_ASSERT(!locked, "momentum operator is locked");
	k_values[0] = point.get_coord(0);
	k_values[1] = point.get_coord(1);
	k_values[2] = point.get_coord(2);
}





// --------------------------------------------------------
// simple base momentum operator for quantized stuff
// --------------------------------------------------------

MomentumOperatorSimpleXxX::MomentumOperatorSimpleXxX(int num_equation_per_node_, short cb_basis_size_, short vb_basis_size_,Geometry& geometry_, MaterialDatabase& material_database_)
:  MomentumOperatorQuantized(num_equation_per_node_, cb_basis_size_, vb_basis_size_, geometry_, material_database_)
{
}
void MomentumOperatorSimpleXxX::calculate_element_matrices(const Element* elem, cplx* lhs, cplx *rhs, int* internal_idx, int &n) const {

	double lmass[Element::max_num_nodes][Element::max_num_nodes]; 	   /* mass matrix */
	int    nnode;              /* number of nodes */
	int    neq;                /* number of kp equations */
	int    lsize; 			   /* size of lhs, rhs */
	int    nsparse;            /* number of sparse interaction matrix elements */

	const vector<cplx>& momentum_tensor = this->momentum_tensors[elem->get_region().get_material().get_id()];
	n     = 0;
	nnode = (signed)elem->get_num_nodes();
	neq   = this->get_num_equations_per_node();
	lsize = nnode * nnode * neq * neq;

	this->get_node_sparsity_pattern(nsparse);

	TDKP_ASSERT(nnode <= Element::max_num_nodes, "nnode <= Element::max_num_nodes (working arrays to small ...)");
	TDKP_ASSERT(nsparse == (signed)momentum_tensor.size(), "nsparse == momentum_tensor.size()");

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
			lmass[ii][jj] = elem->get_element_integral_0th_order(internal_idx[ii],internal_idx[jj]);
		}
	}

	// ------------------------------------------------------------------------------
	// assemble ....
	// ------------------------------------------------------------------------------
	int offset = 0;
	// for all nonzero nodes in element
	for(int ii = 0; ii < n; ii++) {
		for(int jj = 0; jj < n; jj++) {
			offset = (ii * n + jj) * nsparse;
			// build diagonal tensor
			for(int ss = 0; ss < nsparse; ss++) {
				lhs[offset + ss] += lmass[ii][jj] * momentum_tensor[ss];
			}
		}
	}
}


// --------------------------------------------------------
// momentum operator for quantized stuff and kp 6x6
// --------------------------------------------------------
/** the sparsity pattern is diagonal */
const int MomentumOperator6x6::sparsity_pattern[12] = {
	0,0, 1,1, 2,2, 3,3, 4,4, 5,5
};


MomentumOperator6x6::MomentumOperator6x6(Geometry& geometry_, MaterialDatabase& material_database_)
: MomentumOperatorSimpleXxX(6, 1, 6, geometry_, material_database_)
{		
	// -----------------------------------------------
	// interesting! i can't do the next two commands in
	// the base class constructor, because then
	// the passed this has the vtable of the base class and
	// the code terminates when get_node_sparsity_pattern
	// is called (because its virtual and defined in the derived
	// momentum operator)
	// -----------------------------------------------
	this->fem_assembler.create_matrix_structures();
}
MomentumOperator6x6::~MomentumOperator6x6() { }

const int* MomentumOperator6x6::get_node_sparsity_pattern(int& n) const {
	n = 6;
	return sparsity_pattern;
}

/** prepares our object with necessary values */

void MomentumOperator6x6::prepare() {

	cplx i(0.0, 1.0);

	// ---------------------------------------------------------
	// the tensor for quantized stuff is bigger and includes the
	// integrals of element shape functions Ni Nj
	// but the principle is the same as for the bulk case
	// ---------------------------------------------------------
	this->momentum_tensors.resize(this->material_db.get_num_materials());
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		double param = sqrt(constants::m0 / 2.0 * this->material_db.get_material(ii)->get("optical_matrix_element"));
		this->momentum_tensors[ii].assign(6, 0.0);
		this->momentum_tensors[ii][get_operator_direction()]     = i * param;
		this->momentum_tensors[ii][get_operator_direction() + 3] = i * param;
	}

	this->ready = true;

}


void MomentumOperator6x6::eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	TDKP_ASSERT(vb_wave.size() / cb_wave.size() == 6, "vb_wave.size() / cb_wave.size() == 6");
	TDKP_ASSERT(cb_wave.size() == this->geometry.get_num_nodes(), "cb_wave.size() == this->geometry.get_num_nodes()");

	vector<cplx> tmp;
	// --------------------------------------------------
	// multiply with tensor matrix
	// --------------------------------------------------
	fem_assembler.multiply_with_lhs(vb_wave, tmp);

	// --------------------------------------------------
	// assemble spin up and spin down
	// --------------------------------------------------
	cplx spin_up(0), spin_down(0);
	cplx ctrl(0.0);
	for(int ii = 0; ii < (signed)cb_wave.size(); ii++) {
		for(int jj = 0; jj < 3; jj++) {
			spin_up   += conj(cb_wave[ii]) * tmp[ii * 6 + jj];
			spin_down += conj(cb_wave[ii]) * tmp[ii * 6 + 3 + jj];
		}
	}

	// --------------------------------------------------
	// square, abs, sum and return
	// --------------------------------------------------
	results[0] = (spin_up.real()   * spin_up.real()   + spin_up.imag()   * spin_up.imag())
	           + (spin_down.real() * spin_down.real() + spin_down.imag() * spin_down.imag());
	results[1] = spin_up;
	results[2] = spin_down;

}


MomentumOperator6x6WZ::MomentumOperator6x6WZ(Geometry& geometry_, MaterialDatabase& material_database_)
: MomentumOperator6x6(geometry_, material_database_) {
	
}
void MomentumOperator6x6WZ::prepare() {

	cplx i(0.0, 1.0);

	// ---------------------------------------------------------
	// the tensor for quantized stuff is bigger and includes the
	// integrals of element shape functions Ni Nj
	// but the principle is the same as for the bulk case
	// ---------------------------------------------------------
	this->momentum_tensors.resize(this->material_db.get_num_materials());
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		KPMatrix8x8Wurtzite kp_matrix(this->material_db.get_material(ii));
		kp_matrix.calculate();
		double P1, P2;
        	if(this->material_db.get_material(ii)->is_set("optical_momentum_matrix_element_P1")) {
                	P1 = this->material_db.get_material(ii)->get("optical_momentum_matrix_element_P1");
         	} else {
                	P1 = kp_matrix.get_optical_momentum_matrix_element_P1();
        	}
        	if(this->material_db.get_material(ii)->is_set("optical_momentum_matrix_element_P2")) {
                	P2 = this->material_db.get_material(ii)->get("optical_momentum_matrix_element_P2");
        	} else {
                	P2 = kp_matrix.get_optical_momentum_matrix_element_P2();
        	}
		
		// add m0 / hbar		
		P2 *= constants::m0 / constants::hbar;
		P1 *= constants::m0 / constants::hbar;		
		double P2sqrt2 = P2 / constants::sqrt2;
		this->momentum_tensors[ii].assign(6, 0.0);
		TDKP_BOUNDS_ASSERT(!isnan(P1), "");
		TDKP_BOUNDS_ASSERT(!isnan(P2), "");
		switch(get_operator_direction()) {
			// ----------------------------------
			// px
			// ----------------------------------
			case 0:
				this->momentum_tensors[ii][0] = - P2sqrt2;
				this->momentum_tensors[ii][1] =   P2sqrt2;
				this->momentum_tensors[ii][3] =   P2sqrt2;
				this->momentum_tensors[ii][4] = - P2sqrt2;
				break;
			// ----------------------------------
			// py
			// ----------------------------------			
			case 1:
				this->momentum_tensors[ii][0] = - i * P2sqrt2;
				this->momentum_tensors[ii][1] = - i * P2sqrt2;
				this->momentum_tensors[ii][3] = - i * P2sqrt2;
				this->momentum_tensors[ii][4] = - i * P2sqrt2;
				break;
			// ----------------------------------
			// pz
			// ----------------------------------			
			case 2:
				this->momentum_tensors[ii][2] = P1;
				this->momentum_tensors[ii][5] = P1;			
				break;
			default:
				TDKP_GENERAL_EXCEPTION("must not reach that point");			
		}		
	}

	this->ready = true;

	
}



// --------------------------------------------------------
// momentum operator for quantized kp 4x4
// --------------------------------------------------------
/** the sparsity pattern is diagonal, and we assemble both, spin up and spin down |S> */
const int MomentumOperator4x4::sparsity_pattern[16] = {
	0,0, 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7
};


MomentumOperator4x4::MomentumOperator4x4(Geometry& geometry_, MaterialDatabase& material_database_)
: MomentumOperatorSimpleXxX(8, 1, 4, geometry_, material_database_)
{
	num_equations_per_node = 8;
	// -----------------------------------------------
	// interesting! i can't do the next two commands in
	// the base class constructor, because then
	// the passed this has the vtable of the base class and
	// the code terminates when get_node_sparsity_pattern
	// is called (because its virtual and defined in the derived
	// momentum operator)
	// -----------------------------------------------
	this->fem_assembler.create_matrix_structures();
}
MomentumOperator4x4::~MomentumOperator4x4() { }

const int* MomentumOperator4x4::get_node_sparsity_pattern(int& n) const {
	n = 8;
	return sparsity_pattern;
}

/** prepares our object with necessary values */
void MomentumOperator4x4::prepare() {

	cplx i(0.0, 1.0);

	// ---------------------------------------------------------
	// we build the kp 4x4 tensor here as an 1x8 tensor (or better:
	// (8x8 matrix with only diagonal elements), where the first
	// four components are spin up, the other components are spin down
	// ---------------------------------------------------------
	this->momentum_tensors.resize(this->material_db.get_num_materials());
	for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
		double param = sqrt(constants::m0 / 2.0 * this->material_db.get_material(ii)->get("optical_matrix_element"));
		this->momentum_tensors[ii].assign(8, 0.0);
		switch(get_operator_direction()) {
			// ---------------------------------------------------------
			// px
			// ---------------------------------------------------------
			case 0:
				// spin up
				this->momentum_tensors[ii][0] = - 1.0 / constants::sqrt2 * i * param;
				this->momentum_tensors[ii][2] =   1.0 / constants::sqrt6 * i * param;
				// spin down
				this->momentum_tensors[ii][5] = - 1.0 / constants::sqrt6 * i * param;
				this->momentum_tensors[ii][7] =   1.0 / constants::sqrt2 * i * param;
				break;
			// ---------------------------------------------------------
			// py
			// ---------------------------------------------------------
			case 1:
				// spin up
				this->momentum_tensors[ii][0] =   1.0 / constants::sqrt2 * i * param;
				this->momentum_tensors[ii][2] =   1.0 / constants::sqrt6 * i * param;
				// spin down
				this->momentum_tensors[ii][5] =   1.0 / constants::sqrt6 * i * param;
				this->momentum_tensors[ii][7] =   1.0 / constants::sqrt2 * i * param;
				break;
			// ---------------------------------------------------------
			// pz
			// ---------------------------------------------------------
			case 2:
				// spin up
				this->momentum_tensors[ii][1] =   constants::sqrt2 / constants::sqrt3 * i * param;
				// spin down
				this->momentum_tensors[ii][6] =   constants::sqrt2 / constants::sqrt3 * i * param;
				break;
			default:
				TDKP_GENERAL_EXCEPTION("invalid direction");
		}
	}
	this->ready = true;
}


void MomentumOperator4x4::eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	TDKP_ASSERT(vb_wave.size() / cb_wave.size() == 4, "vb_wave.size() / cb_wave.size() == 4");
	TDKP_ASSERT(cb_wave.size() == this->geometry.get_num_nodes(), "cb_wave.size() == this->geometry.get_num_nodes()");

	// --------------------------------------------------
	// the tensor we defined requires me to reorder
	// the valence band wavefunction so that the bands get
	// repeated. that means, the envelopes are expanded from
	// [0,3] to [0,7] where [4,7] is the copy of [0,3]
	// --------------------------------------------------
	vector<cplx> vb_wave_expanded(vb_wave.size() * 2);

	for(int ii = 0; ii < (signed)this->geometry.get_num_nodes(); ii++) {
		for(unsigned int jj = 0; jj < 4; jj++) {
			vb_wave_expanded[ii * 8 + jj]     = vb_wave[ii * 4 + jj];
			vb_wave_expanded[ii * 8 + jj + 4] = vb_wave[ii * 4 + jj];
		}
	}


	vector<cplx> tmp;
	// --------------------------------------------------
	// multiply with tensor matrix
	// --------------------------------------------------
	fem_assembler.multiply_with_lhs(vb_wave_expanded, tmp);

	// --------------------------------------------------
	// assemble spin up and spin down
	// --------------------------------------------------
	cplx spin_up(0), spin_down(0);
	for(int ii = 0; ii < (signed)cb_wave.size(); ii++) {
		for(int jj = 0; jj < 4; jj++) {
			spin_up   += conj(cb_wave[ii]) * tmp[ii * 8 + jj];
			spin_down += conj(cb_wave[ii]) * tmp[ii * 8 + 4 + jj];
		}
	}

	// --------------------------------------------------
	// square, abs, sum and return
	// --------------------------------------------------
	results[0] = (spin_up.real()   * spin_up.real()   + spin_up.imag()   * spin_up.imag())
	           + (spin_down.real() * spin_down.real() + spin_down.imag() * spin_down.imag());
	results[1] = spin_up;
	results[2] = spin_down;

}





// -------------------------------------------------------------
// simple kp 8x8 matrix element calculation
// -------------------------------------------------------------
MomentumOperator8x8Base::MomentumOperator8x8Base(Geometry& geometry_, MaterialDatabase& material_database_)
: MomentumOperatorSimpleXxX(8, 8, 8, geometry_, material_database_)
{

}

void MomentumOperator8x8Base::create_matrix_structures() {
	// -----------------------------------------------
	// build sparsity pattern from base matrix
	// -----------------------------------------------
	KPMatrixBase* tmp = this->get_kp_matrix();
	int n;
	const int* sparse_pat = tmp->get_sparsity_pattern(n);
	sparsity_pattern.assign(2*n, 0);
	for(int ii = 0; ii < 2*n; ii++) {
		sparsity_pattern[ii] = sparse_pat[ii];
	}
	delete tmp;
	this->fem_assembler.create_matrix_structures();	
}

MomentumOperator8x8Base::~MomentumOperator8x8Base() {
	for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
		delete this->kp_matrices[ii];
	}
}

const int* MomentumOperator8x8Base::get_node_sparsity_pattern(int& n) const {
	n = sparsity_pattern.size() / 2;
	return &sparsity_pattern[0];
}

/** prepares our object with necessary values */
void MomentumOperator8x8Base::prepare() {

	// ---------------------------------------------------------
	// the kp 8x8 tensor is built from the first order matrices
	// but be careful, this is a wrong definition of the momentum tensor
	// it the momentum tensor is given at k = 0. at k != 0,
	// we have remote states that mix into the gamma states,
	// so the momentum tensor would also need to incorporate these
	// ---------------------------------------------------------
	if((signed)this->kp_matrices.size() != this->material_db.get_num_materials()) {
		// delete (old) kp matrices
		for(unsigned int ii = 0; ii < this->kp_matrices.size(); ii++) {
			delete this->kp_matrices[ii];
		}
		this->kp_matrices.resize(this->material_db.get_num_materials());
		for(int ii = 0; ii < this->material_db.get_num_materials(); ii++) {
			this->kp_matrices[ii] = this->get_kp_matrix(this->material_db.get_material(ii));
			this->kp_matrices[ii]->calculate();
		}
	}

	// --------------------------------------------------------
	// now, already collect all the information into the momentum tensors
	// --------------------------------------------------------
	StrainTensor zero;
	cplx i(0.0, 1.0);
	vector<cplx> left, right;
	int length;
	this->get_node_sparsity_pattern(length);
	this->momentum_tensors.resize(this->material_db.get_num_materials());
	for(int mm = 0; mm < this->material_db.get_num_materials(); mm++) {
		this->momentum_tensors[mm].assign(length, 0.0e0);
		this->kp_matrices[mm]->build_first_order(KPMatrixBase::op_left,  zero, get_operator_direction(), left);
		this->kp_matrices[mm]->build_first_order(KPMatrixBase::op_right, zero, get_operator_direction(), right);
		for(short ss = 0; ss < length; ss++) {
			this->momentum_tensors[mm][ss] = constants::m0 / constants::hbar * (left[ss] + right[ss]) * i;
		}
	}

	this->ready = true;
}


void MomentumOperator8x8Base::eval(vector<cplx>& results, const vector<cplx>& cb_wave, const vector<cplx>& vb_wave) const {

	TDKP_ASSERT(vb_wave.size() == cb_wave.size(), "vb_wave.size() == cb_wave.size()");
	TDKP_ASSERT(cb_wave.size() / 8 == this->geometry.get_num_nodes(), "cb_wave.size() / 8 == this->geometry.get_num_nodes()");

	vector<cplx> tmp;
	// --------------------------------------------------
	// multiply with tensor matrix
	// --------------------------------------------------
	fem_assembler.multiply_with_lhs(vb_wave, tmp);

	TDKP_ASSERT(tmp.size() == cb_wave.size(), "tmp.size() == cb_wave.size()");

	// --------------------------------------------------
	// assemble
	// --------------------------------------------------
	cplx ret(0.0);
	for(int ii = 0; ii < (signed)cb_wave.size(); ii++) {
		ret += conj(cb_wave[ii]) * tmp[ii];
	}

	// --------------------------------------------------
	// square, abs, sum and return
	// --------------------------------------------------
	results[0] = ret.real() * ret.real() + ret.imag() * ret.imag();
	results[1] = ret;	
}
MomentumOperator8x8::MomentumOperator8x8(Geometry& geometry_, MaterialDatabase& material_database_) 
: MomentumOperator8x8Base(geometry_, material_database_)
{
	this->create_matrix_structures();
}
KPMatrixBase* MomentumOperator8x8::get_kp_matrix() const {
	return new KPMatrix8x8EndersForeman();	
}
KPMatrixBase* MomentumOperator8x8::get_kp_matrix(const Material* material) const {
	return new KPMatrix8x8EndersForeman(material);	
}
MomentumOperator8x8WZ::MomentumOperator8x8WZ(Geometry& geometry_, MaterialDatabase& material_database_) 
: MomentumOperator8x8Base(geometry_, material_database_)
{
	this->create_matrix_structures();
}
KPMatrixBase* MomentumOperator8x8WZ::get_kp_matrix() const {
	return new KPMatrix8x8Wurtzite();	
}
KPMatrixBase* MomentumOperator8x8WZ::get_kp_matrix(const Material* material) const {
	return new KPMatrix8x8Wurtzite(material);	
}

// --------------------------------------------------------
// Quantized Momentum Operator Factory
// --------------------------------------------------------
template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<KP8x81D2D,KP8x81D2D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator8x8(geometry_, material_database_);
}
template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<KP8x81D2DWZ,KP8x81D2DWZ>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator8x8WZ(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP6x61D2D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator6x6(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP6x61D2DWZ>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator6x6WZ(geometry_, material_database_);
}


template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP4x41D2D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator4x4(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, EffectiveMass>(Geometry& geometry_, MaterialDatabase& material_database_) {
	TDKP_GENERAL_EXCEPTION("implement me ... ");
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<KP8x83D, KP8x83D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator8x8(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<KP8x83DWZ, KP8x83DWZ>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator8x8WZ(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP6x63D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator6x6(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP6x63DWZ>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator6x6WZ(geometry_, material_database_);
}

template<>
MomentumOperatorQuantized* MomentumOperatorQuantized::factory<EffectiveMass, KP4x43D>(Geometry& geometry_, MaterialDatabase& material_database_) {
	return new MomentumOperator4x4(geometry_, material_database_);
}

} // end of namespace
