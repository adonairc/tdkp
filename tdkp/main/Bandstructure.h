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

#ifndef BANDSTRUCTURE_H_
#define BANDSTRUCTURE_H_

#include "tdkp/common/all.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/Domain.h"
#include "tdkp/main/EigenSolution.h"

#include <vector>

namespace tdkp {


/** bandstructure container for dispersive bandstructures */
template<class T>
class Bandstructure  {
public:
	Bandstructure(const char* filename);
	Bandstructure(int basis_size, int num_subbands, int solution_length, int num_k_values_);
	Bandstructure(const Bandstructure& copy);
	virtual ~Bandstructure();

	bool operator==(const Bandstructure<T>& rhs) const;
	const Bandstructure& operator=(const Bandstructure<T>& rhs) { TDKP_GENERAL_EXCEPTION("assignment operator of base class is undefined!"); }

	// ------------------------------------------------------
	// solution access
	// ------------------------------------------------------			       
	const EigenSolution<T>& get_eigensolution(int kidx, int band_idx) const throw(Exception*);
	EigenSolution<T>&       get_eigensolution(int kidx, int band_idx) throw(Exception*);
	const EigenSolution<T>* get_eigensolution_ptr(int kidx, int band_idx) const throw(Exception*);	
	EigenSolution<double>*  get_eigensolution_real(int kidx, int band_idx) const throw(Exception*);
	EigenSolution<double>*  get_probability(int kidx, int band_idx) const throw(Exception*);
	const T& get_energy(int kidx, int band_idx) const;
	T&       get_energy(int kidx, int band_idx);
	void     set_energy(int kidx, int band_idx, const T& value);
		 
	// ------------------------------------------------------
	// bandstructure property access
	// ------------------------------------------------------	
	/** return the number of subbands */
	int get_number_of_bands()    const throw(Exception*) { return this->num_bands; }
	/** return the bloch basis size */
	int get_basis_size()         const throw(Exception*) { return this->basis_size; }
	/** return the maximum kidx we have */
	int get_number_of_k_values() const throw(Exception*) { return this->num_k_values; }
	/** return the length of a solution (should be = #nodes * #basis_size) */
	int get_solution_length()    const throw(Exception*) { return this->solution_length; }

	// -------------------------------------------------------
	// io routines, non-portable and ugly routines ...
	// -------------------------------------------------------
	void write_binary(const char* filename) const throw(Exception*) ;
	void read_binary(const char*  filename) throw(Exception*);

	// ------------------------------------------------------
	// resorting of bands, swapping routine
	// ------------------------------------------------------
	void swap(unsigned int kidx1, unsigned int band_idx1, unsigned int kidx2, unsigned int band_idx2);

	// -------------------------------------------------------
	// bandstructure assembly routines, used by problem class
	// -------------------------------------------------------	
	void add_eigensolution(int kidx, int band_idx, EigenSolution<T>* solution) throw(Exception*);
			
protected:
	/** binary file writer: all derived classes implement this and call their parents function upon execution */
	virtual void write_binary_to_stream(ostream& out) const;
	/** binary file reader: all derived classes implement this and call their parents function upon execution */
	virtual void read_binary_from_stream(istream& in);	
	
	Bandstructure();

protected:
	vector<EigenSolution<T>*> eigensolutions;
	vector<T> dispersion_relations;
	void init(int basis_size, int num_bands, int solution_length, int num_k_values);	
				
private:

	int       num_bands;
	int       basis_size;
	int       num_k_values;
	int       solution_length;
		
};

/** bandstructure object with values within a domain */
template<class T>
class BandstructureDomain : public Bandstructure<T>, public XYData<T> {
public:
	BandstructureDomain() {}
	BandstructureDomain(const char* filename);
	BandstructureDomain(int basis_size, int num_bands, int solution_length, const DomainMaster& domain);
	explicit BandstructureDomain(const BandstructureDomain& copy);
	const BandstructureDomain& operator=(const BandstructureDomain<T>& rhs);	
	virtual ~BandstructureDomain();	
	const DomainMaster& get_domain() const { return domain; }

	BandstructureDomain<T>* extract_bands(unsigned int num_bands, unsigned int offset = 0) const;
		
	// ----------------------------------------------
	// interface XYData implementation (for compatibility reasons)
	// ----------------------------------------------
	int    get_x_length()                const { return this->get_number_of_k_values(); }
	int    get_num_x_sets()              const;
	int    get_num_y_sets()              const { return this->get_number_of_bands();    }
	void   get_x(int xidx, vector<T> &x) const;
	void   get_y(int yidx, vector<T>& y) const;
	string get_x_identifier(int xidx)    const;
	string get_y_identifier(int yidx) 	 const;
		
	void reinit(int basis_size, int num_bands, int solution_length, const DomainMaster& domain_);
	
protected:
	/** binary file writer: all derived classes implement this and call their parents function upon execution */
	virtual void write_binary_to_stream(ostream& out) const;
	/** binary file reader: all derived classes implement this and call their parents function upon execution */
	virtual void read_binary_from_stream(istream& in);
	
private:
	DomainMaster domain;
};

/** static pml bandstructure filter function */
void extract_stable_bands( 
	const double& maximum_imag_abs, 
	const BandstructureDomain<complex<double> >& source, 
	BandstructureDomain<complex<double> >& target
);




/** initialize bandstructure from file */
template<class T>
Bandstructure<T>::Bandstructure(const char* filename)
: eigensolutions(0),
  dispersion_relations(0),
  num_bands(0),
  basis_size(0),
  num_k_values(0),
  solution_length(0)
{
	this->read_binary(filename);
}

/** protected constructor 
 * 
 * initializing a bandstructure from a file needs call to read_binary 
 * from derived constructor, after the object was constructed 
 */
template<class T>
Bandstructure<T>::Bandstructure() 
: eigensolutions(0),
  dispersion_relations(0),
  num_bands(0),
  basis_size(0),
  num_k_values(0),
  solution_length(0)
{
	
}

/** copy constructor */
template<class T> 
Bandstructure<T>::Bandstructure(const Bandstructure<T>& copy)
: eigensolutions(0),
  dispersion_relations(0),
  num_bands(copy.num_bands),
  basis_size(copy.basis_size),
  num_k_values(copy.num_k_values),
  solution_length(copy.solution_length)
{
	// copy eigenvalues and eigensolutions
	dispersion_relations = copy.dispersion_relations;
	this->eigensolutions.resize(copy.eigensolutions.size());
	// create eigensolutions
	for(int ii = 0; ii < (signed)copy.eigensolutions.size(); ii++) {		
		this->eigensolutions[ii] = new EigenSolution<T>(*copy.eigensolutions[ii]);
	}			
}

template<class T>
Bandstructure<T>::Bandstructure(int basis_size_, int num_bands_, int solution_length_, int num_k_values_) {
	this->init(basis_size_, num_bands_, solution_length_, num_k_values_);
}

template<class T>
Bandstructure<T>::~Bandstructure() {
	for(typename vector<EigenSolution<T>*>::iterator it = this->eigensolutions.begin(); it != this->eigensolutions.end(); it++) {
		if((*it) != 0) {
			delete (*it);
		}
	}
	this->eigensolutions.resize(0);
}

/** comparison operators 
 * 
 * */
template<class T> 
bool Bandstructure<T>::operator==(const Bandstructure<T>& rhs) const {
	
	if(num_bands != rhs.num_bands) { return false; }
	if(basis_size != rhs.basis_size) { return false; }
	if(num_k_values != rhs.num_k_values) { return false; }
	if(solution_length != rhs.solution_length) { return false; }

	if(dispersion_relations.size() != rhs.dispersion_relations.size()) { return false; }
	if(eigensolutions.size() != rhs.eigensolutions.size()) { return false; }
	// compare dispersion relation
	for(unsigned int ii = 0; ii < dispersion_relations.size(); ii++) {
		if(!tdkp_math::compare(dispersion_relations[ii],dispersion_relations[ii],1.0e-11)) {
			return false;	
		}	
	}
	// compare eigensolutions
	for(unsigned int ii = 0; ii < eigensolutions.size(); ii++) {
		if(!eigensolutions[ii]->compare_probability(*rhs.eigensolutions[ii])) {
			return false;	
		}	
	}	
	return true;
	
}

/** (re-)initialize bandstructure class
 *
 * @param basis_size_      size of kp basis (e.g. Chuang (4x4 KP) == 4, 8x8 kp  == 8 etc)
 * @param num_bands_       number of (sub) bands calculated
 * @param solution_length_ length of one envelope vector
 * @param num_k_values_    the number of k values
 */
template<class T>
void Bandstructure<T>::init(int basis_size_, int num_bands_, int solution_length_, int num_k_values_) {

	this->num_k_values    = num_k_values_;
	this->num_bands       = num_bands_;
	this->basis_size      = basis_size_;
	this->solution_length = solution_length_;
	this->dispersion_relations.resize(this->num_k_values * this->num_bands);
	this->eigensolutions.resize(this->num_k_values * this->num_bands); 
	for(int ii = 0; ii < this->num_k_values * this->num_bands; ii++) {
		this->dispersion_relations[ii] = 0;
		this->eigensolutions[ii]       = 0;		
	}
	
}

/** swap bands (kidx1,band_idx1) <-> (kidx2,band_idx2) */
template<class T>
void Bandstructure<T>::swap(unsigned int kidx1, unsigned int band_idx1, unsigned int kidx2, unsigned int band_idx2) {
	
	
	unsigned int idx1 = kidx1 + this->num_k_values * band_idx1;
	unsigned int idx2 = kidx2 + this->num_k_values * band_idx2;
	
	TDKP_ASSERT(idx1 < dispersion_relations.size(), "");
	TDKP_ASSERT(idx2 < dispersion_relations.size(), "");
	
	// ----------------------------------------
	// swap wavefunction
	// ----------------------------------------
	EigenSolution<T>* tmp = this->eigensolutions[idx1];
	this->eigensolutions[idx1] = this->eigensolutions[idx2];
	this->eigensolutions[idx2] = tmp;
	T tmp_value = dispersion_relations[idx1];
	dispersion_relations[idx1] = dispersion_relations[idx2];
	dispersion_relations[idx2] = tmp_value;
		
}

/** return energy(k) of band i */
template<class T>
const T& Bandstructure<T>::get_energy(int kidx, int band_idx) const {	
	TDKP_BOUNDS_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");
	return this->dispersion_relations[kidx + this->num_k_values * band_idx];
}

/** return energy(k) of band i */
template<class T>
T& Bandstructure<T>::get_energy(int kidx, int band_idx) {
	return const_cast<T&>(static_cast<const Bandstructure<T>&>(*this).get_energy(kidx, band_idx));	
}

/** set energy */
template<class T>
void Bandstructure<T>::set_energy(int kidx, int band_idx, const T& value) {
	this->get_energy(kidx,band_idx) = value;	
}		

/** add eigensolution to bandstructure */
template<class T>
void Bandstructure<T>::add_eigensolution(int kidx, int band_idx, EigenSolution<T>* solution) throw(Exception*) {	
	TDKP_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");
	TDKP_ASSERT(solution->get_length() == this->solution_length, "length of solution is not equal to announced length");
	TDKP_ASSERT(solution->get_num_data_per_node() == this->get_basis_size(), "solutions's basis size is not equal to announced basis size");
	TDKP_ASSERT(this->eigensolutions[kidx + this->num_k_values * band_idx] == 0, "there is already a solution stored! so this would leak huge amounts of memory!");
	this->eigensolutions[kidx + this->num_k_values * band_idx] = solution;
	this->dispersion_relations[kidx + this->num_k_values * band_idx] = solution->get_energy();
}

/** retuns const reference to eigensolution object */
template<class T>
const EigenSolution<T>& Bandstructure<T>::get_eigensolution(int kidx, int band_idx) const throw(Exception*) {
	TDKP_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");
	TDKP_ASSERT(kidx + this->num_k_values * band_idx < (signed)this->eigensolutions.size(), "kidx + this->num_k_values * band_idx < this->eigensolutions.size()");	 
	TDKP_ASSERT(this->eigensolutions[kidx + this->num_k_values * band_idx] != 0, "this->eigensolutions[kidx + this->num_k_values * band_idx] != 0");
	return *this->eigensolutions[kidx + this->num_k_values * band_idx];
}

template<class T>
EigenSolution<T>& Bandstructure<T>::get_eigensolution(int kidx, int band_idx) throw(Exception*) {
	return const_cast<EigenSolution<T>&>(
		static_cast<const Bandstructure<T>&>(*this).get_eigensolution(kidx,band_idx)
	);	
}

/** return pointer to eigensolution (needed for tcl ...) */
template<class T>
const EigenSolution<T>* Bandstructure<T>::get_eigensolution_ptr(int kidx, int band_idx) const throw(Exception*) {
	return &(this->get_eigensolution(kidx,band_idx)); 	
}


/** create a new eigen solution object of the probability distribution */
template<class T>
EigenSolution<double>*  Bandstructure<T>::get_probability(int kidx, int band_idx) const throw(Exception*) {
	TDKP_ASSERT(kidx >= 0 && kidx < this->num_k_values && band_idx >= 0 && band_idx < this->num_bands, "kidx or band idx out of range");		
	int length                   = this->get_solution_length();	
	const EigenSolution<T>& orig = this->get_eigensolution(kidx, band_idx); 
	vector<double> data(length); 

	double absval;
#pragma omp parallel for default(shared) private(absval)	
	for(int ii = 0; ii < length; ii++) {
		data[ii] = 0;
		for(int jj = 0; jj < this->basis_size; jj++) {
			absval = tdkp_math::abs(orig.get_node_value(ii, jj)); 
			data[ii] += absval * absval;		
		}				
	}	
	return new EigenSolution<double>(tdkp_math::only_real(orig.get_energy()), &data[0], length, 1);
}


/** write data to binary file */
template<class T>
void Bandstructure<T>::write_binary(const char* filename) const throw(Exception*) {
	
	ofstream fout(filename, ios::binary);
	const char* binary_file_header = "BandstructureDataBegin";
	const char* binary_file_footer = "BandstructureDataEnd";

	if(fout) {
		// -------------------------------------------------
		// write header data
		// -------------------------------------------------
		fout.write(binary_file_header, strlen(binary_file_header));		
				
		// -------------------------------------------------
		// call virtual functions
		// -------------------------------------------------
		write_binary_to_stream(fout);
		
		// -------------------------------------------------
		// write footer and close
		// -------------------------------------------------
		fout.write(binary_file_footer, strlen(binary_file_footer));
		fout.close();		

	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file");
	}
}

/** write class data to binary stream */
template<class T>
void Bandstructure<T>::write_binary_to_stream(ostream& fout) const {
		
	// ------------------------------------------------ 		
	// copy everything into struct
	// ------------------------------------------------
	struct serialdata {
		int num_bands; int basis_size; int num_k_values; int solution_length;	
	} serial;
	serial.num_k_values    = this->num_k_values;
	serial.num_bands       = this->num_bands;
	serial.basis_size      = this->basis_size;
	serial.solution_length = this->solution_length;
	
	// ------------------------------------------------
	// write struct and eigensolutions
	// ------------------------------------------------
	fout.write((char*)&serial, sizeof(serialdata));		
	for(int ii = 0; ii < this->num_k_values * this->num_bands; ii++) {
		TDKP_ASSERT(this->eigensolutions[ii] != 0, "this->eigensolutions[ii] != 0"); 
		this->eigensolutions[ii]->write_binary(fout);
	}
	
}

/** read data from binary file */
template<class T>
void Bandstructure<T>::read_binary(const char* filename) throw(Exception*) {
	
	const char* binary_file_header = "BandstructureDataBegin";
	const char* binary_file_footer = "BandstructureDataEnd";
		
	ifstream fin(filename, ios::binary);
	if(fin) {
		// -------------------------------------------
		// create buffer for header and footer test
		// -------------------------------------------
		int buflen = max(strlen(binary_file_header), strlen(binary_file_footer)) + 1;
		char* buf = new char[buflen];
		buf[strlen(binary_file_header)] = '\0';

		// -------------------------------------------
		// read header
		// -------------------------------------------
		fin.read(buf, strlen(binary_file_header));
		if(string(buf).compare(binary_file_header) != 0) {
			TDKP_GENERAL_EXCEPTION("header in bandstructure is invalid ... ");
		}		
		
		// -------------------------------------------
		// call specialized reader functions
		// -------------------------------------------
		read_binary_from_stream(fin);
		
		// -------------------------------------------
		// read footer and close file
		// -------------------------------------------
		buf[strlen(binary_file_footer)] = '\0';
		fin.read(buf, strlen(binary_file_footer));
		if(string(buf).compare(binary_file_footer) != 0) {
			ostringstream sout;
			sout << "the bandstructure you read from a file does not end with "
			     << "the expected footer. probably you saved it with a different "
			     << "bandstructure object. you should never do that. but as you did "
			     << "it, i guess you know what you are doing. so i won't quit now. "
			     << "if i segfault, it's your fault!";			      
			Logger::get_instance()->emit(LOG_WARN, sout.str());
		}
		fin.close();
				
	} else {
		TDKP_GENERAL_EXCEPTION("could not open file for reading");
	}
}

/** read class data from binary stream */
template<class T>
void Bandstructure<T>::read_binary_from_stream(istream& fin) {
	
	// ---------------------------------------------
	// delete any EigenSolution objects that exist at the moment
	// ---------------------------------------------
	if(this->eigensolutions.size() > 0) {
		for(unsigned int ii = 0; ii < this->eigensolutions.size(); ii++) {
			if(this->eigensolutions[ii] != 0) {
				delete this->eigensolutions[ii];
				this->eigensolutions[ii] = 0;
			}
		}
	}
	// ---------------------------------------------
	// define struct in order to read all values at one time
	// ---------------------------------------------
	struct serialdata {
		int num_bands; int basis_size; int num_k_values; int solution_length;
	} serial;

	// ---------------------------------------------
	// read data
	// ---------------------------------------------
	fin.read((char*)&serial, sizeof(serialdata));
	this->init(serial.basis_size, serial.num_bands, serial.solution_length, serial.num_k_values);

	// ---------------------------------------------
	// read eigensolutions
	// ---------------------------------------------
	for(int ii = 0; ii < this->num_k_values * this->num_bands; ii++) {
		EigenSolution<T>* tmp = new EigenSolution<T>();
		tmp->read_binary(fin);
		this->eigensolutions[ii] = tmp;
		this->dispersion_relations[ii] = tmp->get_energy();
	}
	
}


/** construct bandstructure domain from stored binary data */
template<class T>
BandstructureDomain<T>::BandstructureDomain(const char* filename) {
	this->read_binary(filename);	
}

/** basic bandstructure domain constructor */
template<class T>
BandstructureDomain<T>::BandstructureDomain(int basis_size_, int num_bands_, int solution_length_, const DomainMaster& domain_)
: Bandstructure<T>(basis_size_, num_bands_, solution_length_, domain_.get_number_of_points()),
  domain(domain_) 
{
	if(!domain.frozen()) {
		domain.freeze();
	}
}

template<class T>
BandstructureDomain<T>::BandstructureDomain(const BandstructureDomain& copy)
: Bandstructure<T>(copy),
  domain(copy.domain)
{
	
}

template<class T>
BandstructureDomain<T>::~BandstructureDomain() {
	
}	
	
/** write bandstructure to stream */	
template<class T>
void BandstructureDomain<T>::write_binary_to_stream(ostream& out) const {
	// ----------------------------------------------
	// first, we add the base class data 
	// ----------------------------------------------
	Bandstructure<T>::write_binary_to_stream(out);
	
	// ----------------------------------------------
	// and now we write the domain master
	// ----------------------------------------------
	domain.write_binary(out);	
}

/** read bandstructure from binary stream */
template<class T>
void BandstructureDomain<T>::read_binary_from_stream(istream& in) {
	// ----------------------------------------------
	// parents data was saved before mine ...
	// ----------------------------------------------
	Bandstructure<T>::read_binary_from_stream(in);
	
	// ----------------------------------------------
	// read domain master
	// ----------------------------------------------
	domain.read_binary(in);
	
}

template<class T>
int BandstructureDomain<T>::get_num_x_sets() const {

	if(domain.radial()) {
		return 1;	
	} else {
		return domain.get_dimension();	
	}
	
}

/** interface implementation for writing xy data files: returning length of vector between 0 and the k point as x-values */
template<class T>
void BandstructureDomain<T>::get_x(int xidx, vector<T> &x) const {
	
	TDKP_ASSERT((signed)domain.get_number_of_points() == this->get_number_of_k_values(), "domain.get_number_of_points() == this->get_number_of_k_values()");	
	x.resize(this->get_number_of_k_values());
	bool radial = domain.radial();
	TDKP_ASSERT(!(xidx > 0 && radial), "!(xidx > 0 && radial)");
	TDKP_ASSERT(xidx < 3 && xidx >= 0, "xidx < 3 && xidx >= 0");  
	for(int ii = 0; ii < this->get_number_of_k_values(); ii++) {
		if(radial) {
			x[ii] = domain.get_point(ii).get_coord_abs();
		} else {
			x[ii] = domain.get_point(ii).get_coord(xidx);	
		}		
	}
	
}

/** interface implementation for writing xy data files: returns y-values */
template<class T>
void BandstructureDomain<T>::get_y(int yidx, vector<T>& y) const {
	TDKP_ASSERT(yidx >= 0 && yidx < this->get_number_of_bands(), "invalid band requested");
	y.resize(this->get_number_of_k_values());
	for(int ii = 0; ii < this->get_number_of_k_values(); ii++) {
		y[ii] = this->get_energy(ii,yidx);
	}
}

/** interface implementation for writing xy data files: returns x-identifier */
template<class T>
string BandstructureDomain<T>::get_x_identifier(int xidx) const {
	ostringstream sout;
	sout << "k" << xidx;
	return sout.str();
}

/** interface implementation for writing xy data files: returns y-identifiers */
template<class T>
string BandstructureDomain<T>::get_y_identifier(int yidx) const {
	TDKP_ASSERT(yidx >= 0 && yidx < this->get_number_of_bands(), "invalid band requested");
	ostringstream sout;
	sout << "band" << yidx;
	return sout.str(); 
}
/** get part of the bands */
template<class T>
BandstructureDomain<T>* BandstructureDomain<T>::extract_bands(unsigned int num_bands, unsigned int offset) const {
	
	TDKP_ASSERT(offset + num_bands <= (unsigned)this->get_number_of_bands(), "too many bands requested!"); 	
	BandstructureDomain<T>* ret = new BandstructureDomain<T>(
		this->get_basis_size(),
		num_bands,
		this->get_solution_length(),
		this->get_domain()
	);
	// copy bands	
	for(int kk = 0; kk < this->get_number_of_k_values(); kk++) {
		for(unsigned int bb = offset; bb < offset + num_bands; bb++) {
			// special treatment for effmass (only copy energies if solution not available)
			if(this->eigensolutions[kk + this->get_number_of_k_values() * bb] != 0) {		
				EigenSolution<T>* sol =  new EigenSolution<T>(this->get_eigensolution(kk, bb));
				ret->add_eigensolution(kk,bb - offset,sol);
			} else {
				ret->get_energy(kk,bb - offset) = this->get_energy(kk,bb);	
			}				
		}
	}
	return ret;
	
}


/** reinits the bandstructure objects and deletes old data */
template<class T>
void BandstructureDomain<T>::reinit(int basis_size_, int num_bands_, int solution_length_, const DomainMaster& domain_) {

	// -----------------------------------------------
	// delete my data 
	// -----------------------------------------------								  
	for(unsigned int ii = 0; ii < this->eigensolutions.size(); ii++) {
		if(this->eigensolutions[ii] != 0) {
			delete this->eigensolutions[ii]; this->eigensolutions[ii] = 0;	
		}	
	}
	// -----------------------------------------------
	// reinit object 
	// -----------------------------------------------
	this->init(
		basis_size_,
		num_bands_,
		solution_length_,
		domain_.get_number_of_points()
	); 
	this->domain = domain_;
	if(!this->domain.frozen()) {
		domain.freeze();	
	}
		
}

/** assignment operator for bandstructure object */
template<class T>
const BandstructureDomain<T>& BandstructureDomain<T>::operator=(const BandstructureDomain<T>& rhs) {
	
	if(this == &rhs) {
		return *this;	
	}
	
	this->reinit(
		rhs.get_basis_size(), 
		rhs.get_number_of_bands(), 
		rhs.get_solution_length(),
		rhs.get_domain()		
	); 
	
	// -----------------------------------------------
	// copy data
	// -----------------------------------------------
	TDKP_ASSERT(rhs.dispersion_relations.size() == this->dispersion_relations.size(), "rhs.dispersion_relations.size() == dispersion_relations.size()");
	TDKP_ASSERT(rhs.eigensolutions.size() == this->eigensolutions.size(), "rhs.eigensolutions.size() == this->eigensolutions.size()");  
	for(unsigned int ii = 0; ii < rhs.dispersion_relations.size(); ii++) {
		this->dispersion_relations[ii] = rhs.dispersion_relations[ii];			
	} 
	
	for(unsigned int ii = 0; ii < rhs.eigensolutions.size(); ii++) {
		if(rhs.eigensolutions[ii] != 0) {
			this->eigensolutions[ii] = new EigenSolution<T>(*rhs.eigensolutions[ii]);	
		}	
	}
	
	return *this;
			
}

}

#endif /*BANDSTRUCTURE_H_*/
