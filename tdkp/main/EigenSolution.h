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

#ifndef EIGENSOLUTION_H_
#define EIGENSOLUTION_H_

#include <vector>
#include <fstream>
#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"

namespace tdkp {

template<class T>
class EigenSolution : public NodeData<T> {
	
public:
	typedef typename vector<T>::const_iterator data_const_iterator;
	typedef typename vector<T>::iterator       data_iterator;

	EigenSolution();
	EigenSolution(int num, int num_per_node = 1);
	EigenSolution(const EigenSolution<T>& copy);
	EigenSolution(T energy_, const T* data, int num, int num_per_node = 1);	
	~EigenSolution() { this->data.clear(); }
	void info() const;
	

	const EigenSolution<T>& add(const T& coeff, const EigenSolution<T>& add);
	
	// ------------------------------------------
	// overloaded operators
	// ------------------------------------------
	T&       operator()(int ii);	 

	const EigenSolution<T>& operator=(const EigenSolution<T>& rhs);

	bool     compare(const EigenSolution<T>& b, double tol = 1.0e-11, bool be_chatty = true) const;
	bool     compare_probability(const EigenSolution<T>& b, double tol = 1.0e-11, bool be_chatty = true) const;
	
	size_t   size() const { return this->data.size(); }
	T        get_energy() const { return this->energy; }	
	void     set_energy(T energy_) { this->energy = energy_; }
	void     read(const char* filename) throw(Exception*);
	void     write(const char* filename) const throw(Exception*);	
	void     read_binary(istream& fin) throw(Exception*);
	void     write_binary(ostream& fout) const;	
	unsigned int get_basis_size() const { return this->num_data_per_node; }
		
	// -------------------------------------------
	// iterator access
	// -------------------------------------------
	data_iterator       get_data_begin()       { return this->data.begin(); }
	data_const_iterator get_data_begin() const { return this->data.begin(); }
	data_const_iterator get_data_end()   const { return this->data.end();   }
	const vector<T>&    get_data() const       { return this->data; }
	
		
	// -------------------------------------------
	// node data functions
	// -------------------------------------------
	const T& get_node_value(unsigned int vidx, unsigned int eqidx = 0) const;
	T&       get_node_value(unsigned int vidx, unsigned int eqidx = 0);
	
	void     set_node_value(unsigned int vidx, unsigned int eqidx, T value) { this->data[vidx * this->num_data_per_node + eqidx] = value; }
	int      get_length() const { return this->data.size() / this->num_data_per_node; }
	void     set_length(int size, int num_data_per_node) { TDKP_GENERAL_EXCEPTION("must not be called"); 	}	
	string   get_identifier(unsigned int idx) const;
	void     set_identifier(unsigned int idx, const string& ident);
	

protected:
	/** read/write binary is implemented different. so calling that function is not allowed ... */
			
private:
	void reset();
	vector<T> data;			
	T energy;
	vector<string> identifier;
};


template<class T> 
class EigenSolutionSet : public NodeData<T>{
public:

	typedef typename vector<EigenSolution<T> * >::iterator       data_iterator;
	typedef typename vector<EigenSolution<T> * >::const_iterator data_const_iterator;

	EigenSolutionSet() {};
	~EigenSolutionSet();
	int 			  get_length() const;
	int               get_num_solutions() const;
	void   			  add(EigenSolution<T>* solution);	
	EigenSolution<T>& get(unsigned int idx);
	const EigenSolution<T>& get(unsigned int idx) const;
	void              write(const char* filename) const;
	void              read(const char* filename);
	void              info() const;
	
	// -------------------------------------------------
	// NodeData implementation
	// -------------------------------------------------
	const T& get_node_value(unsigned int vidx, unsigned int eqidx = 0) const;
	T&       get_node_value(unsigned int vidx, unsigned int eqidx = 0);
	void     set_node_value(unsigned int vidx, unsigned int eqidx, T value); 
	int      get_num_data_per_node() const;
	void     set_length(int size, int num_data_per_node_) { TDKP_GENERAL_EXCEPTION("the EigenSolutionSet can only read data with its own i/o functions"); 	}	
	string   get_identifier(unsigned int idx) const;
		
private:
	vector<EigenSolution<T>* > data;		
};


template<class T>
EigenSolution<T>::EigenSolution() 
: energy(0.0),
  identifier(1)
{
	this->num_data_per_node = 0;
}


/** creates eigensolution object from eigenvalue / eigenvector */
template<class T>
EigenSolution<T>::EigenSolution(T energy_, const T* vdata, int num, int num_per_node)
: data(num * num_per_node),
  energy(energy_),
  identifier(num_per_node)  
{
	this->num_data_per_node = num_per_node;
#ifndef DEADRAT
	#pragma omp parallel for default(shared)
#endif	
	for(int ii = 0; ii < (signed)this->data.size(); ii++) {
		data[ii] = vdata[ii];
	}
}

/** create empty eigensolution object */
template<class T> 
EigenSolution<T>::EigenSolution(int num, int num_per_node)
: data(num * num_per_node),
  energy(0.0),
  identifier(num_per_node)   
{
	this->num_data_per_node = num_per_node;
}

/** copy constructor to perform deep copy */
template<class T> 
EigenSolution<T>::EigenSolution(const EigenSolution<T>& copy)
: NodeData<T>(copy),
  data(copy.data.begin(), copy.data.end()),
  energy(copy.energy),
  identifier(copy.identifier.begin(), copy.identifier.end()) 
{
	
}

/** assignment operator */
template<class T> 
const EigenSolution<T>& EigenSolution<T>::operator=(const EigenSolution<T>& rhs) {
	
	if(this == &rhs) { 
		return *this;	
	}
	this->num_data_per_node = rhs.num_data_per_node;
 	data.assign(rhs.data.begin(), rhs.data.end());
  	energy = rhs.energy,
  	identifier.assign(rhs.identifier.begin(), rhs.identifier.end());
  	return *this;
}

/** add eigensolution |Y> + coeff |add> -> |R> 
 *
 * the energy is kept constant except if eigensolution object
 * was not initialized yet. then the energy of the first added
 * eigensolution is copied 
 */
template<class T>
const EigenSolution<T>& EigenSolution<T>::add(const T& coeff, const EigenSolution<T>& add) {
	
	// ----------------------------------------
	// empty eigensolution is |0> -> copy rhs and scale
	// ----------------------------------------
	if(data.size() == 0) {
		*this = add;
		for(unsigned int ii = 0; ii < data.size(); ii++) {
			data[ii] = coeff * data[ii];	
		}		
	} else {
		// ----------------------------------------
		// add rhs
		// ----------------------------------------
		TDKP_ASSERT(this->num_data_per_node == add.num_data_per_node,"");
		TDKP_ASSERT(this->data.size() == add.data.size(), "");
		for(unsigned int ii = 0; ii < data.size(); ii++) {
			data[ii] = data[ii] + coeff * add.data[ii];
		}
	}
	return *this;
}

/** raw access to solution vector */
template<class T>
T& EigenSolution<T>::operator()(int ii) {
	return this->data[ii];	
}

template<class T>
const T& EigenSolution<T>::get_node_value(unsigned int vidx, unsigned int eqidx) const {
	TDKP_BOUNDS_ASSERT(vidx >= 0  && vidx <  (this->data.size() / this->num_data_per_node), "node idx out of range");
	TDKP_BOUNDS_ASSERT(eqidx >= 0 && eqidx < (unsigned)this->num_data_per_node, "equation idx out of range");	
	return this->data[vidx * this->num_data_per_node + eqidx];  	
}

template<class T>
T& EigenSolution<T>::get_node_value(unsigned int vidx, unsigned int eqidx) {
	return const_cast<T&>(static_cast<const EigenSolution<T>&>(*this).get_node_value(vidx,eqidx)); // use const member function and recast	
}

/** read data from binary stream */
template<class T>
void EigenSolution<T>::read_binary(istream& fin) throw(Exception*) {
	this->reset();
	int num;	
	// read bounding value
 	fin.read((char*)(&num),sizeof(int));
 	if(num != 7041978) {
 		TDKP_GENERAL_EXCEPTION("eigensolution data has invalid header");
 	}	
 	// read data type length
 	fin.read((char*)(&num), sizeof(int));
 	if(num != sizeof(T)) {
 		TDKP_GENERAL_EXCEPTION("data type of stored eigensolution differs from template type"); 		
 	}
 	// read array length
	fin.read((char*)(&num),sizeof(int));
 	fin.read((char*)(&this->energy),sizeof(T));
	fin.read((char*)(&this->num_data_per_node), sizeof(int));	
 	this->data.resize(0);
 	this->data.reserve(num);
 	T tmp;
 	for(int ii = 0; ii < num; ii++) {
 		fin.read((char*)(&tmp), sizeof(T)); 		
		this->data.push_back(tmp);
 	}
 	fin.read((char*)(&num),sizeof(int));
 	if(num != 11092001) {
 		TDKP_GENERAL_EXCEPTION("eigensolution data has invalid end");
 	}
	this->identifier.resize(this->num_data_per_node); 	
}

/** write data to binary stream */
template<class T>
void EigenSolution<T>::write_binary(ostream& fout) const {
	int num = 7041978;
	// store header
	fout.write((char*)(&num),sizeof(int));	
	// store data type
	num = sizeof(T);
	fout.write((char*)(&num),sizeof(int));	
	// store data 	
	num = (int)this->data.size();	
	fout.write((char*)(&num),sizeof(int));
	fout.write((char*)(&this->energy),sizeof(T));
	fout.write((char*)(&this->num_data_per_node), sizeof(int));	
	T tmp;	
	for(int ii = 0; ii < num; ii++) {
		tmp = this->data[ii];
		fout.write((char*)(&tmp),sizeof(T));
	}
	num = 11092001;
	fout.write((char*)(&num),sizeof(int));	
}

/** read binary data from file */
template<class T>
void EigenSolution<T>::read(const char* filename) throw(Exception*){
 	fstream fin(filename, ios::in);
  	if(fin) {
  		this->read_binary(fin);
  		fin.close();
  	} else {
		ostringstream sout;  		
		sout << "could not open file " << filename;
		TDKP_GENERAL_EXCEPTION(sout.str());
  	}
}

/** write binary data to file */
template<class T>
void EigenSolution<T>::write(const char* filename) const  throw(Exception*) {
  	fstream fout(filename, ios::out);
  	if(fout) {
  		this->write_binary(fout);
		fout.close();
  	} else {
  		ostringstream sout;
  		sout << "coult nod write to " << filename;
  		TDKP_GENERAL_EXCEPTION(sout.str());  		
  	}
}

/** compare solution against other eigensolution
 * 
 * energies and vectors are compared against each other. therefore 
 * compare should be called on probability distributions and not 
 * on amplitudes.
 */
template<class T>
bool EigenSolution<T>::compare(const EigenSolution<T>& b, double tol, bool be_chatty) const {
	
	if(!tdkp_math::compare(this->energy, b.energy, tol)) {
		if(be_chatty) {
			ostringstream smsg;
			smsg << "energy of lhs: " << this->energy << " rhs: " << b.energy << "\n";
			Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
		}		
		return false;
	}
	if(this->data.size() != b.data.size()) {
		if(be_chatty) {
			ostringstream smsg;
			smsg << "sizes of vector are not equal!";
			Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
		}		
		return false;	
	}
	for(unsigned int ii = 0; ii < this->data.size(); ii++) {
		if(!tdkp_math::compare(this->data[ii],b.data[ii], tol)) {
			if(be_chatty) {
				ostringstream smsg;
				smsg << "difference at data point: " << ii << " lhs: " 
				     << this->data[ii] << " rhs: " << b.data[ii]
				     << " tol: " << tol;								    
				Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
			}			
			return false;	
		}	
	}
	return true;
}

/** compare probability distributions of eigensolution
 * 
 * due to the degenerate bands we can not compare the 
 * single envelopes, but the probability distributions should
 * be equal up to a certain degree
 * 
 * so if the conversion to a probability eigensolution 
 * should be avoided, one can use this one here 
 */
template<class T> 
bool EigenSolution<T>::compare_probability(const EigenSolution<T>& b, double tol, bool be_chatty) const {
	
	if(!tdkp_math::compare(this->energy, b.energy, tol)) {
		if(be_chatty) {
			ostringstream smsg;
			smsg << "energy of lhs: " << this->energy << " rhs: " << b.energy << "\n";
			Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
		}		
		return false;
	}
	if(this->data.size() != b.data.size()) {
		if(be_chatty) {
			ostringstream smsg;
			smsg << "sizes of vector are not equal!";
			Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
		}		
		return false;	
	}
	const int length = this->get_length();

	T lhs, rhs, tmp1, tmp2;
	for(int ii = 0; ii < length; ii++) {
		lhs = rhs = 0.0;
		// build probability
		for(int ee = 0; ee < this->num_data_per_node; ee++) {
			tmp1 = tdkp_math::abs(this->get_node_value(ii,ee));
			tmp2 = tdkp_math::abs(b.get_node_value(ii,ee));
			lhs += tmp1 * tmp1;
			rhs += tmp2 * tmp2;
		}
		// compare
		if(!tdkp_math::compare(lhs,rhs,tol)) {
			if(be_chatty) {
				ostringstream smsg;
				smsg << "difference at data point: " << ii << " lhs: " 
				     << lhs << " rhs: " << rhs
				     << " tol: " << tol;
				Logger::get_instance()->emit(LOG_INFO_DEVEL1, smsg.str());
			}
			return false;	
		}
	} 
	return true;	
}



template<class F>
void EigenSolution<F>::reset() {
	this->data.clear();
	this->data.resize(0);	
}

/** print information on eigensolution */
template<class T> 
void EigenSolution<T>::info() const {
	if(this->get_length() > 0) {
		ostringstream sout;
		sout << " state with energy " << this->get_energy() << ", defined on " << this->get_length() 
		     << " nodes, " << this->get_num_data_per_node() << " sets per node";
		Logger::get_instance()->emit(LOG_INFO, sout.str());
	} else {
		Logger::get_instance()->emit(LOG_INFO, "eigensolution is uninitialized");	
	}
}

/** return identifier of data */
template<class T> 
string EigenSolution<T>::get_identifier(unsigned int eqidx) const {

	TDKP_ASSERT(eqidx < this->identifier.size(), "identifier index out of range");
	if(this->identifier[eqidx].size() == 0) {
		ostringstream sout;
		sout << "band_" << eqidx;
		return sout.str();	
	} else {
		return this->identifier[eqidx];	
	}	
}

/** set some datas identifer */
template<class T> 
void EigenSolution<T>::set_identifier(unsigned int eqidx, const string& ident) {
	TDKP_ASSERT(eqidx < this->identifier.size(), "identifier index out of range");	
	this->identifier[eqidx] = ident;	
}

/** deletes eigensolutions */
template<class T> 
EigenSolutionSet<T>::~EigenSolutionSet() {
	for(data_iterator it = this->data.begin(); it != this->data.end(); it++) {
		delete (*it);	
	}
};

/** add eigensolution to set */
template<class T> 
void EigenSolutionSet<T>::add(EigenSolution<T>* solution) {
	if(this->data.size() > 0) {
		if(this->data[0]->get_num_data_per_node() != solution->get_num_data_per_node() 
		|| this->data[0]->get_length() != solution->get_length()) {
			TDKP_GENERAL_EXCEPTION("solutions in a set must be of same size");
		}	
	} else {
		this->num_data_per_node = solution->get_num_data_per_node();
	}
	this->data.push_back(solution);	
}

/** get eigensolution from set */
template<class T> 
EigenSolution<T>& EigenSolutionSet<T>::get(unsigned int idx) {
	return const_cast<EigenSolution<T>&>(static_cast<const EigenSolutionSet<T>&>(*this).get(idx));
}

/** get eigensolutions from set */
template<class T>
const EigenSolution<T>& EigenSolutionSet<T>::get(unsigned int idx) const {
	TDKP_ASSERT(idx < this->data.size(), "idx is out of range")
	return *this->data[idx];	 
}
/** display information on solution set */
template<class T>
void EigenSolutionSet<T>::info() const {
	ostringstream sout;
	sout << "EigenSolution set information: " << this->data.size() << " solutions available\n";
	int count = 0;
	for(data_const_iterator it = this->data.begin(); it != this->data.end(); it++) {
		sout << " solution " << count++ << " with energy " << (*it)->get_energy() << " [eV]\n";
	}
	Logger::get_instance()->emit(LOG_INFO, sout.str());
}

/** read eigensolutions from file */
template<class T> 
void EigenSolutionSet<T>::read(const char* filename) {
 	fstream fin(filename, ios::in);
  	if(fin) {
  		while(fin.peek() != EOF) {
  			EigenSolution<T>* tmp = new EigenSolution<T>();
			tmp->read_binary(fin);
			this->data.push_back(tmp);			
  		}
  		fin.close();
  		if(this->data.size() == 0) {
  			TDKP_GENERAL_EXCEPTION("error, there was no valid eigensolution in file");	
  		} else {
  			this->num_data_per_node = this->data[0]->get_num_data_per_node();	
  		}
  	} else {
		ostringstream sout;  		
		sout << "could not open file " << filename;
		TDKP_GENERAL_EXCEPTION(sout.str());
  	}	
}

/** write eigensolutions to file */
template<class T> 
void EigenSolutionSet<T>::write(const char* filename) const {
 	fstream fout(filename, ios::out);
  	if(fout) {  		
		for(data_const_iterator it = this->data.begin(); it != this->data.end(); it++) {		
			(*it)->write_binary(fout);
  		}
  		fout.close();
  	} else {
		ostringstream sout;  		
		sout << "could not open file " << filename;
		TDKP_GENERAL_EXCEPTION(sout.str());
  	}	
}

/** solution set number of data per node
 * 
 * means sum of num_data_per_node of available solutions 
 */
template<class T> 
int EigenSolutionSet<T>::get_num_data_per_node() const {
	if(this->data.size() > 0) {
		return this->data.size() * this->num_data_per_node;
	} else {
		TDKP_GENERAL_EXCEPTION("solution set is uninitialized");	
	}
}

/** return node value of solution set 
 * 
 * (eqidx -> eq = eqidx modulo num_data_per_node)
 * (eqidx -> solution = (eqidx - eq) / num_data_per_node)
 */
template<class T> 
const T& EigenSolutionSet<T>::get_node_value(unsigned int vidx, unsigned int eqidx) const {	
	int eq  = eqidx % this->num_data_per_node;
	int num = (eqidx - eq) / this->num_data_per_node;
	TDKP_BOUNDS_ASSERT(eq >= 0 && eq < this->num_data_per_node, "equation index out of range");
	TDKP_BOUNDS_ASSERT(num >= 0 && num < (signed)this->data.size(), "solution index out of range");
	return this->data[num]->get_node_value(vidx, eq);	
}

template<class T>
T& EigenSolutionSet<T>::get_node_value(unsigned int vidx, unsigned int eqidx) {
	return const_cast<T&>(static_cast<const EigenSolutionSet&>(*this).get_node_value(vidx,eqidx));
}

/** must be implemented for NodeData but must not be called for solution set */
template<class T> 
void EigenSolutionSet<T>::set_node_value(unsigned int vidx, unsigned int eqidx, T value) {
	TDKP_GENERAL_EXCEPTION("EigenSolutionSet can not be initialized via NodeData interface");	
}

/** return length of eigensolutions */
template<class T> 
int EigenSolutionSet<T>::get_length() const {
	TDKP_ASSERT(this->data.size() > 0, "solution set is uninitialized");
	return this->data[0]->get_length();	
}

/** return identifier for eigensolution/band index */
template<class T> 
string EigenSolutionSet<T>::get_identifier(unsigned int eqidx) const {
	int eq  = eqidx % this->num_data_per_node;
	int num = (eqidx - eq) / this->num_data_per_node;
	TDKP_BOUNDS_ASSERT(eq >= 0 && eq < this->num_data_per_node, "equation index out of range");
	TDKP_BOUNDS_ASSERT(num >= 0 && num < (signed)this->data.size(), "solution index out of range");
	ostringstream sout;
	if(this->num_data_per_node > 1) {
		sout << "sol_" << setfill('0') << setw(2) << num << setfill(' ') << "_" << this->data[num]->get_identifier(eq) << "_" << tdkp_math::only_real(this->data[num]->get_energy());
	} else {
		sout << "sol_" << setfill('0') << setw(2) << num << setfill(' ') << "_" << tdkp_math::only_real(this->data[num]->get_energy());
	}
	return sout.str();
}

/** return the number of eigensolutions available */
template<class T> 
int EigenSolutionSet<T>::get_num_solutions() const {
	return this->data.size();	
}

} // end namespace tdkp

#endif /*EIGENSOLUTION_H_*/
