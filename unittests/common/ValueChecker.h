#ifndef VALUECHECKER_H_
#define VALUECHECKER_H_

#include "tdkp/common/all.h"
#include <cxxtest/TestSuite.h>
#include <vector>
#include <string>
#include "tdkp/io/InputParser.h"
#include "tdkp/common/Vector3D.h"
#include "tdkp/main/Bandstructure.h"
#include "tdkp/main/EigenSolution.h"
#include "tdkp/common/ICurves.h"


using namespace std;

namespace tdkp {

class ValueCheckerError {
	public:
		ValueCheckerError(const string& massage) : message(massage) {};
		const string& get_message() const { return this->message; }
	private:
		string message;		
};

/** class to compare result vectors */
template<class T> 
class ValueChecker {
public:
	
	static void compare_vectors(const vector<T>& a, const vector<T>& b, const double& max_dev) throw(ValueCheckerError);
	static void compare_vectors(const T* a, const T* b, unsigned int length, const double& max_dev) throw(ValueCheckerError);
	static void compare_vectors_abs(const vector<T>& a, const vector<T>& b, const double& max_dev) throw(ValueCheckerError);
	static void compare_vectors_abs(const T* a, const T* b, unsigned int length, const double& max_dev) throw(ValueCheckerError);
	static void read_reference_vector(vector<T>& in, const char* filename);
	
	static Vector3D get_random_dir();
	
	static void compare_ICurves(const ICurve& left, const ICurve& right, const double& max_dev);
		
	
private:
	ValueChecker() {};
	struct error_container {
		unsigned int idx;
		double off_by_value;
		T value_a, value_b;
	};
	static string assemble_error_message(const vector<error_container>&, const T* a, const T* b, unsigned int length, const double& max_dev);	 	
};
	
template<class T> 
void ValueChecker<T>::compare_vectors(const vector<T>& a, const vector<T>& b, const double& max_dev) throw(ValueCheckerError) {
	TDKP_ASSERT(a.size() == b.size(), "vectors must be of equal size!");
	compare_vectors(&a[0], &b[0], a.size(), max_dev);
}

template<class T> 
void ValueChecker<T>::compare_vectors(const T* a, const T* b, unsigned int length, const double& max_dev) throw(ValueCheckerError) {
	
	double tmp; 	
	double div;			
	error_container error_data;						
	vector<error_container> errors; 
	
	// check all values
	for(unsigned int ii = 0; ii < length; ii++) {
		// don't divide by zero
		if(tdkp_math::abs(a[ii]) > max_dev) {
			div = tdkp_math::abs(a[ii]);	
		} else if(tdkp_math::abs(b[ii]) > max_dev) {
			div = tdkp_math::abs(b[ii]);
		} else if(a[ii] != 0.0 && b[ii] != 0.0) {
			div = 1.0;
		} else {
			div = 0.0;
		}
		if(div != 0.0) {
			tmp = tdkp_math::abs(a[ii] - b[ii]) / div;
			// collect errors
			if(tdkp_math::abs(tmp) > max_dev) {
				error_data.idx          = ii;
				error_data.off_by_value = tmp;
				error_data.value_a      = a[ii];
				error_data.value_b      = b[ii];
				errors.push_back(error_data);	
			}
		}
	}
	// throw error if necessary
	if(errors.size() > 0) {
		throw ValueCheckerError(assemble_error_message(errors, a, b, length, max_dev));	
	}
}

template<class T> 
void ValueChecker<T>::compare_vectors_abs(const vector<T>& a, const vector<T>& b, const double& max_dev) throw(ValueCheckerError) {
	TDKP_ASSERT(a.size() == b.size(), "vectors must be of equal size!");
	compare_vectors_abs(&a[0], &b[0], a.size(), max_dev);
}

template<class T> 
void ValueChecker<T>::compare_vectors_abs(const T* a, const T* b, unsigned int length, const double& max_dev) throw(ValueCheckerError) {
	// dont care about performance ;-)
	vector<T> tmp_a(length);
	vector<T> tmp_b(length);
	for(unsigned int ii = 0; ii < length; ii++) {
		tmp_a[ii] = tdkp_math::abs(a[ii]);
		tmp_b[ii] = tdkp_math::abs(b[ii]);
	}
	compare_vectors(tmp_a, tmp_b, max_dev);
}

template<class T> 
Vector3D ValueChecker<T>::get_random_dir() {
	timeval time;
	gettimeofday(&time, NULL);
	srand48(time.tv_usec);
	Vector3D dir;
	dir(0) = (0.5 - drand48()) * 2.0;
	dir(1) = (0.5 - drand48()) * 2.0;
	dir(2) = (0.5 - drand48()) * 2.0;
	// don't make problems if you are too small ... 
	if(dir.norm() < 1.0e-10) {
		return get_random_dir();	
	} else {
		dir.normalize();
		return dir;
	} 
}

template<class T> 
string ValueChecker<T>::assemble_error_message(const vector<error_container>& errors, const T* a, const T* b, unsigned int length, const double& max_dev) {
	
	// ----------------------------------------
	// assemble statistics
	// ----------------------------------------
	double value_a_average        = 0.0;
	double value_a_square_average = 0.0;
	double value_b_average        = 0.0;
	double value_b_square_average = 0.0;
	double error_average          = 0.0;
	double error_square_average   = 0.0;

	double abs_a, abs_b, diff_abs;
	for(unsigned int ii = 0; ii < length; ii++) {
		abs_a    = tdkp_math::abs(a[ii]);
		abs_b    = tdkp_math::abs(b[ii]);
		diff_abs = tdkp_math::abs(a[ii] - b[ii]);
		value_a_average        += abs_a;
		value_b_average        += abs_b;
		error_average          += diff_abs;
		value_a_square_average += abs_a * abs_a; 
		value_b_square_average += abs_b * abs_b;
		error_square_average   += diff_abs * diff_abs;
		cout << "a: \t" << a[ii] << " \tb: \t"<< b[ii] << "\n";
	}
	
	cout << "LENGTH: " << length << endl;
	
	// -----------------------------------------
	// find max error
	// -----------------------------------------
	typename vector<error_container>::const_iterator it = errors.begin();	
	unsigned int idx, max_idx;
	max_idx = idx = 0; 
	while(idx < errors.size()) {
		if(errors[max_idx].off_by_value < errors[idx].off_by_value) {
			max_idx = idx;	
		}
		idx++;
	}
	
	unsigned int field_size = 23;
	
	ostringstream sout;
	sout << " ----------------------------------------------------- \n"
	     << "  tdkp ValueChecker found some differences in two vecs \n"
	     << "    errors: " << errors.size() << " out of " << length << " data points\n"
	     << "\n"
	     << setw(field_size) << " " 
	     << setw(field_size) << "vector a "  
	     << setw(field_size) << "vector b " << "\n"
	     << setw(field_size) << "abs average: "
	     << setw(field_size) << value_a_average / double(length)
		 << setw(field_size) << value_b_average / double(length) << "\n"
	     << setw(field_size) << "deviation: "
	     << setw(field_size) << sqrt(fabs((value_a_square_average - value_a_average * value_a_average) / double(length))) 
	     << setw(field_size) << sqrt(fabs((value_b_square_average - value_b_average * value_b_average) / double(length)))
	     << "\n"
	     << setw(field_size) << "error average: " 
	     << error_average / double(length) << "\n"
	     << setw(field_size) << "error deviation: " 
	     << sqrt(fabs((error_square_average - error_average*error_average) / double(length)))
	     << "\n"
	     << "\n"
	     << " max error at: " << errors[max_idx].idx << "," 
	     << " error: " << errors[max_idx].off_by_value << ", "
	     << " left: "  << errors[max_idx].value_a << ","
	     << " right: " << errors[max_idx].value_b << "\n";
	     
	return sout.str();
	     
}


template<class T>
void ValueChecker<T>::compare_ICurves(const ICurve& left, const ICurve& right, const double& threshold) {
	vector<double> ldata, rdata;
	
	const vector<double>& grid = left.GetGridReference();
	left.EvalAt(grid, ldata);
	right.EvalAt(grid, rdata);
	compare_vectors(&ldata[0], &rdata[0], grid.size(), threshold);
		
}

class ResultTestBase {
protected:
	
	MaterialDatabase matdb;
	std::fstream fout;		
	bool write_mode;	
	const double kmin, kmax;
	const unsigned int knum;
	const Vector3D one_one_one;
	const Vector3D one_one_zero;	
	const Vector3D one_zero_zero;
	const string my_directory;
	
	bool compare_bandstructures(const char* message, const char* filename, const BandstructureDomain<complex<double> >& values) const;	
	bool compare_bandstructures(const char* message, const BandstructureDomain<complex<double> >& reference, const BandstructureDomain<complex<double> >& values) const;
	
	ResultTestBase(const char* my_directory_);
	
};


// ---------------------------------------------------------
// ResultTestBase implementation
// ---------------------------------------------------------
ResultTestBase::ResultTestBase(const char* my_directory_)
: write_mode(false),
  kmin(0.0),
  kmax(2.0*constants::pi / 0.5),   
  knum(60),
  one_one_one  (1.0, 1.0, 1.0),
  one_one_zero (1.0, 1.0, 0.0),
  one_zero_zero(1.0, 0.0, 0.0),
  my_directory(my_directory_) 
{
	matdb.add_search_path(my_directory.c_str());
	matdb.load_material("GaAs");
	matdb.load_material("InAs");
	
	if(write_mode) {
		ostringstream sout;
		sout << my_directory << "data/DATAUPDATE";
		ofstream fout(sout.str().c_str(), ios::out);
		if(fout) {
			time_t my_time = time(NULL);
			fout << "The reference data in this directory was updated:\n"
			     << ctime(&my_time) << endl;
			fout.close(); 
		} else {
			TDKP_GENERAL_EXCEPTION("can not write to DATAUPDATE file!");	
		}
	}		
}

bool ResultTestBase::compare_bandstructures(const char* message, const BandstructureDomain<complex<double> >& reference, const BandstructureDomain<complex<double> >& values) const {
		
	bool is_equal = true;		
		
	vector<complex<double> > lhs;
	vector<complex<double> > rhs;
	// -------------------------------------------
	// check general data
	// -------------------------------------------
	TSM_ASSERT_EQUALS(message, values.get_number_of_bands(), reference.get_number_of_bands());
	TSM_ASSERT_EQUALS(message, values.get_basis_size(), reference.get_basis_size());
	TSM_ASSERT_EQUALS(message, values.get_number_of_k_values(), reference.get_number_of_k_values());
	TSM_ASSERT_EQUALS(message, values.get_solution_length(), reference.get_solution_length());
	
				
	// -------------------------------------------
	// compare dispersion relations (but only the the first 2/3
	// -------------------------------------------
	const int num_bands_to_compare = int(ceil(2.0/3.0*reference.get_number_of_bands()));
	TDKP_ASSERT(num_bands_to_compare > 0 && num_bands_to_compare < reference.get_number_of_bands(), "num_bands_to_compare > 0 && num_bands_to_compare < reference.get_number_of_bands()"); 		
	for(int nn = 0; nn < num_bands_to_compare; nn++) {
		lhs.resize(reference.get_number_of_k_values());
		rhs.resize(reference.get_number_of_k_values());
		
		for(int kk = 1; kk < reference.get_number_of_k_values(); kk++) {				
			lhs[kk] = values.get_energy(kk, nn);
			rhs[kk] = reference.get_energy(kk,nn);		
			TS_ASSERT_EQUALS(isnan(lhs[kk].real()), 0);
			TS_ASSERT_EQUALS(isnan(rhs[kk].real()), 0);
			TS_ASSERT_EQUALS(isnan(lhs[kk].imag()), 0);
			TS_ASSERT_EQUALS(isnan(rhs[kk].imag()), 0);		
			if(!values.get_eigensolution(kk,nn).compare_probability(reference.get_eigensolution(kk,nn))) {
				ostringstream sout;
				sout << message << " probability differs for band " 
				     << nn << " at kk " << kk; 				     
				TS_FAIL(sout.str());
				is_equal = false;				
			}			
		} 	
	
		try {			
			ValueChecker<complex<double> >::compare_vectors(lhs,rhs,1.0e-6);
		} catch(ValueCheckerError& e) {
			ostringstream sout;
			sout << message << " bandstructure differs for band " 
				 << nn << "\n"
				 << e.get_message();
			//delete e;					 					     
			TS_FAIL(sout.str());
			Logger::get_instance()->emit(LOG_INFO, sout.str());
			is_equal = false;	
		}				
	}
	if(!is_equal) {
		ostringstream ssfname, sout;
		ssfname << message << "_my_xy.dat";
		sout << "writing bandstructure (and ref) to files " 
		     << ssfname.str() << " for inspection"; 		 
		TS_WARN(sout.str());
		InputParser parser;
		parser.write_ascii_real(values, ssfname.str().c_str());
		ssfname.str("");
		ssfname << message << "_ref_xy.dat";
		parser.write_ascii_real(reference, ssfname.str().c_str());
	}
	
	return is_equal;	
}

bool ResultTestBase::compare_bandstructures(const char* message, const char* filename, const BandstructureDomain<complex<double> >& values) const {
		
	// reference data file
	ostringstream ssfname;
	ssfname << my_directory << "data/" << filename; 
	// write
	if(write_mode) {		
		values.write_binary(ssfname.str().c_str()); 
		cout << "updated file " << ssfname.str() << "\n";		
		return true;
	// compare		
	} else {
		// read reference data
		BandstructureDomain<complex<double> > reference(ssfname.str().c_str());
		return compare_bandstructures(message, reference, values);
	}
	
}



} // end of namespace

#endif /*VALUECHECKER_H_*/
