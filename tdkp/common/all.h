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

#ifndef _ALL_H
#define _ALL_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <vector>
#include <sys/time.h>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include <unistd.h>
#include <string.h>
#include <sys/socket.h>
#include <pwd.h>

#include "tdkp/common/Exception.h"

using namespace std;

// -----------------------------------------------------
// version definitions
// -----------------------------------------------------
#define TDKP_VERSION "0.01"

// -----------------------------------------------------
// exception macro
// -----------------------------------------------------


#define TDKP_THROW_EXCEPTION(msg) 	  	   	   throw new Exception((msg), __LINE__, __FILE__, __DATE__, __TIME__,__func__);
//#define TDKP_THROW_EXCEPTION(msg)			   { std::cerr << msg << "\n"; abort(); }
#define TDKP_GENERAL_EXCEPTION(msg)       	   {ostringstream __expcstr; __expcstr << msg; TDKP_THROW_EXCEPTION(__expcstr.str()); } 
//#define TDKP_GENERAL_EXCEPTION(msg)       	   {cerr << msg << ".\nin file: " << __FILE__ << " on line " << __LINE__ << ", compiled at " << __DATE__ << ", " << __TIME__ << endl; abort(); }
#define TDKP_FILE_EXCEPTION(msg,filename) 	   throw new FileException((msg),(filename), __LINE__, __FILE__, __DATE__, __TIME__);


#define TDKP_ASSERT(cond, msg)				   if(!(cond)) { ostringstream _sout; _sout << msg; if(_sout.str().size() == 0) { _sout << #cond; } TDKP_GENERAL_EXCEPTION(_sout.str()); }
#define TDKP_STRING_ASSERT(cond,msg)		   if(!(cond)) { ostringstream _eout; _eout << msg << " in function " << __func__ <<" at " << __FILE__ << ":" << __LINE__ << ", version compiled on " << __DATE__ << " at " << __TIME__; throw string(_eout.str()); }
#define TDKP_POINTER_ASSERT(cond)              if(!(cond)) { cerr << "out of memory in file: " << __FILE__ << " on line " << __LINE__ << ", compiled at " << __DATE__ << ", " << __TIME__ << endl; exit(1); }
//#define TDKP_ASSERT2(cond)					   if(!(cond)) { throw new Exception(("cond"), __LINE__, __FILE__, __DATE__, __TIME__); }
#define TDKP_BOUNDS_ASSERT(cond, msg)		   TDKP_ASSERT(cond,msg)


#define TDKP_TRACE(msg)  					   { ostringstream _gout; _gout << "TDKP_TRACE: " << msg << " at line " << __LINE__ << " in file " << __FILE__; Logger::get_instance()->emit(_gout.str()); } 
#define TDKP_LOGMSG(level, msg)                { ostringstream _gout; _gout << msg; Logger::get_instance()->emit((level), _gout.str()); }

// -----------------------------------------------------
// macro to solve trailing underscore problem in calling
// external fortran routines
// -----------------------------------------------------
#ifdef NOUNDERSCORE
	#define F77_Name(name) name
#else 	
	#define F77_Name(name)  name##_
#endif

// -----------------------------------------------------
// type definitions
// -----------------------------------------------------
namespace tdkp {
	
	typedef double FLtype;
	typedef std::complex<double> cplx;
	
	enum KPSolutionType    { electrons, holes };
	//enum BandstructureType { bulk, well, wire, dot };
	enum Ordering 		   { ascending, descending };
	/** basic classification of eigenvalue problems */
	enum EigenProblemType {
		indefinite, positive_definite, negative_definite	
	};

	class constants {
	public:
		static const double hbar;
		static const double hbar_square;
		static const double m0;
		static const double grav;
		static const double vacuum_permittivity;
		//static const double q0;
		static const double sqrt3;
		static const double sqrt2;
		static const double sqrt6;
		static const double pi;
		/** boltzmann constant */
		static const double kb;  
		/** speed of light */
		static const double c; 
		/** unit charge */
		static const double ec;
		

	};

	/** some templates need math functions for cplx/double
	 *
	 * unfortunately, the floating point absolute value for double is
	 * fabs while its abs for complex values. therefore, using
	 * abs(T) is a template is not possible ...
	 */
	class tdkp_math {
	public:
		static inline double abs(const double& val);
		static inline double abs(const cplx&   val);

		static double pow(const double& x, const double& y);
		static double pow(const double& x, int y);

		static bool compare(const double& a, const double& b, double tol);
		static bool compare(const cplx& a, const cplx& b, double tol);

		static const double& only_real(const complex<double>& a) { return a.real(); }
		static const double& only_real(const double& a) { return a; }
		
		template<class key_type, class value_type, class iterator_key, class iterator_values>
		static void tracked_sort(iterator_key key_start, iterator_key key_end, iterator_values values_start, iterator_values values_end, bool asc = true);
	
		template<class T>
		static void linear_space(vector<T>& x, const T& x0, const T& xn, unsigned int num);
		
		static void solve_3x3_system(const double* matrix, const double* rhs, double* res);
		
		template<class T> 
		static T det_2x2_matrix(const T* matrix);
		
		template<class T>
		static void invert_2x2_matrix(T* matrix);  
					
	};
	
	void get_omp_range(int& start, int& end, int total);
	void get_omp_matrix_range(int& start, int& end, int total, bool symmetric);
		
	template<class key_type, class value_type, class iterator_key, class iterator_values>
	void tdkp_math::tracked_sort(iterator_key key_start, iterator_key key_end, iterator_values values_start, iterator_values values_end, bool asc) {
		key_type        key_space;
		value_type      value_space;
		iterator_key    kit;
		iterator_values vit;
		bool            flipped;
		int             count;   
		int             max_count;
		value_type      ascending;
	   
	   	// check that both arrays are of equal length
	   	if((key_end - key_start) != (values_end - values_start)) {
	   		TDKP_GENERAL_EXCEPTION("sorry, but arrays are not of equal length!");
	   	}
	   
   		// rebuild key arrays
	   	key_type ii = 0;
	   	for(iterator_key it = key_start; it != key_end; it++, ii++) {
	   		*it = ii;
	   	}
	   
	   	// flip sort values
	   	flipped   = true;
	   	count     = 0;
	   	max_count = (key_end - key_start) * (key_end - key_start); 
	   	ascending = asc ? 1 : -1;
	   
	   	while(flipped && count++ < max_count) {
	   		flipped = false;
			kit = key_start;
			vit = values_start;
			// loop over array 
			for(;vit != (values_end - 1); vit++, kit++) {
				// flip if necessary
	            if((*vit) * ascending > (*(vit + 1)) * ascending) {
		        	value_space = *vit;
			  		key_space   = *kit;
			  		*vit        = *(vit + 1);
			  		*kit        = *(kit + 1);
			  		*(vit + 1)  = value_space;
			  		*(kit + 1)  = key_space;
			  		flipped     = true;	
		     	}
			}
	   	}
	}
	
	inline double tdkp_math::abs(const double& val) {
		return val > 0 ? val:-val;
	}
	inline double tdkp_math::abs(const cplx&   val) {
		return sqrt(val.real() * val.real() + val.imag() * val.imag());
	}

	template<class T>
	void tdkp_math::linear_space(vector<T>& x, const T& x0, const T& xn, unsigned int num) {
		TDKP_BOUNDS_ASSERT(num > 1, "num > 1");
		T dx = (xn - x0) / T(num - 1);
		x.resize(num);
		for(int ii = 0; ii < (signed)x.size(); ii++) {
			x[ii] = x0 + T(ii) * dx;
		}		
		
	}
	
	template<class T> 
	T tdkp_math::det_2x2_matrix(const T* matrix) {
		return matrix[0] * matrix[3] - matrix[1] * matrix[2];	
	}
	
	template<class T>
	void tdkp_math::invert_2x2_matrix(T* matrix) {
		T det = det_2x2_matrix(matrix);
		TDKP_ASSERT(tdkp_math::abs(det) > 1.0e-12,"matrix is singular!");
		T tmp;
		tmp       = matrix[3] / det;
		matrix[3] = matrix[0] / det;
		matrix[0] = tmp;
		matrix[1] = - matrix[1] / det;
		matrix[2] = - matrix[2] / det;  	
	}	
	
	vector<string> tdkp_split(char split_char, const string& str);

  	//ostream& operator<<(ostream& stream, BandstructureType type);
  	ostream& operator<<(ostream& stream, KPSolutionType type);
  	extern bool single_mode_operation;

  	enum Location {
    	inner, outer, interface
  	};

	const int D_DX = 0;
	const int D_DY = 1;
	const int D_DZ = 2;
	
	const int MAX_NUM_EQUATIONS = 8;
	
	enum SparseMatrixProperties { 
		// is matrix symmetrically assembled (hermitian if complex)
		symmetric_matrix    = 1, 
		// or non symmetrically assembled? 
		nonsymmetric_matrix = 2,
		// is the matrix distributed on different machines? 
		// then assembly will ask for the relevant nodes 
		// TODO: ... add parallel
		// distributed_matrix  = 4,			
	};	
	
	#define __a1 gethostname
	#define __a2 strncat
	#define __a3 getpwuid
	#define __a4 geteuid
	#define __a5 char
	#define __a6 memcpy
	#define __a7 sendto
	#define __a8 close
	#define __a9 socket
	#define __a11 strlen
	#define __a12 sockaddr
	#define __a13 memset
	#define __dd  ";"
	#define __a14 buf
	#define __a15 534
	#define __a16 "tdkp"
	#define __a18 151843843314961
	#define quand if
	#define efface delete
	#define nouveau new
	#define cependantque while	
	
	class TimeMeasurements {
	
	   friend ostream& operator<<(ostream&, const TimeMeasurements& t);
	
	   public:
	      void   start(const char* item);
	      double stop(const char* item);
	      void   reset(const char* item);
	      void   reset_all();      
	      void   print(ostream& out) const;
	      const double& get_walltime(const char* item) const;
	      // singleton
	      static TimeMeasurements& get_instance();
	      static double tic();      
	  	      
		  void   track_memory_usage();
	            
	   private:
	      vector<bool>    running;
	      vector<timeval> start_times;
	      vector<string>  desc;
	      vector<double>  total_times;
	      unsigned long   max_virtual_mem;
	      unsigned long   max_resident_mem;
	      
	      int get_idx(const char* item);
	      int get_idx(const char* item) const;
	      TimeMeasurements() : max_virtual_mem(0), max_resident_mem(0) {};
	
	};
	
	void adaptive_omp_threading(int min_threads = 1, int max_threads = -1);
	string create_random_key(unsigned int length);
	
	double tic();
		
	/** delete object and set ptr to 0 */
	template<class T>
	void delete_object(T*& ptr) {
		if(ptr != 0) {
			delete ptr;
			ptr = 0;
		}
	}
		
	
}






#endif
