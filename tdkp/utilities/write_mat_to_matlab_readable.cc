#include <sys/time.h>
#include <math.h>
#include <stdio.h>	
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <complex>

#include "tdkp/common/all.h"
#include "tdkp/main/CSRMatrix.h"

#ifndef SHOW_NONZEROS
#define SHOW_NONZEROS true
#endif

using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace tdkp;
using namespace std;

string readable(cplx a, int precision) {
  ostringstream sout;
  sout.precision(precision);
  sout << a.real() << " " << a.imag();
  return string(sout.str());
}

string readable(double a, int precision) {
  ostringstream sout;
  sout.precision(precision);
  sout << a;
  return string(sout.str());
}

double myabs(cplx a) {
  return abs(a);
}

double myabs(double a) {
  return fabs(a);
}

template<class T>
void convert_matrix(const string& filename_in, const string& filename_out, int precision ) {
 
  	CSRMatrix<T> matrix(filename_in.c_str());  	
  	int size = (int)matrix.get_size();
  	int fidx = (int)matrix.get_fidx();
  	int* prow   = matrix.get_prow();
  	int* icol   = matrix.get_icol();
  	T* nonzeros = matrix.get_nonzeros();
  	fstream fout(filename_out.c_str(), ios::out);
  	fout.precision(precision);
	int from, until;
  	for(int ii = 0; ii < size; ii++) {
		until    = prow[ii + 1] - fidx;
		from     = prow[ii] - fidx;						
		for(int jj = from; jj < until; jj++) {	
			if(tdkp_math::abs(nonzeros[jj]) != 0.0 || SHOW_NONZEROS) { 
				fout << ii << " " << icol[jj] - fidx << " " << readable(nonzeros[jj], precision) << "\n";
			}
			
      	}
    	
	}
  	fout.close(); 

}

template<class T>
void convert_matrix_csr(const string& filename_in, const string& filename_out, int precision ) {
 
  	CSRMatrix<T> matrix(filename_in.c_str());
  	T* nonzeros = matrix.get_nonzeros();
  	int* icol   = matrix.get_icol();
  	int* prow   = matrix.get_prow(); 
	int  fidx   = matrix.get_fidx();
//  T* full = matrix.create_full_matrix();
  	int size    = (int)matrix.get_size();
  	fstream fout(filename_out.c_str(), ios::out);
  	fout.precision(precision);
//  fout << "huhu\n";
	
  	int until, from; 

  	if(matrix.symmetric()) {
  		cout << " WARNING: MATRIX IS SYMMETRIC, ONLY GIVING THE UPPER TRIDIAGONAL\n";	
  	}		
	for(int ii = 0; ii < size; ii++) {
		until    = prow[ii + 1] - fidx;
		from     = prow[ii] - fidx;						
		for(int jj = from; jj < until; jj++) {	
			if(tdkp_math::abs(nonzeros[jj]) != 0.0) { 		
				fout << ii << " " << (icol[jj] - fidx) << " " << readable(nonzeros[jj], precision) << "\n";
			}			
		}	
	}	
	

  	fout.close(); 

}



void matlabcode() {
  cout << "\nmatlabcode: \n% complex: \n"
       << "load -ascii conv_lhs.mat\n"
       << "n = max(max(conv_lhs(:,1)),max(conv_lhs(:,2))) + 1;\n"
       << "A = sparse(n,n);\n"
       << "for ii = 1:size(conv_lhs,1)\n"
       << "  A(conv_lhs(ii,1) + 1, conv_lhs(ii,2) + 1) = conv_lhs(ii,3) + i * conv_lhs(ii,4);\n"
       << "end\n";                

}

int main(int argc, char* argv[]) {

  if(argc != 4) {
    cout << "usage: sparse2full <precision> <c|r> filename\n";
    matlabcode(); 
    return 0;
  }
    
  string filename_in(argv[3]);
  string type(argv[2]);
  int precision = lexical_cast<int>(argv[1]);
  cout << "using " << precision << " precision\n";
  if(type.compare("c") == 0) {
    ostringstream sout;
    sout << "conv_" << filename_in;
    convert_matrix<cplx>(filename_in, sout.str(), precision);
  } else if(type.compare("r") == 0) {
    ostringstream sout;
    sout << "conv_" << filename_in;
    convert_matrix<double>(filename_in, sout.str(), precision);
  } else {
    cout << "r = real, c = complex required\n";
    return -1;
  }


}
