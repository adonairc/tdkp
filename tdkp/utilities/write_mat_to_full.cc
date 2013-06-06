#include <sys/time.h>
#include <math.h>
#include <stdio.h>	
#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "tdkp/main/CSRMatrix.h"

using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace std;
using namespace tdkp;

void convert_matrix_complex(const string& filename_in, const string& filename_out, int precision ) {

 
  CSRMatrix<complex<double> > matrix(filename_in.c_str());
  complex<double>* full = matrix.create_full_matrix();
  int size = matrix.get_size();
  // write real part
  ostringstream sout;
  sout << "real_" << filename_out;  
  fstream fout(sout.str().c_str(), ios::out);
  fout.precision(precision);
  for(int ii = 0; ii < matrix.get_size(); ii++) {
    for(int jj = 0; jj < matrix.get_size(); jj++) {
      fout << full[ii * matrix.get_size() + jj].real() << " ";
    }
    fout << "\n";
  }
  fout.close();
  sout.str("");
  sout << "imag_" << filename_out;
  fout.open(sout.str().c_str(), ios::out);
  fout.precision(precision + 1);
  for(int ii = 0; ii < matrix.get_size(); ii++) {
    for(int jj = 0; jj < matrix.get_size(); jj++) {
      fout << full[ii * matrix.get_size() + jj].imag() << " ";
    }
    fout << "\n";
  }
  fout.close();
  
  delete[] full;
 
};

void convert_matrix_real(const string& filename_in, const string& filename_out, int precision) {

  CSRMatrix<double> matrix(filename_in.c_str());
  double* full = matrix.create_full_matrix();
  int size = matrix.get_size();
  // write real part
  fstream fout(filename_out.c_str(), ios::out);
  fout.precision(precision);
  for(int ii = 0; ii < matrix.get_size(); ii++) {
    for(int jj = 0; jj < matrix.get_size(); jj++) {
      fout << full[ii * matrix.get_size() + jj] << " ";
    }
    fout << "\n";
  }
  fout.close();

}

int main(int argc, char* argv[]) {

  if(argc != 4) {
    cout << "usage: sparse2full <precision> <c|r> filename\n";
    return 0;
  }
    
  string filename_in(argv[3]);
  string type(argv[2]);
  int precision = lexical_cast<int>(argv[1]);
  cout << "using " << precision << " precision\n";
  if(type.compare("c") == 0) {
    ostringstream sout;
    sout << "full_" << filename_in;
    convert_matrix_complex(filename_in, sout.str(), precision);
  } else if(type.compare("r") == 0) {
    ostringstream sout;
    sout << "full_" << filename_in;
    convert_matrix_real(filename_in, sout.str(), precision);
  } else {
    cout << "r = real, c = complex required\n";
    return -1;
  }


}
