#include <sys/time.h>
#include <complex>
#include <math.h>

using namespace std;

complex<double>* create_random_hermitian(unsigned int size) {
	
	double fill = 0.1 / 2.0;
	
	timeval tp;
	gettimeofday(&tp, NULL);
	srand48(tp.tv_usec);
	
	complex<double>* ret = new complex<double>[size * size];
	
	for(unsigned int ii = 0; ii < size; ii++) {
		ret[ii * size + ii] = 5.0 * drand48();
		for(unsigned int jj = ii + 1; jj < size; jj++) {
			if(drand48() < fill) {
				ret[ii * size + jj] = complex<double>((drand48() - 0.5), drand48() - 0.5);
				ret[jj * size + ii] = conj(ret[ii * size + jj]); 	
			}	
		}			
	} 
	return ret;	
}


complex<double>* create_laplace_hermitian(unsigned int n) {		
	complex<double>* matrix = new complex<double>[n * n];
	for(unsigned int ii = 0; ii < n*n; ii++) {
		matrix[ii] = 0;	
	}
	complex<double> i(0.0, 1.0);
	for(unsigned int ii = 0; ii < n; ii++) {
		matrix[ii * n + ii] = 2.0;			
		if(ii < n - 1) {
			// jj * size + ii
			matrix[ii * n + ii + 1]   =   i;
			matrix[(ii + 1) * n + ii] = - i;	
		}
	}
	return matrix;
}

double* create_laplace(unsigned int n) {
	double* matrix = new double[n * n];
	for(unsigned int ii = 0; ii < n*n; ii++) {
		matrix[ii] = 0.0;	
	}
	for(unsigned int ii = 0; ii < n; ii++) {
		matrix[ii * n + ii] = 2.0;			
		if(ii < n - 1) {
			// jj * size + ii
			matrix[ii * n + ii + 1]   =   1;
			matrix[(ii + 1) * n + ii] =   1;	
		}
	}
	return matrix;	
}

template<class T>
void print_matrix(T* matrix, unsigned int size) {
	cout.precision(6);
	cout << "\n";
	for(unsigned int ii = 0; ii < size; ii++) {
		for(unsigned int jj = 0; jj < size; jj++) {
			cout << setw(10) << matrix[ii * size + jj] << " "; 
		}	
		cout << "\n";
	}	
}

template<class T>
T* get_random_vector(unsigned int size) {
	T* ret = new T[size];
	for(uint ii = 0; ii < size; ii++) {
		ret[ii] = drand48();	
	}	
	return ret;
}

template<class T> 
void reorder_for_fortran(T* matrix, unsigned int size) {
	T tmp;
	// copy
	for(unsigned int ii = 0; ii < size; ii++) {
		for(unsigned int jj = ii + 1; jj < size; jj++) {
			tmp = matrix[ii * size + jj];
			matrix[ii * size + jj] = matrix[jj * size + ii];
			matrix[jj * size + ii] = tmp; 		
		}	
	}
}
