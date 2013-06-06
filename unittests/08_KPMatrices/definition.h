
#define TDKP_UNITTEST

// ------------------------------------------
// standard includes
// ------------------------------------------
#include <cxxtest/TestSuite.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
// ------------------------------------------
// standard tdkp includes and logger def
// ------------------------------------------
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"
#include "tdkp/common/Exception.h"

// ------------------------------------------
// class includes
// ------------------------------------------
#include <map>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <vector>
#include "tdkp/common/Vector3D.h"
#include "tdkp/main/PropertyContainer.h"
#include "tdkp/main/MaterialDatabase.h"
#include "tdkp/main/Fields.h"
#include "tdkp/kpmatrices/KPMatrix4x4EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix6x6EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix8x8EndersForeman.h"
#include "tdkp/kpmatrices/KPMatrix6x6Wurtzite.h"
#include "tdkp/kpmatrices/KPMatrix8x8Wurtzite.h"
#include "ValueChecker.h"


using namespace tdkp;

extern "C" {
void zheevx_(char* jobz, char* range, char* uplo, int* n, cplx* A, int lda, double* VL,double* VU,
			 int* IL, int* IU, double* ABSTOL, int* M, double* W,cplx* Z, int* LDZ, cplx* WORK, int* LWORK, double* RWORK,
             int* IWORK, int* IFAIL, int* INFO );	     
}

class KPMatrixTest : public CxxTest::TestSuite {
public:	
		
	KPMatrixTest() { 		
		fout.open("08_KPMatrices_output.log", ios::app | ios::out);
		Logger::get_instance()->set_level(LOG_INFO_DEVEL2);		
		Logger::get_instance()->add_listener(&fout);		
		Logger::get_instance()->del_listener(&std::cout);	
		this->material_db = new MaterialDatabase();
		this->material_db->add_search_path("./08_KPMatrices/");
		try {
			this->material_db->load_material("GaAs");
			this->material_db->load_material("InAs");
			this->material_db->load_material("GaN");
		} catch (Exception*e) {
			cout << e->get_reason() << endl;
			throw e;	
		}				
	}	
	
	~KPMatrixTest() {
		fout.close();
		delete this->material_db;				
	}	
		
	MaterialDatabase* material_db;						
	std::fstream fout;	
					
	void test_sparsity_matrix_diagonal() {
		int diagonal_pattern[] = {
			0,0, 
				1,1,
					2,2,
						3,3,
							4,4,
		};
		vector<double> test;
		test.resize(5,1.0);
		RMatrix<double> matrix(test,diagonal_pattern,5);		
		TS_ASSERT_EQUALS((signed)matrix.rows(), 5);
		for(int ii = 0; ii < 5; ii++) {
			TS_ASSERT_EQUALS(matrix(ii,ii),1.0);
			for(int jj = ii + 1; jj < 5; jj++) {
				TS_ASSERT_EQUALS(matrix(ii,jj),0.0);		
				TS_ASSERT_EQUALS(matrix(jj,ii),0.0);
			}	
		}					
	}
	
	void test_sparsity_matrix_full() {

		int pattern[] = {
			0,0,    0,2,
				1,1,    1,3,
			2,0,	2,2,    2,4,
				3,1,	3,3,
					4,2,	4,4,
		};		
		double values[] = {
			1.0, 2.0, 3.0, 4.0, 2.0, 5.0, 6.0, 4.0, 7.0, 6.0, 8.0
		};
		vector<double> test;
		for(int ii = 0; ii < 11; ii++) {
			test.push_back(values[ii]);	
		}
		RMatrix<double> matrix(test, pattern, 11);
		TS_ASSERT_EQUALS((signed)matrix.rows(), 5);
		TS_ASSERT_EQUALS((signed)matrix.cols(), 5);
		for(int ii = 0; ii < 5; ii++) {
			for(int jj = ii + 1; jj < 5; jj++) {
				TS_ASSERT_EQUALS(matrix(ii,jj),matrix(jj,ii));		
			}	
		}					
	}
	
	StrainTensor get_strain_tensor() {
		StrainTensor strain;
		strain.set(iexx, 0.05);
		strain.set(iexx, 0.06);		
		strain.set(iexx, 0.03);		
		strain.set(iexy, -0.01);
		strain.set(iexz, -0.02);		
		strain.set(ieyz, 0.04);			
		return strain;
	}


	void test_kp_matrix_4x4_strainless_gaas() {
		KPMatrix4x4EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		perform_strainless_test(matrix, "KPMatrix4x4EndersForeman GaAs");		
	}

	void test_kp_matrix_4x4_final_gaas() {
		KPMatrix4x4EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		this->perform_final_matrix_test(matrix, "KPMatrix4x4EndersForeman GaAs");		
	}		
	
	void test_kp_matrix_4x4_straindependent_gaas() {
		KPMatrix4x4EndersForeman matrix;				
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		this->perform_straindependent_matrix_test(matrix, "KPMatrix4x4EndersForeman GaAs");								
	}	
	
	void test_kp_matrix_4x4_final_strain_gaas() {
		KPMatrix4x4EndersForeman matrix;							
		StrainTensor strain = this->get_strain_tensor();					
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(strain));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix4x4EndersForeman GaAs");		
	}

	void test_kp_matrix_6x6_wurtztite_strainless_gan() {
		KPMatrix6x6Wurtzite matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		perform_strainless_test(matrix, "KPMatrix6x6Wurtzite GaN");		
	}

	void test_kp_matrix_6x6_wurtztite_final_gan() {
		KPMatrix6x6Wurtzite matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix6x6Wurtzite GaN");		
	}
			
	void test_kp_matrix_6x6_wurtzite_straindependent_gan() {
		KPMatrix6x6Wurtzite  matrix;				
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(this->get_strain_tensor()));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());					
		this->perform_straindependent_matrix_test(matrix, "KPMatrix6x6Wurtzite  GaN");								
	}	
	
	void test_kp_matrix_6x6_wurtzite_final_strain_gan() {
		KPMatrix6x6Wurtzite  matrix;							
		StrainTensor strain = this->get_strain_tensor();					
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(strain));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix6x6Wurtzite GaN");		
	}	
	
	void test_kp_matrix_8x8_wurtztite_strainless_gan() {
		KPMatrix8x8Wurtzite matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		perform_strainless_test(matrix, "KPMatrix8x8Wurtzite GaN");		
	}

	void test_kp_matrix_8x8_wurtztite_final_gan() {
		KPMatrix8x8Wurtzite matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaN")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix8x8Wurtzite GaN");		
	}	



	void test_kp_matrix_6x6_strainless_gaas() {
		KPMatrix6x6EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		perform_strainless_test(matrix, "KPMatrix6x6EndersForeman GaAs");		
	}

	void test_kp_matrix_6x6_final_gaas() {
		KPMatrix6x6EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix6x6EndersForeman GaAs");		
	}		
	
	void test_kp_matrix_6x6_straindependent_gaas() {
		KPMatrix6x6EndersForeman matrix;				
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(this->get_strain_tensor()));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());					
		this->perform_straindependent_matrix_test(matrix, "KPMatrix6x6EndersForeman GaAs");								
	}	
	
	void test_kp_matrix_6x6_final_strain_gaas() {
		KPMatrix6x6EndersForeman matrix;							
		StrainTensor strain = this->get_strain_tensor();					
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(strain));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix6x6EndersForeman GaAs");		
	}	
	
	void test_kp_matrix_8x8_strainless_gaas() {
		KPMatrix8x8EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		perform_strainless_test(matrix, "KPMatrix8x8EndersForeman GaAs");		
	}

	void test_kp_matrix_8x8_final_gaas() {
		KPMatrix8x8EndersForeman matrix;								
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix8x8EndersForeman GaAs");		
	}		
	
	void test_kp_matrix_8x8_straindependent_gaas() {
		KPMatrix8x8EndersForeman matrix;				
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		this->perform_straindependent_matrix_test(matrix, "KPMatrix8x8EndersForeman GaAs");								
	}	
	
	void nnotest_kp_matrix_8x8_final_strain_gaas() {
		KPMatrix8x8EndersForeman matrix;							
		StrainTensor strain = this->get_strain_tensor();					
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());			
		TS_ASSERT_THROWS_NOTHING(matrix.set_strains(strain));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		this->perform_final_matrix_test(matrix, "KPMatrix8x8EndersForeman GaAs");		
	}	
	
	void perform_strainless_test(KPMatrixBase& matrix, const char* classname) {
	
		RMatrix<cplx>* tmp = 0;
		RMatrix<cplx>* tmp2 = 0;		
										
		ostringstream sout;		
				
		cplx kvalues[][3] = {
			{0.0, 0.0, 0.0},			
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0},
			{1.0, 1.0, 0.0},
			{1.0, 0.0, 1.0},
			{0.0, 1.0, 1.0},
			{1.0, 1.0, 1.0}
		};
		
		int   num_k_sets = 8;
		int   size = matrix.get_number_of_bands();
		cplx  i(0.0, 1.0);
		cplx* tmp_matrix = new cplx[size * size];
		
		// matrices are dependent on operator not on k value
		// so we map the k value to the operator -> d/dx = - kx / i
		for(int kk = 0; kk < num_k_sets; kk++) {
			for(int ii = 0; ii < 3; ii++) {
				kvalues[kk][ii] = - kvalues[kk][ii] / i;
			} 	
		}
				
		// -----------------------------------------------------	
		// assemble strainless matrix for different k value sets
		// -----------------------------------------------------
		for(int kk = 0; kk < num_k_sets; kk++) {
			// build test string so we may track what failed
			sout.str("");
			sout << "error for " << classname << " (" 
			     << (-i) *kvalues[kk][0] << ", "
			     << (-i) *kvalues[kk][1] << ", "
			     << (-i) *kvalues[kk][2] << ") strainless "; 
			
			// reset temporary matrix
			for(int ii = 0; ii < size*size; ii++) {
				tmp_matrix[ii] = 0.0;	
			}
			// assemble temporary matrix
			for(int ii = 0; ii < 3; ii++) {
				// first order terms
				tmp  = matrix.get_first_order_strainless_matrix(KPMatrixBase::op_left,  ii);
				tmp2 = matrix.get_first_order_strainless_matrix(KPMatrixBase::op_right, ii);
				for(int aa = 0; aa < size; aa++) {
					for(int bb = 0; bb < size; bb++) {
						tmp_matrix[aa * size + bb] += ((*tmp)(aa,bb) + (*tmp2)(aa,bb)) * kvalues[kk][ii];	
					}	
				}
				delete tmp; 
				delete tmp2;			
				// second order terms	
				for(int jj = 0; jj < 3; jj++) {
					tmp = matrix.get_second_order_strainless_matrix(ii,jj);
					for(int aa = 0; aa < size; aa++) {
						for(int bb = 0; bb < size; bb++) {
							tmp_matrix[aa * size + bb] += (*tmp)(aa,bb) * kvalues[kk][ii] * kvalues[kk][jj];	
						}	
					}		
					delete tmp;							
				}				
			} 	
			// zero order terms
			tmp = matrix.get_zero_order_strainless_matrix();
			for(int aa = 0; aa < size; aa++) {
				for(int bb = 0; bb < size; bb++) {
					tmp_matrix[aa * size + bb] += (*tmp)(aa,bb);	
				}	
			}
			delete tmp;
			
			TS_ASSERT(is_hermitian(tmp_matrix, size, sout.str()));
		}
	}			
	

	
	// checks for hermiticity	
	void perform_final_matrix_test(KPMatrixBase& matrix, const char* classname) {		
		
		RMatrix<cplx>* tmp = 0;
		RMatrix<cplx>* tmp2 = 0;	
										
		ostringstream sout;		
				
		cplx kvalues[][3] = {
			{0.0, 0.0, 0.0},			
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0},
			{1.0, 1.0, 0.0},
			{1.0, 0.0, 1.0},
			{0.0, 1.0, 1.0},
			{1.0, 1.0, 1.0}
		};
		
		int   num_k_sets = 8;
		int   size = matrix.get_number_of_bands();
		cplx  i(0.0, 1.0);
		cplx* tmp_matrix = new cplx[size * size];
		
		// matrices are dependent on operator not on k value
		// so we map the k value to the operator -> d/dx = - kx / i
		for(int kk = 0; kk < num_k_sets; kk++) {
			for(int ii = 0; ii < 3; ii++) {
				kvalues[kk][ii] = - kvalues[kk][ii] / i;
			} 	
		}
				
		// -----------------------------------------------------	
		// assemble strainless matrix for different k value sets
		// -----------------------------------------------------
		for(int kk = 0; kk < num_k_sets; kk++) {
			// build test string so we may track what failed
			sout.str("");
			sout << "error for " << classname << " (" 
			     << (-i) *kvalues[kk][0] << ", "
			     << (-i) *kvalues[kk][1] << ", "
			     << (-i) *kvalues[kk][2] << ") final "; 
			
			// reset temporary matrix
			for(int ii = 0; ii < size*size; ii++) {
				tmp_matrix[ii] = 0.0;	
			}
			// assemble temporary matrix
			for(int ii = 0; ii < 3; ii++) {
				// first order terms
				tmp = matrix.get_first_order_matrix(KPMatrixBase::op_left,   ii);
				tmp2 = matrix.get_first_order_matrix(KPMatrixBase::op_right, ii);
				for(int aa = 0; aa < size; aa++) {
					for(int bb = 0; bb < size; bb++) {
						tmp_matrix[aa * size + bb] += ((*tmp)(aa,bb) + (*tmp2)(aa,bb)) * kvalues[kk][ii];	
					}	
				}
				delete tmp; delete tmp2;
				// second order terms	
				for(int jj = 0; jj < 3; jj++) {
					tmp = matrix.get_second_order_matrix(ii,jj);
					for(int aa = 0; aa < size; aa++) {
						for(int bb = 0; bb < size; bb++) {
							tmp_matrix[aa * size + bb] += (*tmp)(aa,bb) * kvalues[kk][ii] * kvalues[kk][jj];	
						}	
					}		
					delete tmp;							
				}				
			} 	
			// zero order terms
			tmp = matrix.get_zero_order_matrix();
			for(int aa = 0; aa < size; aa++) {
				for(int bb = 0; bb < size; bb++) {
					tmp_matrix[aa * size + bb] += (*tmp)(aa,bb);	
				}	
			}
			delete tmp;
			
			TS_ASSERT(is_hermitian(tmp_matrix, size, sout.str()));
		}
	}			


	// checks for hermiticity	
	void perform_straindependent_matrix_test(KPMatrixBase& matrix, const char* classname) {		
		
		RMatrix<cplx>* tmp = 0;	

		ostringstream sout;		
				
		cplx kvalues[][3] = {
			{0.0, 0.0, 0.0},			
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0},
			{1.0, 1.0, 0.0},
			{1.0, 0.0, 1.0},
			{0.0, 1.0, 1.0},
			{1.0, 1.0, 1.0}
		};
		
		int   num_k_sets = 8;

		
		int   size = matrix.get_number_of_bands();
		cplx  i(0.0, 1.0);
		cplx* tmp_matrix = new cplx[size * size];
		
		// matrices are dependent on operator not on k value
		// so we map the k value to the operator -> d/dx = - kx / i
		for(int kk = 0; kk < num_k_sets; kk++) {
			for(int ii = 0; ii < 3; ii++) {
				kvalues[kk][ii] = - kvalues[kk][ii] / i;
			} 	
		}
				
		// -----------------------------------------------------	
		// assemble strain dependent matrix (only zero order left) 
		// -----------------------------------------------------
		// for all strain components
		// reset temporary matrix
		for(int ii = 0; ii < size*size; ii++) {
			tmp_matrix[ii] = 0.0;	
		}		
		// build test string so we may track what failed
		sout.str("");
		sout << "error for " << classname << " zero order strain dependent";		
		for(int ee = 0; ee < 3; ee++) {
			for(int ff = 0; ff < 3; ff++) {	 					
				// zero order terms
				tmp = matrix.get_zero_order_strain_dependent_matrix(ee,ff);
				for(int aa = 0; aa < size; aa++) {
					for(int bb = 0; bb < size; bb++) {
						tmp_matrix[aa * size + bb] += (*tmp)(aa,bb);	
					}	
				}
				delete tmp;			
			}								
		}
		TS_ASSERT(is_hermitian(tmp_matrix, size, sout.str()));
	}		
	
	bool is_hermitian(cplx* matrix, int size, const string& message) const;
	void perform_rotation_test(const string& message, KPMatrixBase& matrix, KPMatrixBase& ref) const;
	//void compare_kp_matrices(const string& message,   const KPMatrixBase& matrix, const KPMatrixBase& ref) const;	
	void assemble_first_order(const Vector3D& dir, vector<complex<double> >& first_order, const KPMatrixBase& ref) const;
	void assemble_second_order(const Vector3D& dir, vector<complex<double> >& second_order, const KPMatrixBase& ref) const;

	void test_kp_matrix_4x4_rotation_gaas() {
		KPMatrix4x4EndersForeman matrix;								
		KPMatrix4x4EndersForeman reference;
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		TS_ASSERT_THROWS_NOTHING(reference.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(reference.calculate());		
		this->perform_rotation_test("kp4x4_gaas", matrix, reference);			
	}	
	
	void test_kp_matrix_6x6_rotation_gaas() {
		KPMatrix6x6EndersForeman matrix;								
		KPMatrix6x6EndersForeman reference;
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		TS_ASSERT_THROWS_NOTHING(reference.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(reference.calculate());		
		this->perform_rotation_test("kp6x6_gaas", matrix, reference);			
	}		

	void test_kp_matrix_8x8_rotation_gaas() {
		KPMatrix8x8EndersForeman matrix;								
		KPMatrix8x8EndersForeman reference;
		TS_ASSERT_THROWS_NOTHING(matrix.set_material(this->material_db->get_material("GaAs")));		
		TS_ASSERT_THROWS_NOTHING(matrix.calculate());	
		TS_ASSERT_THROWS_NOTHING(reference.set_material(this->material_db->get_material("GaAs")));
		TS_ASSERT_THROWS_NOTHING(reference.calculate());		
		this->perform_rotation_test("kp8x8_gaas", matrix, reference);			
	}		
	
};



/** check if kp matrix is hermitian
 */
bool KPMatrixTest::is_hermitian(cplx* matrix, int size, const string& message) const {
	bool hermitian = true;
	ostringstream sout;
	double rl, rr, il, ir, c, d;
	

			
	for(int aa = 0; aa < size; aa++) {
		for(int bb = aa + 1; bb < size; bb++) {
			rl = matrix[aa * size + bb].real();
			rr = matrix[bb * size + aa].real();
			il = matrix[aa * size + bb].imag();
			ir = matrix[bb * size + aa].imag();
			c  = fabs(rl);
			if(c < 1.0e-15) {
				c = 1.0;	
			}
			d = fabs(il);
			if(d < 1.0e-15) {
				d = 1.0;	
			}
			if(fabs((rl - rr) / c) > 1.0e-14 || fabs((il + ir)/d) > 1.0e-14) {
				sout.str("");
				sout << message << " nonhermitian at " << aa << "/" << bb
				     << " " << matrix[aa * size + bb] << " vs. " << matrix[bb * size + aa]; 
				TS_FAIL(sout.str());
				hermitian = false;
			}
		}
	}
	if(!hermitian) {
		cout << " ----------------------------------------------------\n"
		     << message << "\n";
		for(int ii = 0; ii < size; ii++) {
			for(int jj = 0; jj < size; jj++) {
				cout << setw(20) << matrix[ii * size + jj] << "  ";	
			}	
			cout << "\n";
		}		     
		cout << " ----------------------------------------------------\n";     			
	}
		
	return hermitian;
}
	
/** rotate matrix subsequently and compare against reference (unrotated matrix)
 */
void KPMatrixTest::perform_rotation_test(const string& message, KPMatrixBase& matrix, KPMatrixBase& ref) const {
	
	// try to rotate matrix subsequently by some rotations
	// e.g. 4 times 90 deg, 3 times 120 deg around each
	// axis. at the end, it should be equal to the copy
	vector<double> degrees;
	degrees.push_back(30);
	degrees.push_back(45);
	degrees.push_back(60);
	degrees.push_back(90);
	degrees.push_back(120);
	degrees.push_back(180);
	degrees.push_back(33.33333333);
	degrees.push_back(133.2334151);
	

	const int num_bands = matrix.get_number_of_bands();
	
	// get test direction 
	Vector3D        dir;
	RMatrix<double> rotate(3,3);
	ostringstream   sout1, sout2, sout3;

	// initialize matrices iwth values
	matrix.calculate(); 
	ref.calculate();
			
	// create matrices to store the first order terms 
	vector<complex<double> > first_order_matrix(num_bands * num_bands);
	vector<complex<double> > first_order_reference(num_bands * num_bands);	
	// create matrices to store the second order terms
	vector<complex<double> > second_order_matrix(num_bands * num_bands);
	vector<complex<double> > second_order_reference(num_bands * num_bands);
	
	try {

		// do for all degrees
		for(unsigned dd = 0; dd < degrees.size(); dd++) {

			// message string for TS 
			sout1.str("");
			sout1 << message << "_x_rot_by_" << degrees[dd] << "_deg";
			
			// do for different rotation directions
			for(unsigned int rr = 0; rr < 3; rr++) {
				
				// message for TS
				sout2.str(""); sout2 << sout1.str() << "_rot_" << rr;
										
				// set rotation dir
				Vector3D rot_dir;
				switch(rr) {
					case 0:
						rot_dir = Vector3D(1.0, 0.0, 0.0);
						break;
					case 1:
						rot_dir = Vector3D(0.0, 1.0, 0.0);
						break;
					case 2:
						rot_dir = Vector3D(0.0, 0.0, 1.0);
						break;						
				}
			
				// get rotation matrix
				rotate = RMatrix<double>::get_rotation_matrix(rot_dir, cos(2.0*constants::pi * degrees[dd] / 360));								
										
				// set rotation					
				matrix.set_rotation(rotate);
				matrix.calculate();
				
				// loop over different (random) directions
				for(unsigned vv = 0; vv < 10; vv++) {
					//
					switch(vv) {
						case 0:
							dir = ex;
							break;
						case 1:
							dir = ey;
							break;
						case 2:
							dir = ez;
							break;
						default:
							dir = ValueChecker<double>::get_random_dir();											
					}
					// msg TS
					sout3.str(""); sout3 << sout2.str();
					sout3 << dir;	
				
					// assemble first order 
					// first_order_matrix and first_order_reference should now be equal
					// (as we transformed dir to the same coordinate system)
					assemble_first_order(rotate * dir, first_order_matrix, matrix); 
					assemble_first_order(         dir, first_order_reference, ref);
				
					// compare first order
					try { 
						ValueChecker<complex<double> >::compare_vectors(first_order_matrix, first_order_reference, 1.0e-8);
					} catch (ValueCheckerError& e) {
						ostringstream smsg;
						smsg << sout3.str() << " first order failed for " << dir;
						TS_FAIL(smsg.str());
						smsg << "\n" <<  e.get_message();
						Logger::get_instance()->emit(LOG_INFO, smsg.str()); 
						TDKP_GENERAL_EXCEPTION(smsg.str());				     	
					}
		
					
					// assemble second order 
					assemble_second_order(rotate * dir, second_order_matrix, matrix); 
					assemble_second_order(         dir, second_order_reference, ref);			
					try { 
						ValueChecker<complex<double> >::compare_vectors(second_order_matrix, second_order_reference, 1.0e-6);
					} catch (ValueCheckerError& e) {
						ostringstream smsg;
						smsg << sout3.str() << "second order failed for " << dir;
						TS_FAIL(smsg.str());
						smsg << "\n" <<  e.get_message();
						Logger::get_instance()->emit(LOG_INFO, smsg.str());
						TDKP_GENERAL_EXCEPTION(smsg.str()); 				     	
					}			
				}
			}			
		}
	} catch(Exception* e) {
		sout1.str("");
		sout1 << "rotation test failed for " << message 
		     << ". WILL NOT CONTINUE. check log file\n";
		TS_FAIL(sout1.str());
		TS_WARN("kp matrices dumped to files");
		sout1.str("");
		sout1 << message << "_rotated.dat";
		matrix.dump(sout1.str().c_str());
		sout1.str("");
		sout1 << message << "_reference.dat";
		ref.dump(sout1.str().c_str());						
	}
		 		
}

void KPMatrixTest::assemble_first_order(const Vector3D& dir, vector<complex<double> >& first_order, const KPMatrixBase& ref) const {
	
	const int num_bands = ref.get_number_of_bands();
	int sparsity_num;
	ref.get_sparsity_pattern(sparsity_num);
	first_order.assign(num_bands * num_bands, 0.0);
	
	for(int ii = 0; ii < 3; ii++) {
		for(short gg = 0; gg < 2; gg++) {
			const vector<complex<double> >& first = ref.get_first_order(gg, ii);
			for(int ss = 0; ss < sparsity_num; ss++) {
				first_order[ss] += dir(ii) * first[ss];	
			}
		}

	}
		
}

void KPMatrixTest::assemble_second_order(const Vector3D& dir, vector<complex<double> >& second_order, const KPMatrixBase& ref) const {
	
	const int num_bands = ref.get_number_of_bands();
	int sparsity_num;
	ref.get_sparsity_pattern(sparsity_num);
	second_order.assign(num_bands * num_bands, 0.0);
	
	for(int ii = 0; ii < 3; ii++) {
		for(int jj = 0; jj < 3; jj++) {
			const vector<complex<double> >& second = ref.get_second_order(ii,jj);
			for(int ss = 0; ss < sparsity_num; ss++) {
				second_order[ss] += dir(ii) * dir(jj) * second[ss];	
			}
		}
	}
	
}




		
