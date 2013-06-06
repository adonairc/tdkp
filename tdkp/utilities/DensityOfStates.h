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

#ifndef DENSITYOFSTATES_H_
#define DENSITYOFSTATES_H_

#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/common/ICurves.h"
#include "tdkp/main/Bandstructure.h"

namespace tdkp {

/** calculate density of states using the radial approximation 
 * 
 * we do NOT include spin degeneracy in the DOS! (as kp does also lead
 * to two degenerate solutions) 
 */
class DensityOfStates : public XYData<double> {
public:
	DensityOfStates();
	DensityOfStates(unsigned int dim);
	virtual ~DensityOfStates();	
	
	// --------------------------------------------
	// initialization and assignment
	// --------------------------------------------
	void set_dimension(unsigned int dim);
	void set_bandstructure(
		unsigned int number_of_subbands, 
		const BandstructureDomain<complex<double> >& bands
	);		
	void set_bandstructure(
		const double& kmin, 
		const double& kmax, 
		unsigned int num_k_values, 
		unsigned int number_of_subbands, 
		const vector<double>& bandstructure
	);
	void set_bandstructure(
		const double& kmin, 
		const double& kmax, 
		unsigned int num_k_values, 
		unsigned int number_of_subbands, 
		const vector<cplx>& bandstructure
	);	
	void calculate();  
			
	// ---------------------------------------------
	// result access
	// ---------------------------------------------
	const ICurve& get_dos(unsigned int subband) const;
	const ICurve& get_full_dos() const;			
			
	// ---------------------------------------------
	// XY data interface
	// ---------------------------------------------						
	int    get_x_length()                const;
	int    get_num_y_sets()              const;
	int    get_num_x_sets()              const;
	void   get_x(int xidx, vector<double> &x) const;
	void   get_y(int yidx, vector<double>& y) const; 
	string get_x_identifier(int xidx)    const;
	string get_y_identifier(int yidx) 	 const;	
	
private:
	unsigned int         dimension;
	vector<double> 		 k_values;
	vector<double> 		 bandstructure;
	unsigned int 		 number_of_subbands;
	vector<LinearCurve>  band_density_of_states;
	LinearCurve			 total_density_of_states;
	
	 
	
	class DOSSplineCurve : public SplineCurve {
	public:		
		DOSSplineCurve(BoundaryCondition abc);
  		DOSSplineCurve(const ICurve& curve);
  		DOSSplineCurve(const GridPoints& x, const GridValues& y);
  		  	
		void prepare_segments();	  		  		
		void find_x_values(const double& y0_value, vector<double>& xi_values) const;
		
		void get_minmax_all(double& miny, double& maxy);
		  		
	private:
		struct MinMaxValue {
			double min_value;
			double max_value;			
		};
		vector<MinMaxValue> segment_ranges;
			
	};
	
	void merge_sort_and_unique(const vector<double>& one, const vector<double>& two, vector<double>& target) const;
	
};

}

#endif /*DENSITYOFSTATES_H_*/
