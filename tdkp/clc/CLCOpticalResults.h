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

#ifndef CLCOPTICALRESULTS_H_
#define CLCOPTICALRESULTS_H_

#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"
#include "tdkp/clc/CLCIngredients.h"

namespace tdkp {
	
/** result class for the optical suszeptibility calculation */
class CLCOpticalResults : public XYData<double> {
public:
	CLCOpticalResults(
		const CLCIngredients& ingredients,
		const double& omega_min,
		const double& omega_max,
		const vector<cplx>& suszeptibility,
		const vector<double>& spont_emission,
		const double& cb_fermi,
		const double& vb_fermi,
		const double& temperature
	);		
	const vector<double>& get_absorption() const;
	const vector<double>& get_spont_emission() const;
	const vector<double>& get_spont_emission_kms() const;
	const vector<double>& get_delta_n() const;
	const vector<cplx>&   get_suszeptibility() const;
	
	const double& get_omega_min() const;
	const double& get_omega_max() const;
	unsigned int  get_omega_num() const;
	
	// -----------------------------------------------
	// interface XYData 
	// -----------------------------------------------
	virtual int         get_x_length()                const;
	virtual int         get_num_y_sets()              const;
	virtual int         get_num_x_sets()              const;
	virtual void        get_x(int xidx, std::vector<double> &x) const;
	virtual void        get_y(int yidx, std::vector<double>& y) const; 
	virtual std::string get_x_identifier(int xidx)   const;
	virtual std::string get_y_identifier(int yidx) 	 const;		
	
private:
	void evaluate_spont_emission_kms();

	vector<double> absorption;
	vector<double> spont_emission;
	vector<double> spont_emission_kms;
	vector<double> delta_n;
	vector<cplx>   suszeptibility;
	
	double omega_min;
	double omega_max;
	
	const CLCIngredients& ingredients;
	
	double cb_fermi;
	double vb_fermi;
	double temperature;
			
};
	
	
}

#endif /*CLCOPTICALRESULTS_H_*/
