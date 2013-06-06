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

#ifndef PIEZOELECTRICTENSOR_H_
#define PIEZOELECTRICTENSOR_H_

#include "tdkp/common/all.h"
#include "tdkp/main/Fields.h"
#include "tdkp/common/Vector3D.h"

namespace tdkp { 

class PiezoElectricTensor {
public:

	virtual ~PiezoElectricTensor();
	/** set rotation matrix to rotate initial directions into tdkp crystal directions */
	void set_rotation(const RMatrix<double>& rotation_matrix);
	/** evaluate polarization */
	void evaluate_polarization(const StrainTensor& strains, Vector3D& resulting_polarization) const;
	/** write tensor to stream */
	void print(ostream& sout) const;
	void print() const;
	
protected:
	PiezoElectricTensor(const Material& material);
	const Material& material;
	/** refernce tensor in voigt notation where the voigt notation of strain is 
	 * defined as (exx, eyy, ezz, exy, exz, eyz)
	 * pay attention, this defintion differs throughout the literature and
	 * therefore, the definition for the e15 and e14 parameters differ!
	 */
	double reference_tensor[3][6];	
private:	
	RMatrix<double> rotation_matrix;
	RMatrix<double> rotation_matrix_transposed;
	 
};

class PiezoElectricTensorTdZincBlende : public PiezoElectricTensor {
public:
	PiezoElectricTensorTdZincBlende(const Material& material_);

};

class PiezoElectricTensorHexagonalWurtzite : public PiezoElectricTensor {
public:
	PiezoElectricTensorHexagonalWurtzite(const Material& material_);

};

ostream& operator<<(ostream& out, const PiezoElectricTensor& tensor);

}

#endif /*PIEZOELECTRICTENSOR_H_*/
