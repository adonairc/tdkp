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

#include "tdkp/utilities/PiezoElectricTensor.h"

namespace tdkp
{

ostream& operator<<(ostream& out, const PiezoElectricTensor& tensor) {
	tensor.print(out);
	return out;	
}

void PiezoElectricTensor::print(ostream& sout) const {
	// print reference tensor in voigt notation
	sout.precision(5);
	for(unsigned int ii = 0; ii < 3; ii++) {
		for(unsigned int jj = 0; jj < 6; jj++) {
			sout << setw(8) << reference_tensor[ii][jj] << " ";			
		}
		sout << "\n"; 	
	}		
}

void PiezoElectricTensor::print() const {
	ostringstream sout;
	print(sout);
	TDKP_LOGMSG(LOG_INFO, sout.str());	
}

PiezoElectricTensor::PiezoElectricTensor(const Material& material_)
: material(material_) 
{
	// ----------------------------------------------
	// initially we have the identity matrix as rotation
	// ----------------------------------------------
	this->rotation_matrix            = identity_matrix;
	this->rotation_matrix_transposed = identity_matrix;
	
	// ----------------------------------------------
	// init tensor to piezotensor to zero
	// ----------------------------------------------
	for(unsigned int ii = 0; ii < 3; ii++) {
		for(unsigned int jj = 0; jj < 6; jj++) {
			reference_tensor[ii][jj] = 0.0;	
		}	
	}	
}
PiezoElectricTensor::~PiezoElectricTensor() {}
/** set rotation matrix to rotate initial directions into tdkp crystal directions */
void PiezoElectricTensor::set_rotation(const RMatrix<double>& rotation_matrix_) {
	this->rotation_matrix = rotation_matrix_;
	this->rotation_matrix_transposed = rotation_matrix.get_transpose(); 
}
/** evaluate polarization 
 * 
 * strain should be given in system coordinates k' and the rotation matrix maps
 * the principal kp crystal axes to the system axes via k' = Rk, so the strain is 
 * backmapped to the principal axes system via E = R^t E' R
 * 
 * the piezo polarization is given by P = TE and is then backmapped to the system axes
 * via P' = RP 
 */
void PiezoElectricTensor::evaluate_polarization(const StrainTensor& strains, Vector3D& resulting_polarization) const {
		
	// ---------------------------------------
	// init polarization to 0
	// ---------------------------------------
	resulting_polarization = Vector3D(0.0, 0.0, 0.0);
	
	// ---------------------------------------
	// rotate strain tensor into reference coordinates
	// E = R^t E R
	// ---------------------------------------
	StrainTensor reference_strain(strains);
	reference_strain.rotate(rotation_matrix_transposed);
	
	// ---------------------------------------
	// evaluate polarization
	// ---------------------------------------
	for(unsigned int ii = 0; ii < 3; ii++) {
		for(short jj = 0; jj < 6; jj++) {
			resulting_polarization(ii) += reference_tensor[ii][jj] * reference_strain.get(jj);	
		}	
	}
	
	// ---------------------------------------
	// rotate polarization into system coordinates
	// P' = RP
	// ---------------------------------------
	resulting_polarization = rotation_matrix * resulting_polarization;
		
}

PiezoElectricTensorTdZincBlende::PiezoElectricTensorTdZincBlende(const Material& material_)
: PiezoElectricTensor(material_)
{
	// -------------------------------------------
	// attention, the factor of 2 comes from the problem
	// with the voigt notation as people sometimes use voigt
	// notation including the factor of two in the shear strains.
	// some people (e.g. vurgaftman) do that different etc. etc.
	// so it's a huge mess, but i stick to the definition that 
	// the voigt notation used here has no factor of two and that the factor 
	// of two is not condensed into the coefficent of e14
	// therefore i have to add it explicitly.
	// -------------------------------------------
	const double& e14 = material.get("piezo_coefficient_e14");
	reference_tensor[0][5] = 2.0 * e14; // Px = 2.0 * e14 * eyz 
	reference_tensor[1][4] = 2.0 * e14; // Py = 2.0 * e14 * exz
	reference_tensor[2][3] = 2.0 * e14; // Pz = 2.0 * e14 * exy
}

PiezoElectricTensorHexagonalWurtzite::PiezoElectricTensorHexagonalWurtzite(const Material& material_)
: PiezoElectricTensor(material_)
{
	const double& e31 = material.get("piezo_coefficient_e31");
	const double& e33 = material.get("piezo_coefficient_e33");
	const double& e15 = material.get("piezo_coefficient_e15");
	reference_tensor[0][4] = 2.0 * e15; // Px = 2 e15 exz
	reference_tensor[1][5] = 2.0 * e15; // Py = 2 e15 eyz
	reference_tensor[2][0] = e31; // Pz = e31 exx + e31 eyy + e33 ezz
	reference_tensor[2][1] = e31; // Pz = e31 exx + e31 eyy + e33 ezz
	reference_tensor[2][2] = e33; // Pz = e31 exx + e31 eyy + e33 ezz 
}


}
