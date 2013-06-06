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


#include <stdio.h>
#include "tdkp/geometry/Region.h"

namespace tdkp {


/** create region 
 * 
 * @param name_ name of region
 */	
Region::Region(const char* name_) : name(name_), num_elements(0), is_enabled(true) {
	this->mat = 0;		
}
Region::Region(const string& name_) : name(name_), num_elements(0), is_enabled(true) {
	this->mat = 0;
}
Region::~Region() {
	this->mat = 0;	
}
Region::Region() 
: name(""), material_name(""), mat(0), num_elements(0), is_enabled(true) {
	
}
				
//const char* Region::get_name() const {
//	return this->name.c_str();	
//}

const string& Region::get_name() const {
	return this->name;	
}

/** set material name
 * 
 * during grid file reading, only the material name will be set to the region. later, 
 * defined in ProblemDefinition, the material with that name is loaded from the material
 * database and then assigned to the region
 */ 
void Region::set_material_name(const char* key) {
	this->material_name = key;		
}
const string& Region::get_material_name() const {
	return this->material_name;
}

void Region::set_material(const Material* mat_) {
	this->mat = mat_;	
}

const Material& Region::get_material() const {
	TDKP_ASSERT(this->mat != 0, "");
	return *this->mat;	
}

/** number of elements assigned to that region
 * 
 * the region has now information about the elements (and does not need them) but
 * the elements have a pointer to their region. for some df-ise file writings, we
 * need the number of elements in advance. 
 */
void Region::set_num_elements(unsigned int nelem) {
	this->num_elements = nelem;	
}
unsigned int Region::get_num_elements() const {
	return this->num_elements;	
}

} // end of namespace tdkp
