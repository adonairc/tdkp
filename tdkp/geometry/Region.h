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

#ifndef _REGION_H
#define _REGION_H

#include <string>
#include "tdkp/main/MaterialDatabase.h"

namespace tdkp {

	class Region {
	public:
		Region();
		Region(const char* name);
		Region(const string& name);
		~Region();
					
		void  set_enabled(bool is_enabled_) { this->is_enabled = is_enabled_; }
		bool  enabled() const { return this->is_enabled; }
		const string& get_name() const;				
		void  set_material_name(const char* name);
		const string& get_material_name() const;
		void  set_material(const Material* mat);
		const Material& get_material() const;
		bool  material_set() const { return this->mat != 0; }
		void set_num_elements(unsigned int nelem);
		unsigned int get_num_elements() const;
		
		int  get_index_global()   const { return index_global; }
		void set_index_global(unsigned int index_global_) { index_global = index_global_; }		
						
	private:
		int index_global;
		string name;
		string material_name;
		const Material* mat;					
		unsigned int num_elements;
		bool is_enabled;
		bool is_pml;			
	};

}		


#endif
