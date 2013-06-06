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

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_


#include <string>
#include <map>
#include "tdkp/common/all.h"
#include "tdkp/main/PropertyContainer.h"

#ifdef isset
#undef isset
#endif

namespace tdkp {

/** basic program configuration class for main tdkp operation
 * 
 * holds program parameters as doubles
 * should only be used for debugging and tuning and not for daily work
 * forwards subproblem (e.g. eigensolver) configuration classes
 */
class Configuration {

public:

	static Configuration* get_instance(); 
	
	// -------------------------------------------------------
	// individual property container access (e.g. eigensolver)
	// -------------------------------------------------------
	const  PropertyContainer<double>* get_config(const char* filename) throw(Exception*);	
	const  PropertyContainer<double>* get_config(const std::string& filename) throw(Exception*);

	double get_config_value(const char* filename, const char* key) throw(Exception*);
	double get_config_value(const std::string& filename, const std::string& key) throw(Exception*);
	
	// ----------------------------------------------------------------
	// tdkp property access (controlled and validated via controls.cnf)
	// ----------------------------------------------------------------
	void   load (const char*   key);
	void   set  (const char*   key, double value) throw(Exception*);
	void   set  (const string& key, double value) throw(Exception*);	
	bool   isset(const char*   key) const throw(Exception*);
	bool   isset(const string& key) const throw(Exception*);	
	double get  (const char*   key) const throw(Exception*);
	double get  (const string& key) const throw(Exception*);
	
	
protected:
	std::string add_conf_path_to_filename(const char* filename) throw(Exception*);
	std::string add_conf_path_to_filename(const std::string& filename) throw(Exception*);							
		
	map<string,PropertyContainer<double>*> cache;
	typedef map<string,PropertyContainer<double>*>::iterator cache_iterator;
	
	/** general program controls */
	PropertyContainer<double> controls;
	
private:	
	Configuration();
	~Configuration();
	static Configuration* singleton;	
	void test_capabilities();
};

}

#endif /*CONFIGURATION_H_*/
