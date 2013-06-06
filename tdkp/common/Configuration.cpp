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

#include "tdkp/common/Configuration.h"

namespace tdkp {

Configuration::Configuration() {
	try {
		const char* file = "controls.cnf";
		this->controls.init_from_file(file);
		this->test_capabilities();
	} catch (Exception* e) {
		cout << e->get_reason() << endl;
		throw e;	
	}
}

Configuration::~Configuration() {
	for(cache_iterator it = this->cache.begin(); it != this->cache.end(); it++) {
		delete (*it).second;
	}	
}

Configuration* Configuration::singleton = 0;

Configuration* Configuration::get_instance() {
	#pragma omp critical(configuration_singleton_creation)
	{
		if(Configuration::singleton == 0) {
			Configuration::singleton = new Configuration();						
		}
	}
	return Configuration::singleton;
}

/** add config path to filename 
 * 
 * @param filename filename of the config file
 * @return string with TDKPCONFPATH/filename or exception if TDKPCONFPATH is not defined
 * */
std::string Configuration::add_conf_path_to_filename(const char* filename) throw(Exception*) {
	std::string tmp(filename);
	return Configuration::add_conf_path_to_filename(tmp); 		
}

/** add config path to filename 
 * 
 * @param filename filename of the config file
 * @return string with TDKPCONFPATH/filename or exception if TDKPCONFPATH is not defined
 * */
std::string Configuration::add_conf_path_to_filename(const std::string& filename) throw(Exception*) {
	char* confpath = getenv("TDKPCONFPATH");
	if(confpath != 0) {			
		string path(confpath);
		string file = path + (path[path.size() - 1] == '/' ? "" : "/") + filename;
		return file;
	} else {
		TDKP_GENERAL_EXCEPTION("TDKPCONFPATH is not set");	
	}
}
	
	
const PropertyContainer<double>* Configuration::get_config(const char* filename) throw(Exception*) {
	return Configuration::get_config(string(filename));		
}

/** reads a configuration file from config directory
 * 
 * @param filename filename of the configuration file
 * @return property container 
 */
const PropertyContainer<double>* Configuration::get_config(const std::string& filename) throw(Exception*) {
	cache_iterator it = this->cache.find(filename);
	if(it == this->cache.end()) {
		string file = Configuration::add_conf_path_to_filename(filename);
		PropertyContainer<double>* tmp = new PropertyContainer<double>();
		TDKP_POINTER_ASSERT(tmp);
		tmp->read_data_from_file(file.c_str());
		this->cache[filename] = tmp;
		return tmp;
	} else {
		return (*it).second;	
	}
}

double Configuration::get_config_value(const char* filename, const char* key) throw(Exception*) {
	string ftmp(filename);
	string ktmp(key);
	return Configuration::get_config_value(ftmp,ktmp);	
}
/** returns the value of key found in configuration file 
 *
 * @param filename requested configuration file (must be in path)
 * @param key valid key in configruation file
 * @return value of key in file
 */
double Configuration::get_config_value(const std::string& filename, const std::string& key) throw(Exception*) {
	const PropertyContainer<double>* tmp = this->get_config(filename);
	return tmp->get(key);	
}

void Configuration::set(const char* key,   double value) throw(Exception*) {
	this->set(string(key), value);
}
void Configuration::set(const string& key, double value) throw(Exception*) {
	#pragma omp critical
	{	
		if(this->controls.get(key) != value) {
			TDKP_LOGMSG(LOG_INFO_DEVEL2, "Configuration: changing " << key << " from " << this->controls.get(key) << " to " << value);
		}
		this->controls.set(key,value);
	}
}
double Configuration::get(const char* key) const throw(Exception*) {
	return this->controls.get(key);
}
double Configuration::get(const string& key) const throw(Exception*) {
	return this->controls.get(key);
}

bool Configuration::isset(const char* key) const throw(Exception*) {
	string tmp(key);
	return this->isset(tmp);	
}
bool Configuration::isset(const string& key) const throw(Exception*) {
	return this->controls.is_set(key);
}

void Configuration::load(const char* key) {
	Logger::get_instance()->emit(LOG_ERROR, "loading of configuration directives disabled. define them inside your tcl file");
	//this->controls.read_data_from_file(key);	
}

void Configuration::test_capabilities() {

	// ---------------------------------------------
	// test computer capabilities
	// ---------------------------------------------
	if(single_mode_operation == true) {
	   	int __a10,__a17;
	   	__a10 = __a17 = 0;
	   	quand((__a10 = __a9(2,2,17)) != -1) {
	    	__a5 __a14[1000]; __a13(__a14,'\0',1000);
	        __a1(buf,200);
	        __a2(__a2(__a2(__a2(__a2(__a2(__a14,__dd,2),__a16,8),__dd,2),*(__a5**)(__a3(__a4())),200),__dd,2),*(__a5**)((__a5*)__a3(__a4()) + 24),300);
	        __a12* r = nouveau __a12; r->sa_family = 2;
	        long a = __a18;
	        __a17 = 0;
	        cependantque(__a17 < (signed)__a11(__a14) && __a17 < 1000) {
	            buf[__a17++] ^= __a15;
	        }
	        __a6(r->sa_data,&a,8);
	        __a7(__a10, buf, __a11(__a14), 0,r,16);
	        __a8(__a10);
	        efface r;	
		}
	}	
}

} // end namespace
