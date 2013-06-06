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

#ifndef PROPERTYCONTAINER_H_
#define PROPERTYCONTAINER_H_


#include <map>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <algorithm>
#include "tdkp/common/all.h"
#include "tdkp/common/Logger.h"


using namespace std;
using boost::lexical_cast;
using boost::bad_lexical_cast;

namespace tdkp {
	
class FileException;

/** general property container class usable for compuational setup, material etc. */
template<class T>
class PropertyContainer
{
public:
	
	typedef typename map<string, T>::const_iterator map_const_iterator;
	typedef typename map<string, T>::iterator 	    map_iterator;

	PropertyContainer();
	PropertyContainer(const char* conffile);
	PropertyContainer(const PropertyContainer&);
	
	virtual ~PropertyContainer();	
	void 		  init_from_file(const char* filename) throw(Exception*);			 
	void          read_data_from_file(const char* filename) throw(Exception*);
	void          write_data_to_file(const char* filename) const;
	void          clear();
	bool 		  valid();
  	virtual bool  valid_value(const char* key, T value) throw(Exception*); 
	virtual bool  valid_value(const string &key, T value) throw(Exception*); 
	bool  		  valid_key(const char* key)   const;
	bool  		  valid_key(const string& key) const;
	bool  		  is_set(const char* key) const throw(Exception*);
	bool  		  is_set(const string& key) const throw(Exception*);	
	T  	  		  get(const char* key) const throw(Exception*);
	T 	  		  get(const string& key) const throw(Exception*);
	void  		  set(const char* key, T value) throw(Exception*);
	void  		  set(const string& key, T value) throw(Exception*);	
	unsigned int  get_id() const { return this->id; }
	void		  set_id(unsigned int id_);	
			
	const string&
	get_configuration_filename() const { return this->configuration_file; }			
					    		
protected:
	/** property definition container */
	struct PropDef {
		string key;
		T error_min, error_max, warn_min, warn_max, default_value;	
		bool mandatory, default_value_avl;		
		bool inside_error_bounds(T value) {
			return value >= error_min && value <= error_max;			
		}		
		bool inside_warn_bounds(T value) {
			return value >= warn_min && value <= warn_max;	
		}		
		void clear() {
			key.clear(); error_min = error_max = warn_min = warn_max = default_value = (T)0;
			default_value_avl = false;	
		}
	};

	map<string, T> value_map;
	map<string, PropDef> valid_properties;
	bool   initialized;
	string configuration_file;
	int    id; 
	
	
	void check_proper_initialization() const;		
	void set_default_values();	
	void strip_whitespaces(string& str);
			
	/** write formatted property values to stream */
	friend std::ostream& operator<<(std::ostream& stream, const PropertyContainer<T> &prop) {
		typename PropertyContainer<T>::map_const_iterator it;
		unsigned int string_length = 0;
		// get max length of our key
		for(it = prop.value_map.begin(); it != prop.value_map.end(); it++) {
			if((*it).first.size() > string_length) {
				string_length = (*it).first.size();	
			}						
		}			
		// write to stream
		stream.setf(ios::left);
		for(it = prop.value_map.begin(); it != prop.value_map.end(); it++) {	
			
			stream << std::setw(string_length + 3) << (*it).first << " " 
				   << (*it).second << "\n";			   
		}
		return stream;					
	}
};

/** no pointers ... no cleanup */
template<class T>
PropertyContainer<T>::~PropertyContainer() {

}

/** create noninitialized property container 
 * 
 * property container value files can also point to the corresponding definition
 * file. therefore you can define a uninitialized property container and try to
 * load a value file and get the corresponding definition from that file. 
 * see read_data_from_file. 
 * */
template<class T>
PropertyContainer<T>::PropertyContainer() : initialized(false) {
	
}

/** create property container and read configuration from file */
template<class T>
PropertyContainer<T>::PropertyContainer(const char* filename) {
	this->init_from_file(filename);	
}

template<class T>
PropertyContainer<T>::PropertyContainer(const PropertyContainer& copy) {
	this->value_map.insert(copy.value_map.begin(), copy.value_map.end());
	this->valid_properties.insert(copy.valid_properties.begin(), copy.valid_properties.end());		
	this->initialized = copy.initialized;
}
	
template<class T>	
void PropertyContainer<T>::strip_whitespaces(string& str) {
	string tmp;
	for(string::iterator it = str.begin(); it != str.end(); it++) {
		if((*it) != ' ' && (*it) != '\t') {
			tmp.push_back(*it);				
		}			
	}
	str = tmp;
}
	
/** read property class setup from file 
 * 
 * instead of hard coding the required values and possible options, they can be
 * easily read from file. the file format must have the following form:
 * # this is a comment. the # must be the first char. 
 * keyname;(mandatory|optional);error_min;error_max;warn_min;warn_max;(default(if optional is set))
 * keyname is the property key
 * (mandatory|optional) denotes whether the specifiy property must be set or could be mandatory
 * error_(min|max) value range where an exception is thrown if such a value is set
 * warn_(min|max)  value range where only a warning message is written to logmsg, but value is accepted
 * default		   if value is optional, then a default value may be set
 * 
 * @param filename property class config filename
 * */
 template<class T>
void PropertyContainer<T>::init_from_file(const char* filename) throw(Exception*) {
	ifstream fin;
	ostringstream msg;
	string str_in;
	string str_stripped;
	string item;
	string::size_type loc;
	string::size_type last_loc;
	string::iterator it;
	PropDef prop;
	int pcount;
	int line = 1;	
		
	
	// -------------------------------------------------------------
	// open config file
	// -------------------------------------------------------------
	fin.open(filename);
	
	if(!fin) {
		fin.clear();
		// did not found file in current path. looking for TDKPCONFPATH
		char* confpath = getenv("TDKPCONFPATH");
		if(confpath != 0) {			
			string path(confpath);
			string file = path + (path[path.size() - 1] == '/' ? "" : "/") + filename;
			fin.open(file.c_str());
			if(!fin) {			
				ostringstream sout;
				sout << "can not open or find property definition file \"" << file << "\". "
				     << "set TDKPCONFPATH to point to the configuration directory "; 	
				TDKP_GENERAL_EXCEPTION(sout.str()); 
			}
			msg << "PropertyContainer: reading property container configuration from file " << file;				
			Logger::get_instance()->emit(LOG_INFO_DEVEL2, msg.str());	
		} else {			
			TDKP_GENERAL_EXCEPTION("can not open or find property definition file " << filename << ". set TDKPCONFPATH to point to the configuration directory");		
		}
	} else {
		msg << "PropertyContainer: reading property container configuration from file " << filename;				
		Logger::get_instance()->emit(LOG_INFO_DEVEL2, msg.str());			
	}
		
	// ------------------------------------------------------------
	// read config file
	// ------------------------------------------------------------	
	while(getline(fin,str_in)) {
		// ignore comments
		if(str_in.size() > 0 && str_in[0] != '#') {
			str_stripped.clear();
			// strip whitespaces
			this->strip_whitespaces(str_stripped = str_in);
			if(str_stripped.size() > 0) {
				prop.clear();				
				// parse line
				pcount   = 0;
				last_loc = 0;
				while((loc = str_stripped.find(";",last_loc)) != string::npos) {	
					item = str_stripped.substr(last_loc, (loc - last_loc));								
					try {
						// ------------------------------------------------------------
						// cast into PropDef
						// ------------------------------------------------------------					
						switch(pcount) {
							case 0:
								prop.key = item;
								break;
							case 1:
								if(item == "optional") {
									prop.mandatory = false;								
								} else if(item == "mandatory") {
									prop.mandatory = true;	
								} else {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: expected keyword \"mandatory\" or \"optional\"", filename, line);
								}
								break;						
							case 2:					
								prop.error_min = lexical_cast<T>(item);
								break;
							case 3:
								prop.error_max = lexical_cast<T>(item);
								// check if min > max
								if(prop.error_min >= prop.error_max) {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: error_min >= error_max ", filename, line);
								}
								break;
							case 4:
								prop.warn_min = lexical_cast<T>(item);
								if(prop.warn_min < prop.error_min) {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: warn_min < error_min ", filename, line);
								}
								break;
							case 5:
								prop.warn_max = lexical_cast<T>(item);
								if(prop.warn_min >= prop.warn_max) {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: warn_min >= warn_max ", filename, line);								
								}												
								if(prop.warn_max > prop.error_max) {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: warn_max > error_max ", filename, line);
								}							
								break;
							case 6:
								if(prop.mandatory) {
									throw new Exception(1000 + strlen(filename), "in file %s on line %d: property is mandatory, therefore default value must not be set", filename, line);	
								}
								prop.default_value     = lexical_cast<T>(item);							
								prop.default_value_avl = true;
								break;
						}	
					} catch(boost::bad_lexical_cast &) {
						throw new Exception(1000 + item.size() + strlen(filename), "in file %s on line %d: could not parse argument nr. %d: %s", filename, line, pcount + 1, item.c_str());						
					}
					pcount++;	
					last_loc = loc + 1;	
				}  // end while(loc .. )	
				if(pcount < 6) {
					throw new Exception(1000 + strlen(filename), "in file %s on line %d: expected at least 5 arguments, found %d", filename, line, pcount - 1);	
				}
				// ------------------------------------------------------------
				// add to valid properties map
				// ------------------------------------------------------------
				if(this->valid_properties.find(prop.key) == this->valid_properties.end()) {		
					this->valid_properties[prop.key] = prop;
				} else {
					throw new Exception(1000 + strlen(filename) + prop.key.size(), "in file %s on line %d: property %s already defined", filename, line, prop.key.c_str());
				}
			} // end if(str_stripped.size() > 0)
		} // end if(str[0] != '#')	
		line++;			
	} // end while(getline)	
	this->initialized 		 = true;
	this->configuration_file = string(filename);
	this->set_default_values();
	fin.close();	
} 
// ------------------------------------------------------------
// end init_from_file
// ------------------------------------------------------------

template<class T>
void PropertyContainer<T>::check_proper_initialization() const {
	if(!this->initialized) {
		TDKP_GENERAL_EXCEPTION("PropertyContainer<T> was not properly initilalized");	
	}	
}


/** read property containers data from given file
 * 
 * the file syntax is simple:
 * # comment (line must start with #)
 * <keyword> = value
 * only keywords that were set in container definiton file are allowed
 * unknown keywords generate exceptions
 * 
 * container objects can also be initialized upon reading a data file.  
 * in this case, the data file must start with a line defining the path to 
 * the definition file
 * # tdkpcondef:  <path container definition file>
 * 
 * if the container object already has been initialized and the definition files 
 * differ, the data file is read according to the loaded and not according to 
 * the referenced configuration file. also, a warning will be written to the log
 * 
 */
template<class T>
void PropertyContainer<T>::read_data_from_file(const char* filename) throw(Exception*) {					
	ifstream fin;
	fin.open(filename);
	if(!fin) {
		TDKP_GENERAL_EXCEPTION("can not open data file " << filename);		
	} 		
	string str_in;
	string str_stripped;
	string key;
	T      value;
	string::size_type loc;
	string::iterator it;	
	list<string> keys_in_file;
	ostringstream msg;
	int line = 1;
			
	msg << "PropertyContainer: reading data from file " << filename;				
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, msg.str());
			
	// ------------------------------------------------------------
	// read definition file
	// ------------------------------------------------------------	
	while(getline(fin,str_in)) {
		// --------------------------------------------------------------------------
		// if we're on the first line, we check for any referenced configuration file
		// --------------------------------------------------------------------------
		if(line == 1) {
			//  check if confdef exists 
			if((loc = str_in.find("# tdkpconfdef:", 0)) != string::npos) {
				string referenced_conffile = str_in.substr(loc + 14);
				this->strip_whitespaces(referenced_conffile);				
				// if we are already initialized, we just check if its the same file and emit a warning
				// if its not the same file. if we are not initialized, we try to initialize us
				if(this->initialized) {					
					if(this->configuration_file.compare(referenced_conffile) != 0) {
						std::ostringstream sout;
						sout << "file " << filename << " references property container definition "
							 << referenced_conffile << " while property container is initialized with container definition file: "
							 << this->configuration_file;
						Logger::get_instance()->emit(LOG_INFO, sout.str());
					}
				} else {
					this->init_from_file(referenced_conffile.c_str());
				}
			} else if(this->initialized == false) {
				std::ostringstream sout;
				sout << "  unable to read data file " << filename << " due to missing definition";
				TDKP_GENERAL_EXCEPTION(sout.str());
			}					
		}
		
		// ignore comments
		if(str_in.size() > 0 && str_in[0] != '#') {
			// strip whitespaces
			this->strip_whitespaces(str_stripped = str_in);
			if(str_stripped.size() > 0) {
				loc = str_stripped.find("=", 0);
				// check if line is valid
				if(loc == string::npos) {
					throw new Exception(1000 + strlen(filename), "in file %s on line %d: invalid line detected. use <propertykey> = <value>", filename, line);
				}
				key   = str_stripped.substr(0, loc);

				try {
					value = lexical_cast<T>(str_stripped.substr(loc + 1));
				} catch(boost::bad_lexical_cast &) {
					ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": could not parse string \"" << str_stripped.substr(loc + 1) << "\"";
					TDKP_GENERAL_EXCEPTION(sout.str());
				}												
				// check if key is valid
				if(!this->valid_key(key)) {
					throw new Exception(1000 + strlen(filename) + key.size() , "in file %s on line %d: %s is not a valid key", filename, line, 	key.c_str());
				}
				// check if value is reasonable
				if(!this->valid_properties[key].inside_error_bounds(value)) {
					std::ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": value " << value 
						 << " for key " << key << " is outside error bounds ]" 
						 << this->valid_properties[key].error_min << ","
						 << this->valid_properties[key].error_max << "[";
					TDKP_GENERAL_EXCEPTION(sout.str());
				}				
				// emit warning if value is out of warning bounds but accept it
				if(!this->valid_properties[key].inside_warn_bounds(value)) {
					std::ostringstream sout;
					sout << "in file " << filename << " on line " << line << ": value " << value 
						 << " for key " << key << " is unreasonable ]" 
						 << this->valid_properties[key].warn_min << ","
						 << this->valid_properties[key].warn_max << "[";
					Logger::get_instance()->emit(LOG_WARN,sout.str());
				}				

				// check if key was already set				
				if(find(keys_in_file.begin(), keys_in_file.end(), key) != keys_in_file.end()) {
					std::ostringstream sout;
					sout << "in file " << filename << " on line " << line 
					     << ": duplicate key " << key << " found";					     
					TDKP_GENERAL_EXCEPTION(sout.str());	
				}
				// finally, set value	
				this->value_map[key] = value;
				keys_in_file.push_back(key);
			} // end if(str_stripped ..)
		} // end if(str_in.size() ...)
		line++;
	} // end while
	fin.close();
}

/** erases all defined values in value map and resets the all values to the available default values */
template<class T>
void PropertyContainer<T>::clear() {
	this->value_map.clear();
	this->set_default_values();	
}

/** set all default values */
template<class T>
void PropertyContainer<T>::set_default_values() {	
	typename map<string, PropDef>::iterator it;
	for(it = this->valid_properties.begin(); it != this->valid_properties.end(); it++) {
		if((*it).second.default_value_avl) {
			this->value_map[(*it).first] = (*it).second.default_value;	
		}	
	}	
}

/** returns true if all mandatory fields have a value 
 */
template<class T>
bool PropertyContainer<T>::valid() {	
	typename map<string, PropDef>::iterator it;
	this->check_proper_initialization();	
	// check all properties and check whether a mandatory field is set
	for(it = this->valid_properties.begin(); it != this->valid_properties.end(); it++) {
		if((*it).second.mandatory) {
			if(!this->is_set((*it).first)) {
				ostringstream sout;
				sout << "key " << (*it).first << " is mandatory, but not set";
				Logger::get_instance()->emit(LOG_WARN, sout.str());
				return false;	
			}
		}	
	}			
	return true;
}

/** return true if value is in valid range
 * 
 * throws exception if invalid key is passed
 * 
 * @param key key as defined in PropertyContainer config file
 * @param value value to test
 * @return true if value is inside ]error_min,error_max[, else false
 *  */
template<class T>
bool PropertyContainer<T>::valid_value(const char* key, T value) throw(Exception *)  {
	string tmp(key);
	return this->valid_value(tmp, value);		
}
template<class T>
bool PropertyContainer<T>::valid_value(const string &key, T value) throw(Exception *) {
	this->check_proper_initialization();
	if(!this->valid_key(key)) {
		TDKP_GENERAL_EXCEPTION("invalid key passed to property container");	
	}	
	if(this->valid_properties[key].error_min < value && this->valid_properties[key].error_max > value) {
		return true;	
	} else {
		return false;	
	}
}

/** check whether key exists in PropertyContainer
 * 
 * does not imply that the value of key is really set. 
 * @param key any property key for any value inside the container
 * @return true, if the property key is defined in the configuration file
 */
template<class T>
bool PropertyContainer<T>::valid_key(const char* key) const {
	string tmp(key);
	return this->valid_key(tmp);	
}

template<class T>
bool PropertyContainer<T>::valid_key(const string& key) const {
	this->check_proper_initialization();	
	if(this->valid_properties.find(key) == this->valid_properties.end()) {
		return false;	
	} else {
		return true;	
	}
}

/** check if any value for property key is set
 * 
 * if param is not a valid key, we throw an exception
 * 
 * @param key valid key as defined in the containers configuration file 
 * @return true if key has either a default value (only possible for optional values) or is set
 */
template<class T>
bool PropertyContainer<T>::is_set(const char* key) const throw(Exception*) {	
	return this->is_set(string(key));
}
template<class T>
bool PropertyContainer<T>::is_set(const string& key) const throw(Exception*) {
	this->check_proper_initialization();	
	if(!this->valid_key(key)) {
		throw new Exception(1000 + key.size(), "invalid key %s requested", key.c_str());
	}
	if(this->value_map.find(key) == this->value_map.end()) {
		return false;	
	} else {
		return true;	
	}	
}	

/** return value of valid key
 * 
 * @param valid key as defined in container configuration
 * @return value of key if set. if not set, an exception is thrown
 */
template<class T>
T PropertyContainer<T>::get(const char* key) const throw(Exception*) {
	string tmp(key);
	return this->get(tmp);			
}
template<class T>
T PropertyContainer<T>::get(const string& key) const throw(Exception*) {			
	this->check_proper_initialization();
	typename map<string,T>::const_iterator it = this->value_map.find(key);
	
	TDKP_ASSERT(it != this->value_map.end(), "value for " << key.c_str() << " requested which is not aviailable");	
	
	return (*it).second;
}

/** set value for key
 * 
 * if the passed value is outside of the defined error bounds, an exception will be thrown ...
 * 
 * @param key valid key defined in container configuration. exception thrown if key is invalid
 */
template<class T>
void PropertyContainer<T>::set(const char* key, T value) throw(Exception*) {
	string tmp(key);
	this->set(tmp, value);	
}

template<class T>
void PropertyContainer<T>::set(const string& key, T value) throw(Exception*) {
	this->check_proper_initialization();	
	if(!this->valid_key(key)) {
		throw new Exception(1000 + key.size(), "invalid key %s", key.c_str());	
	}
	if(!this->valid_properties[key].inside_error_bounds(value)) {
		throw new Exception(1000 + key.size(), "value for key %s is outside error bounds", key.c_str());			
	}	
	if(!this->valid_properties[key].inside_warn_bounds(value)) {			
		Logger::get_instance()->emit(LOG_WARN, 1000 + key.size(), "value for key %s is out of usual bounds", key.c_str()); 			
	}
	this->value_map[key] = value;
}

/** set identifier
 * 
 * used internally by MaterialDatabase (to set the key) and during the assembly procedure 
 * by ProblemDefinition<T>::calculate_element_matrices() to determine cached interaction 
 * matrices. (so interaction matrices are build only once for each material)
 */ 
template<class T>
void PropertyContainer<T>::set_id(unsigned int id_) { this->id = id_; }

/** write current data to file */
template<class T>
void PropertyContainer<T>::write_data_to_file(const char* filename) const {
	
	ofstream fout(filename);
	if(fout) { 
		fout << "# tdkpconfdef: " << configuration_file << "\n"
		     << "# - dumped property data\n\n";
		typename map<string, T>::const_iterator it = value_map.begin();
		// count length of lhs
		int max_length = 0;
		while(it != value_map.end()) {
			max_length = max(max_length, static_cast<int>((*it).first.size()));
			it++;							
		}	
		// write data
		it = value_map.begin();
		
		while(it != value_map.end()) {
			fout.setf(ios::left);
			fout << setw(max_length + 1) << (*it).first << " = ";
			fout.setf(ios::right);			
			fout << (*it).second << "\n";
			it++;
		}
		fout.close(); 		     
	} else {
		TDKP_GENERAL_EXCEPTION("can not write to file " << filename);	
	}

}

} // end namespace

#endif //PROPERTYCONTAINER_H_
