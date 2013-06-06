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


#include <iostream>
#include <sstream>
#include "tdkp/main/MaterialDatabase.h"
using namespace std;


namespace tdkp {
	
	
	
// ------------------------------------------------------
// implementation of Property Container Wrapper
// simply a decorator 	
// ------------------------------------------------------		
MaterialPropertyContainer::MaterialPropertyContainer(PropertyContainer<double>* container_)
:container(container_) 
{
	TDKP_ASSERT(container != 0, "container != 0");	
}
MaterialPropertyContainer::~MaterialPropertyContainer() {
	delete container;
	container = 0;	
}
double MaterialPropertyContainer::get(const char* value) const throw(Exception*) {
	return container->get(value);	
}
void   MaterialPropertyContainer::set(const char* key, double value) throw(Exception*) {
	container->set(key, value);	
}    
bool   MaterialPropertyContainer::valid_key(const char* key) const {
	return container->valid_key(key);	
}
bool   MaterialPropertyContainer::is_set(const char* key) const {
	return container->is_set(key);	
}	
bool   MaterialPropertyContainer::valid() const {
	return container->valid();	
}
unsigned int MaterialPropertyContainer::get_id() const {
	return container->get_id();	
}
void   MaterialPropertyContainer::set_id(unsigned int id_) {
	return container->set_id(id_);	
}	
void MaterialPropertyContainer::dump() const {	
	TDKP_LOGMSG(LOG_INFO, *container) ;		
}

/** write container data to file */
void MaterialPropertyContainer::write_to_file(const char* filename) const {
	container->write_data_to_file(filename);	
}


/** create a material property container object and loads the data definitions 
 * 
 * if you like to define your materials inside tcl files, you may use this one
 * to create a corresponding material object
 */
MaterialPropertyContainer* MaterialPropertyContainer::create_zincblende_material() {
	return new MaterialPropertyContainer(new PropertyContainer<double>("material_zincblende.cnf"));	
}
/** material property container for wurtzite crystals */
MaterialPropertyContainer* MaterialPropertyContainer::create_wurtzite_material() {
	return new MaterialPropertyContainer(new PropertyContainer<double>("material_wurtzite.cnf"));		
}
/** material property container for gas */
MaterialPropertyContainer* MaterialPropertyContainer::create_gas_material() {
	return new MaterialPropertyContainer(new PropertyContainer<double>("material_gas.cnf"));		
}

// ------------------------------------------------------
// implementation of material database 	
// ------------------------------------------------------	

/** main constructor adds current directory to search path per default */
MaterialDatabase::MaterialDatabase() {
	paths.push_back(".");	
}


MaterialDatabase::~MaterialDatabase() {
	for(material_iterator it = materials.begin(); it != materials.end(); it++) {
		if((*it).second != 0) {
			delete (*it).second;
			(*it).second = 0;	
		}	
	}		
}

/** add path to search path for search material files **/
void MaterialDatabase::add_search_path(const char* cpath_) {
	string path(cpath_);
	if(path.size() == 0) {
		TDKP_GENERAL_EXCEPTION("tried to add empty path to search path");			
	}
	paths.push_back(path);
}
		
/** get value for some materials property
 * 
 * @param mat material identifier. 
 * @param key property key as defined in a property container config file
 * @return the property value. if either material parameter can not be found or key is invalid, an exception is thrown
 */ 
double MaterialDatabase::get(const char* cmat,   const char* ckey) throw(Exception*) {
	string mat(cmat);
	string key(ckey);
	if(this->property_is_set(cmat, ckey))	{
		return this->materials[mat]->get(key.c_str());
	}	
	ostringstream sout; sout << "material property " << key << " for " << mat << " is not set";
	TDKP_GENERAL_EXCEPTION(sout.str());	
}


void MaterialDatabase::set(const char* mat, const char* key, const double& value) throw(Exception*) {
	if(this->material_exists(mat)) {
		if(this->materials[mat]->valid_key(key)) {
			this->materials[mat]->set(key,value); 
		} else {
			ostringstream sout;
			sout << "invalid key " << key << " for material " << mat;
			TDKP_GENERAL_EXCEPTION(sout.str());	
		}
	} else {
		ostringstream sout;
		sout << "material " << mat << " is not loaded";
		TDKP_GENERAL_EXCEPTION(sout.str());		 
	}	
}
/** check whether material properties for key are available 
 * 
 * @param key material key as referenced in grid file
 * */
bool MaterialDatabase::material_exists(const char* ckey) {
	string key(ckey);
	material_iterator it;
	if((it = this->materials.find(key)) != this->materials.end()) {
		return true;	
	} else {
		return false;
	}		
}

/** check whether material has all mandatory values set and all values are inside given bounds 
 * 
 * @param key the material key as defined in grid file
 */
bool MaterialDatabase::material_is_valid(const char* key) {
	if(!this->material_exists(key)) {
		std::ostringstream sout; sout << "unloaded material " << key << " requested";
		TDKP_GENERAL_EXCEPTION(sout.str());
		return false;	
	}
	return this->materials[string(key)]->valid();
}

/** check if a property of a material is set
 * @param mat material key
 * @param key property key
 * @return true if property is set, false if not set, exception is thrown if mat is not available or 
 *         key is not a valid key (as defined in the configuration file)
 */
bool MaterialDatabase::property_is_set(const char* mat, const char* key)     throw(Exception*) {
	if(this->material_exists(mat)) {
		return this->get_material(mat)->is_set(key);		
	} else {
		std::ostringstream sout; sout << "request for " << key << " of material " << mat 
		<< " is invalid as material was never loaded";
		TDKP_GENERAL_EXCEPTION(sout.str());
	}
}

/** return pointer to the property container having the informatin on the material
 * 
 * @return the pointer. if material key is invalid, 0 is returned
 */
const Material* MaterialDatabase::get_material(const char* ckey) throw (Exception*) {
	string key(ckey);
	material_iterator it;
	if((it = this->materials.find(key)) != this->materials.end()) {
		return (*it).second;
	} else {
		TDKP_GENERAL_EXCEPTION("invalid material requested");
	}		
}			

	
const Material* MaterialDatabase::get_material(const int &id) throw(Exception*) {
	material_const_iterator it;
	for(it = this->materials.begin(); it != this->materials.end(); it++) {
		if((signed)(*it).second->get_id() == id) {
			return (*it).second;	
		}	
	}		
	cout << "available materials: \n";
	for(it = this->materials.begin(); it != this->materials.end(); it++) {
		cout << (*it).first << " with id " << (*it).second->get_id() << "\n";
	}		
	abort(); // create core dump
	TDKP_GENERAL_EXCEPTION("invalid material with id " << id << " requested");
}

/** return name of material with id */
const string& MaterialDatabase::get_material_name(int id) const throw(Exception*) {
	material_const_iterator it;
	for(it = this->materials.begin(); it != this->materials.end(); it++) {
		if((signed)(*it).second->get_id() == id) {
			return (*it).first;	
		}	
	}	
	TDKP_GENERAL_EXCEPTION("invalid material requested");	
}

/** tries to load the material from filesystem if it does not exist in the db
 *  
 */				
void MaterialDatabase::load_material(const char* cname) throw (Exception*) {
	
	if(this->material_exists(cname)) {
		ostringstream sout;
		sout << "material " << cname << " already exists in the material database. skipping loading.";
		Logger::get_instance()->emit(LOG_WARN, sout.str());
		return;	
	}
	
	string name(cname);
	list<string>::const_iterator it;
	string file;
	ifstream fin;
	ostringstream sout;
	for(it = paths.begin(); it != paths.end(); it++) {
		file = (*it) + ((*it)[(*it).size() - 1] == '/' ? "":"/") + name + ".mat";
		fin.open(file.c_str());
		if(fin) {
			sout << "reading material parameters for " << name << " from file " << file;
			Logger::get_instance()->emit(LOG_INFO, sout.str());
			fin.close();
			PropertyContainer<double>* prop = new PropertyContainer<double>();
			prop->read_data_from_file(file.c_str());
			prop->valid();
			this->add_material(name, new MaterialPropertyContainer(prop));
			return;				
		}			
		fin.close(); fin.clear(); // needed due to bug in libg++					
	}
	sout << "unable to locate material file " << name << ".mat in search path:\n";
	for(it = paths.begin(); it != paths.end(); it++) {
		sout << (*it) << "\n";	
	}
	TDKP_GENERAL_EXCEPTION(sout.str());	
}


void MaterialDatabase::add_material(const string& name, Material* prop) throw(Exception*) {
	if(this->material_exists(name.c_str())) {
		ostringstream sout; sout << "material " << name << " already defined";
		TDKP_GENERAL_EXCEPTION(sout.str());
		return;	
	}
	prop->set_id(this->materials.size());
	this->materials[name] = prop;
	return;
}

void MaterialDatabase::add_material(const char* name, Material* material) throw(Exception*) {
	add_material(string(name), material);	
}


} // end namepsace
