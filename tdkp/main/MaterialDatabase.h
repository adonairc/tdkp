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

#ifndef MATERIALDATABASE_H_
#define MATERIALDATABASE_H_

#include <list>
#include <map>
#include <string>
#include "tdkp/common/all.h"
#include "tdkp/common/Exception.h"
#include "tdkp/main/PropertyContainer.h"


using namespace std;

namespace tdkp
{
	
class Material {
public:	
	Material() {}
	virtual ~Material() {}
	virtual double get(const char* value) const throw(Exception*) = 0;
	virtual void   set(const char* key, double value) throw(Exception*) = 0;    
	virtual bool   valid_key(const char* key) const = 0;
	virtual bool   is_set(const char* key) const = 0;	
	virtual bool   valid() const = 0;
	virtual unsigned int  get_id() const = 0;
	virtual void   set_id(unsigned int id_) = 0;
	virtual void   write_to_file(const char* filename) const = 0;		
	
private:
		
};	

/** decorator for property container
 * 
 * initially, a material class was simply a typedef for the
 * property container. but as sebi and i added the tdkp interface
 * to aqua it became necessary to link my material database to his ...
 * 
 * therefore i added the abstract material class and wrapped 
 * a decorator around the property container.
 *  
 */ 
class MaterialPropertyContainer : public Material {

public:	
	MaterialPropertyContainer(PropertyContainer<double>* container);
	virtual ~MaterialPropertyContainer();
	virtual double get(const char* value) const throw(Exception*);
	virtual void   set(const char* key, double value) throw(Exception*);    
	virtual bool   valid_key(const char* key) const;
	virtual bool   is_set(const char* key) const;	
	virtual bool   valid() const;
	virtual unsigned int  get_id() const;
	virtual void   set_id(unsigned int id_);
	virtual void   write_to_file(const char* filename) const;
	
	void dump() const;
	
	static MaterialPropertyContainer* create_zincblende_material();
	static MaterialPropertyContainer* create_wurtzite_material();
	static MaterialPropertyContainer* create_gas_material();
		
private:	
	PropertyContainer<double>* container;
	
};
	
	
	
	
/** material properties container
 * 
 * the material properties container holds the information on all materials avaiable.
 * materials in the calculation are defined by their name (identification string). if 
 * a certain property of a material referenced in the grid file is requested, the 
 * material database tries to find it internally. if it is not yet loaded, it
 * performs a lookup in all path directories to locate a file called
 * <identificationstring>.mat which should hold the material properties. 
 */
class MaterialDatabase
{
public:

	typedef map<string,	Material*>::iterator material_iterator;
	typedef map<string,	Material*>::const_iterator material_const_iterator;

	MaterialDatabase();
	virtual ~MaterialDatabase();
		
	void   add_search_path(const char* path);
	
	double get(const char* mat,   const char* key)   throw(Exception*);	
	void   set(const char*   mat, const char*   key, const double& value) throw(Exception*);
		
	bool   material_exists(const char* key);	
	bool   material_is_valid(const char* key);
			
	bool   property_is_set(const char* mat, const char* key)     throw(Exception*);		
	virtual void load_material(const char* name)   throw(Exception* );		
	
	void   add_material(const string& name, Material* material) throw(Exception*);
	void   add_material(const char*   name, Material* material) throw(Exception*);
	
	const Material* get_material(const char*   key) throw(Exception* );	
	const Material* get_material(const int &id)  throw(Exception* );	
	int   get_num_materials() const { return this->materials.size(); }

	const string& get_material_name(int id) const throw(Exception*);	
							
protected:
	list<string> paths;
	map<string,	Material*>	materials;
	
};



}

#endif /*MATERIALDATABASE_H_*/
