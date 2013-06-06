/*
 * MEDDataIO.h
 *
 *  Created on: Jul 25, 2010
 *      Author: vepi
 */

#ifndef MEDDATAIO_H_
#define MEDDATAIO_H_

#include "tdkp/io/BaseDataIO.h"
#include "tdkp/geometry/Geometry.h"

namespace tdkp {

class MEDDataIOPimpl;

/** writes data into .med files
 *
 * note, when we write to .med files, we simply add all the
 * data into the same file which was selected using the constructor
 */
class MEDDataIO : public BaseDataIO {
public:
	/** write med data based on a geometry object */
	MEDDataIO(const Geometry& geometry, const char* outputfile);
	/** finalizes med file (calls write) and deletes med object */
	virtual ~MEDDataIO();
	/** append element data to file */
	template<class T>
	void write_element(const ElementData<T>& data, const char* varname);
	/** append node data to file */
	template<class T>
	void write_nodal(const NodeData<T>& data, const char* varname);

private:

	void write_element_cplx(const ElementData<cplx>& data, const char* varname);
	void write_nodal_cplx(const NodeData<cplx>& data, const char* varname);
	void write_element_dbl(const ElementData<double>& data, const char* varname);
	void write_nodal_dbl(const NodeData<double>& data, const char* varname);

	MEDDataIOPimpl* impl;

};


/** append element data to file */
template<>
void MEDDataIO::write_element(const ElementData<complex<double> >& data, const char* varname) {
	write_element_cplx(data, varname);
}
/** append node data to file */
template<>
void MEDDataIO::write_nodal(const NodeData<complex<double> >& data, const char* varname) {
	write_nodal_cplx(data, varname);
}

/** append element data to file */
template<>
void MEDDataIO::write_element(const ElementData<double>& data, const char* varname) {
	write_element_dbl(data, varname);
}
/** append node data to file */
template<>
void MEDDataIO::write_nodal(const NodeData<double>& data, const char* varname) {
	write_nodal_dbl(data, varname);
}



}

#endif /* MEDDATAIO_H_ */
