/*
 * TCLDataIOWrapper.h
 *
 *  Created on: Jul 25, 2010
 *      Author: vepi
 */

#ifndef TCLDATAIOWRAPPER_H_
#define TCLDATAIOWRAPPER_H_

namespace tdkp {

/** templated wrapper to use data io objects within tcl wrapper
 *
 * this wrapper class reduces the template definition hassle in tdkpshell.i
 * and provides the tcl interface with some logical function and function names
 */
template<class T>
class TCLDataIOWrapper {
public:
	TCLDataIOWrapper(T* io_obj);
	virtual ~TCLDataIOWrapper();

	/** write real part of x/y data to file (compatibility to cplx version) */
	void write_real(const XYData<double>& data, const char* filename, int precision = 12) throw(Exception*);
	/** write x/y data to file */
	void write(const XYData<double>& data, const char* filename, int precision = 12) throw(Exception*);
	/** write element data to file */
	void write(const ElementData<double>& data, const char* filename);
	/** write node data to file */
	void write(const NodeData<double>& data, const char* filename);
	/** read element data from file */
	void read(ElementData<double>& data, const char* filename);
	/** read node data from file */
	void read(NodeData<double>& data, const char* filename);

	/** write real part of x/y data to file (compatibility to cplx version) */
	void write_real(const XYData<complex<double> >& data, const char* filename, int precision = 12) throw(Exception*);
	/** write x/y data to file */
	void write(const XYData<complex<double> >& data, const char* filename, int precision = 12) throw(Exception*);
	/** write element data to file */
	void write(const ElementData<complex<double> >& data, const char* filename);
	/** write node data to file */
	void write(const NodeData<complex<double> >& data, const char* filename);
	/** read element data from file */
	void read(ElementData<complex<double> >& data, const char* filename);
	/** read node data from file */
	void read(NodeData<complex<double> >& data, const char* filename);

private:
	T* io_obj;
};

template<class T>
TCLDataIOWrapper<T>::TCLDataIOWrapper(T* io_obj_) :
	io_obj(io_obj_)
{

}
template<class T>
TCLDataIOWrapper<T>::~TCLDataIOWrapper() {
	delete_object(io_obj);
}

/** write x/y data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const XYData<double>& data, const char* filename, int precision) throw(Exception*) {
	io_obj->write_xy(data, filename, precision);
}
/** write element data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const ElementData<double>& data, const char* filename) {
	io_obj->write_element(data, filename);
}
/** write node data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const NodeData<double>& data, const char* filename) {
	io_obj->write_nodal(data, filename);
}
/** read element data from file */
template<class T>
void TCLDataIOWrapper<T>::read(ElementData<double>& data, const char* filename) {
	io_obj->read_element(data, filename);
}
/** read node data from file */
template<class T>
void TCLDataIOWrapper<T>::read(NodeData<double>& data, const char* filename) {
	io_obj->read_nodal(data, filename);
}

/** write x/y data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const XYData<complex<double> >& data, const char* filename, int precision) throw(Exception*) {
	io_obj->write_xy(data, filename, precision);
}
/** write element data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const ElementData<complex<double> >& data, const char* filename) {
	io_obj->write_element(data, filename);
}
/** write node data to file */
template<class T>
void TCLDataIOWrapper<T>::write(const NodeData<complex<double> >& data, const char* filename) {
	io_obj->write_nodal(data, filename);
}
/** read element data from file */
template<class T>
void TCLDataIOWrapper<T>::read(ElementData<complex<double> >& data, const char* filename) {
	io_obj->read_element(data, filename);
}
/** read node data from file */
template<class T>
void TCLDataIOWrapper<T>::read(NodeData<complex<double> >& data, const char* filename) {
	io_obj->read_nodal(data, filename);
}

/** write real part of x/y data to file (compatibility to cplx version) */
template<class T>
void TCLDataIOWrapper<T>::write_real(const XYData<double>& data, const char* filename, int precision) throw(Exception*) {
	io_obj->write_xy(data, filename, precision);
}

/** write real part of x/y data to file (compatibility to cplx version) */
template<class T>
void TCLDataIOWrapper<T>::write_real(const XYData<complex<double> >& data, const char* filename, int precision) throw(Exception*) {
	XYDataReal<cplx> real(data);
	io_obj->write_xy(real, filename, precision);
}

}

#endif /* TCLDATAIOWRAPPER_H_ */
