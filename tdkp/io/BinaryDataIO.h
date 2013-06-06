/*
 * BinaryDataIO.h
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */

#ifndef BINARYDATAIO_H_
#define BINARYDATAIO_H_

#include "BaseDataIO.h"

namespace tdkp {

/** binary data writer and reader */
class BinaryDataIO: public BaseDataIO {
public:
	BinaryDataIO();
	virtual ~BinaryDataIO();

	/** write element data to file */
	template<class T>
	void write_element(const ElementData<T>& data, const char* filename) const;
	/** write node data to file */
	template<class T>
	void write_nodal(const NodeData<T>& data, const char* filename) const;
	/** read element data from file */
	template<class T>
	void read_element(ElementData<T>& data, const char* filename) const;
	/** read node data from file */
	template<class T>
	void read_nodal(NodeData<T>& data, const char* filename) const;

	// ----------------------------------------------
	// binary specific functions
	// ----------------------------------------------
	/** read node data from stream */
	template<class T>
	void read_nodal(NodeData<T>&  data, istream& fin) const throw(Exception*);
	/** write node data to stream */
	template<class T>
	void write_nodal(const NodeData<T>& data, ostream& fout) const throw(Exception*);
	/** read element data from stream */
	template<class T>
	void read_element(ElementData<T>&  data, istream& fin) const throw(Exception*);
	/** write element data to stream */
	template<class T>
	void write_element(const ElementData<T>& data, ostream& fout) const throw(Exception*);


private:
	/** binary data file header and footer definition, attention functions reading
	 *  have only a buffer of length 50 to read headers/footer. adjust if you change definitions here */
	static const char* NodeDataHeader;
	static const char* NodeDataFooter;
	static const char* ElementDataHeader;
	static const char* ElementDataFooter;

};


/** read node data from stream */
template<class T>
void BinaryDataIO::read_nodal(NodeData<T>&  data, istream& fin) const throw(Exception*) {
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "BinaryDataIO: reading data from input stream");
	char header[50];
	fin.read(header, strlen(NodeDataHeader));
	header[strlen(NodeDataHeader)] = '\0';
	TDKP_ASSERT(string(header).compare(NodeDataHeader) == 0, "invalid header");
	int typesize;
	int size;
	int numeq;
	fin.read((char*)(&typesize), sizeof(int));
	fin.read((char*)(&size),     sizeof(int));
	fin.read((char*)(&numeq),    sizeof(int));
	TDKP_ASSERT(typesize == sizeof(T), "datatype size in file is not equal to datatype-size in object");
	TDKP_ASSERT(size > 0 && numeq > 0, "invalid size information");
	data.set_length(size, numeq);
	T* tmp = new T[numeq];
	for(int ii = 0; ii < size && !fin.eof(); ii++) {
		fin.read((char*)tmp, sizeof(T) * numeq);
		for(int jj = 0; jj < numeq; jj++) {
			data.set_node_value(ii,jj, tmp[jj]);
		}
	}
	delete[] tmp;
	fin.read(header, strlen(NodeDataFooter));
	header[strlen(NodeDataFooter)] = '\0';
	TDKP_ASSERT(string(header).compare(NodeDataFooter) == 0, "datafile did not end with correct footer. so there is something wrong.");
}
/** write node data to stream */
template<class T>
void BinaryDataIO::write_nodal(const NodeData<T>& data, ostream& fout) const throw(Exception*) {

	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "BinaryDataIO: writing data to output stream");

	fout.write(NodeDataHeader, strlen(NodeDataHeader));

	int size     	 		= data.get_length();
	int type_size 			= sizeof(T);
	int num_data_per_node = data.get_num_data_per_node();
	fout.write((char*)(&type_size), sizeof(int));
	fout.write((char*)(&size),  sizeof(int));
	fout.write((char*)(&num_data_per_node), sizeof(int));
	T* tmp = new T[num_data_per_node];
	for(int ii = 0; ii < size && !fout.eof(); ii++) {
		for(int jj = 0; jj < num_data_per_node; jj++) {
			tmp[jj] = data.get_node_value(ii,jj);
		}
		fout.write((char*)tmp, sizeof(T) * num_data_per_node);
	}
	fout.write(NodeDataFooter, strlen(NodeDataFooter));
	delete[] tmp;
}
/** read element data from stream */
template<class T>
void BinaryDataIO::read_element(ElementData<T>&  data, istream& fin) const throw(Exception*) {
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "BinaryDataIO: reading data from input stream");

	char header[50];
	fin.read(header, strlen(ElementDataHeader));
	header[strlen(ElementDataHeader)] = '\0';
	TDKP_ASSERT(string(header).compare(ElementDataHeader) == 0, "invalid header");
	int typesize;
	int size;
	int numeq;
	fin.read((char*)(&typesize), sizeof(int));
	fin.read((char*)(&size),  sizeof(int));
	fin.read((char*)(&numeq), sizeof(int));
	TDKP_ASSERT(typesize == sizeof(T), "datatype size in file is not equal to datatype-size in object");
	TDKP_ASSERT(size > 0 && numeq > 0, "invalid size information");
	data.set_length(size, numeq);
	T* tmp = new T[numeq];
	for(int ii = 0; ii < size && !fin.eof(); ii++) {
		fin.read((char*)tmp, sizeof(T) * numeq);
		for(int jj = 0; jj < numeq; jj++) {
			data.set_element_value(ii,jj, tmp[jj]);
		}
	}
	delete[] tmp;
	fin.read(header, strlen(ElementDataFooter));
	header[strlen(ElementDataFooter)] = '\0';
	TDKP_ASSERT(string(header).compare(ElementDataFooter) == 0, "datafile did not end with correct footer. so there is something wrong.");

}
/** write element data to stream */
template<class T>
void BinaryDataIO::write_element(const ElementData<T>& data, ostream& fout) const throw(Exception*) {

	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "BinaryDataIO: writing data to output stream");

	fout.write(ElementDataHeader, strlen(ElementDataHeader));

	int size = data.get_length();
	int num_data_per_element = data.get_num_data_per_element();
	int typesize = sizeof(T);
	fout.write((char*)(&typesize), sizeof(int));
	fout.write((char*)(&size),     sizeof(int));
	fout.write((char*)(&num_data_per_element), sizeof(int));

	T* tmp = new T[num_data_per_element];
	for(int ii = 0; ii < size && !fout.eof(); ii++) {
		for(int jj = 0; jj < num_data_per_element; jj++) {
			tmp[jj] = data.get_element_value(ii,jj);
		}
		fout.write((char*)tmp, sizeof(T) * num_data_per_element);
	}
	fout.write(ElementDataFooter, strlen(ElementDataHeader));
	delete[] tmp;
}


/** write element data to file */
template<class T>
void BinaryDataIO::write_element(const ElementData<T>& data, const char* cfilename) const {
	string filename = get_elementdata_filename(cfilename);
	ofstream fout(filename.c_str(), ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "BinaryDataIO: writing element data to file " << filename);
	if(fout) {
		this->write_element(data, fout);
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);
	}
}
/** write node data to file */
template<class T>
void BinaryDataIO::write_nodal(const NodeData<T>& data, const char* cfilename) const {
	string filename = get_nodedata_filename(cfilename);
	ofstream fout(filename.c_str(), ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "BinaryDataIO: writing node data to file " << filename);
	if(fout) {
		this->write_nodal(data, fout);
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);
	}
}
/** read element data from file */
template<class T>
void BinaryDataIO::read_element(ElementData<T>& data, const char* cfilename) const {
	string filename = get_elementdata_filename(cfilename);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "BinaryDataIO: reading element data from file " << filename);
	ifstream fin(filename.c_str(), ios::binary);
	if(fin) {
		this->read_element(data, fin);
		fin.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not read from file " << filename);
	}
}
/** read node data from file */
template<class T>
void BinaryDataIO::read_nodal(NodeData<T>& data, const char* cfilename) const {
	string filename = get_nodedata_filename(cfilename);
	ifstream fin(filename.c_str(), ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "BinaryDataIO: reading node data from file " << filename);
	if(fin) {
		this->read_nodal(data, fin);
		fin.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not read from file " << filename);
	}
}


}

#endif /* BINARYDATAIO_H_ */
