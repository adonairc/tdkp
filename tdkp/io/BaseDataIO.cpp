/*
 * BaseDataIO.cpp
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */


#include "tdkp/io/BaseDataIO.h"

namespace tdkp {

BaseDataIO::BaseDataIO() :
	xy_data_ending(".dat"), element_data_ending(".dat"),
	node_data_ending(".dat")
{

}
BaseDataIO::~BaseDataIO() {

}

/** returns adjusted xydata file name */
string BaseDataIO::get_xydata_filename(const char* filename) const {
	return string(filename) + xy_data_ending;
}
/** returns adjusted node data file name */
string BaseDataIO::get_nodedata_filename(const char* filename) const {
	return string(filename) + node_data_ending;
}
/** returns adjusted element data file name */
string BaseDataIO::get_elementdata_filename(const char* filename) const {
	return string(filename) + element_data_ending;
}

}
