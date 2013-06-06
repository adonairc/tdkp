/*
 * AsciiDataIO.cpp
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */

#include "tdkp/io/AsciiDataIO.h"

namespace tdkp {

/** constructor for plotting of geometry based data */
AsciiDataIO::AsciiDataIO(const Geometry* geometry_) :
	geometry(geometry_)
{
	this->node_data_ending    = "_ascii.dat";
	this->element_data_ending = "_ascii.dat";
}



/** empty constructor for 0D plotting */
AsciiDataIO::AsciiDataIO() :
	geometry(0)
{

}

AsciiDataIO::~AsciiDataIO() {

}


}
