/*
 * BinaryDataIO.cpp
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */

#include "BinaryDataIO.h"

namespace tdkp {

const char* BinaryDataIO::NodeDataHeader  = "NodeDataBegin";
const char* BinaryDataIO::NodeDataFooter  = "NodeDataEnd";
const char* BinaryDataIO::ElementDataHeader = "ElementDataBegin";
const char* BinaryDataIO::ElementDataFooter = "ElementDataEnd";

BinaryDataIO::BinaryDataIO()
{
	this->node_data_ending    = ".bin";
	this->element_data_ending = ".bin";
}

BinaryDataIO::~BinaryDataIO()
{
}



}
