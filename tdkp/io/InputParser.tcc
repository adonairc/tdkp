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



namespace tdkp {

template<class T> 
void InputParser::read_binary(NodeData<T>&  data, const char* filename) const throw(Exception*) {
	ifstream fin(filename, ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "InputParser: reading node data from file " << filename);
	if(fin) {
		this->read_binary(data, fin);
		fin.close();	
	} else {
		TDKP_GENERAL_EXCEPTION("could not read from file " << filename);	
	}
}

template<class T> 		
void InputParser::write_binary(const NodeData<T>& data, const char* filename) const throw(Exception*) {
	ofstream fout(filename, ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "InputParser: writing node data to file " << filename); 
	if(fout) {
		this->write_binary(data, fout); 
		fout.close(); 	
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);	
	}
}

template<class T> 
void InputParser::read_binary(ElementData<T>&  data, const char* filename) const throw(Exception*) {
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "InputParser: reading element data from file " << filename);
	ifstream fin(filename, ios::binary);
	if(fin) {
		this->read_binary(data, fin);
		fin.close();	
	} else {
		TDKP_GENERAL_EXCEPTION("could not read from file " << filename);	
	}	
}

template<class T> 		
void InputParser::write_binary(const ElementData<T>& data, const char* filename) const throw(Exception*) {
	ofstream fout(filename, ios::binary);
	TDKP_LOGMSG(LOG_INFO_DEVEL1, "InputParser: writing element data to file " << filename); 
	if(fout) {
		this->write_binary(data, fout); 
		fout.close(); 	
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);	
	}	
}
	
	
template<class T>
void InputParser::read_binary(NodeData<T>& data, istream& fin) const throw(Exception*) {
	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "InputParser: reading data from input stream");	
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

template<class T> 
void InputParser::write_binary(const NodeData<T>& data, ostream& fout) const throw(Exception*) {
	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "InputParser: writing data to output stream");	
	
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

template<class T>
void InputParser::read_binary(ElementData<T>& data, istream& fin)  const throw(Exception*) {
	
	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "InputParser: reading data from input stream");	

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
template<class T> 
void InputParser::write_binary(const ElementData<T>& data, ostream& fout)  const throw(Exception*) {

	Logger::get_instance()->emit(LOG_INFO_DEVEL2, "InputParser: writing data to output stream");	
	
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
	
template<class T> 
void InputParser::write_ascii_real(const XYData<T>& data, const char* filename, int precision) const throw(Exception*) {
	XYDataReal<T> real(data); 
	this->write_ascii(real, filename, precision); 	
}

template<class T> 
void InputParser::write_ascii(const XYData<T>& data, const char* filename, int precision) const throw(Exception*) {
		
	ofstream fout(filename); 
	if(fout) {
		// get data
		const int x_length = data.get_x_length();
		const int y_sets   = data.get_num_y_sets();
		const int x_sets   = data.get_num_x_sets();
		
		
		vector<vector<T> > tmpdat(y_sets + x_sets);
		vector<string>     ident(y_sets + x_sets);
		int width = 0;
		for(int ii = 0; ii < x_sets; ii++) {
			data.get_x(ii, tmpdat[ii]);		
			ident[ii] = data.get_x_identifier(ii);
			width = max(width, static_cast<int>(ident[ii].size()));
		}		
		
		
		// get data and strings
		for(int ii = 0; ii < y_sets; ii++) {
			data.get_y(ii,tmpdat[ii+x_sets]);
			ident[ii + x_sets] = data.get_y_identifier(ii);			
			width = max(width, int(ident[ii + x_sets].size()));
		}
		if(width < precision + 3) {
			width = precision + 3;	
		}
		// write head
		fout.precision(precision); 
		fout << "#";
		for(int ii = 0; ii < y_sets + x_sets; ii++) {
			fout << setw(width) << ident[ii] << " ";		
		}
		fout << "\n";
		fout.setf(ios::scientific | ios::right); 		
		for(int ii = 0; ii < x_length; ii++) {
			for(int jj = 0; jj < y_sets + x_sets; jj++) {
				fout << setw(width) << tmpdat[jj][ii] << " ";
			}
			fout << "\n";
		}
		fout.close(); 					
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file");	
	}
}


/** write ascii output of element data
 * 
 * writes (coord) (data) (if data is complex -> data -> data.real() \tdata.imag())
 * (coords) are the average of the node coordinates
 */ 
template<class T>
void InputParser::write_ascii(const Geometry& geometry, const ElementData<T>& data, const char* filename) const {

	double coordinates[3] = {0.0, 0.0, 0.0};
	ofstream fout(filename, ios::out);
	const int width = 18;  	
	if(fout) {
		for(unsigned int ee = 0; ee < geometry.get_num_elements(); ee++) {
			coordinates[0] = coordinates[1] = coordinates[2] = 0.0e0;
			const Element& elem = geometry.get_element(ee);
			if(elem.enabled()) {
				for(unsigned int vv = 0; vv < elem.get_num_nodes(); vv++) {
					const Node& node = elem.get_node(vv);
					for(unsigned int cc = 0; cc < geometry.get_dimension(); cc++) {
						coordinates[cc] += node.get_coord(cc);	
					}		
				}
				for(unsigned int cc = 0; cc < geometry.get_dimension(); cc++) {
					coordinates[cc] /= double(elem.get_num_nodes());
					fout << setw(width) << coordinates[cc] << " \t";
				}
				for(int ii = 0; ii < data.get_num_data_per_element(); ii++) {
					fout << setw(width);
					write_ascii_value(fout, data.get_element_value(ee,ii));
					fout << " \t";	
				} 
				fout << "\n";
			}				
		}
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);	
	} 
		
}

/** write ascii output of node data
 * 
 * writes (coord) (data) (if data is complex -> data -> data.real() \tdata.imag())
 */ 
template<class T>
void InputParser::write_ascii( const Geometry& geometry, const NodeData<T>& data, const char* filename) const {

	ofstream fout(filename, ios::out);
	const int width = 18;  	
	if(fout) {
		for(unsigned int vv = 0; vv < geometry.get_num_nodes(); vv++) {
			const Node& node = 	geometry.get_node(vv);
			for(unsigned int cc = 0; cc < geometry.get_dimension(); cc++) {				
				fout << setw(width) << node.get_coord(cc) << " \t";
			}
			for(int ii = 0; ii < data.get_num_data_per_node(); ii++) {
				write_ascii_value(fout, data.get_node_value(vv,ii));
				fout << " \t";
			} 
			fout << "\n";				
		}
		fout.close();
	} else {
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);	
	}
}

} // end of namespace tdkp
