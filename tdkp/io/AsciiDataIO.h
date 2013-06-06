/*
 * AsciiDataIO.h
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */

#ifndef ASCIIDATAIO_H_
#define ASCIIDATAIO_H_

#include "BaseDataIO.h"

namespace tdkp {

/** writes (and possibly reads) node and element data to ascii files */
class AsciiDataIO : public BaseDataIO {
public:
	/** constructor for plotting of geometry based data */
	AsciiDataIO(const Geometry* geometry_);
	/** empty constructor for 0D plotting */
	AsciiDataIO();
	virtual ~AsciiDataIO();

	/** write element data to file */
	template<class T>
	void write_element(const ElementData<T>& data, const char* filename) const;
	/** write node data to file */
	template<class T>
	void write_nodal(const NodeData<T>& data, const char* filename) const;

private:

	void write_ascii_value(ofstream& stream, const double& value) const {
		stream << value;
	}

	void write_ascii_value(ofstream& stream, const complex<double> & value) const {
		stream << value.real() << " \t" << value.imag();
	}

	const Geometry* geometry;
};


/** write ascii output of element data
 *
 * writes (coord) (data) (if data is complex -> data -> data.real() \tdata.imag())
 * (coords) are the average of the node coordinates
 */
template<class T>
void AsciiDataIO::write_element(const ElementData<T>& data, const char* filename) const {
	double coordinates[3] = {0.0, 0.0, 0.0};
	ofstream fout(get_elementdata_filename(filename).c_str(), ios::out);
	const int width = 18;
	if(fout) {
		for(unsigned int ee = 0; ee < geometry->get_num_elements(); ee++) {
			coordinates[0] = coordinates[1] = coordinates[2] = 0.0e0;
			const Element& elem = geometry->get_element(ee);
			if(elem.enabled()) {
				for(unsigned int vv = 0; vv < elem.get_num_nodes(); vv++) {
					const Node& node = elem.get_node(vv);
					for(unsigned int cc = 0; cc < geometry->get_dimension(); cc++) {
						coordinates[cc] += node.get_coord(cc);
					}
				}
				for(unsigned int cc = 0; cc < geometry->get_dimension(); cc++) {
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
		TDKP_GENERAL_EXCEPTION("could not write to file " << get_elementdata_filename(filename));
	}
}

/** write ascii output of node data
 *
 * writes (coord) (data) (if data is complex -> data -> data.real() \tdata.imag())
 */
template<class T>
void AsciiDataIO::write_nodal(const NodeData<T>& data, const char* filename) const {
	ofstream fout(get_nodedata_filename(filename).c_str(), ios::out);
	const int width = 18;
	if(fout) {
		for(unsigned int vv = 0; vv < geometry->get_num_nodes(); vv++) {
			const Node& node = 	geometry->get_node(vv);
			for(unsigned int cc = 0; cc < geometry->get_dimension(); cc++) {
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
		TDKP_GENERAL_EXCEPTION("could not write to file " << get_nodedata_filename(filename));
	}
}




}

#endif /* ASCIIDATAIO_H_ */
