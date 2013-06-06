/*
 * BaseDataIO.h
 *
 *  Created on: Jul 23, 2010
 *      Author: vepi
 */

#ifndef BASEDATAIO_H_
#define BASEDATAIO_H_

#include "tdkp/common/all.h"
#include "tdkp/common/DataTypes.h"


namespace tdkp {

/** basic almost abstract output writer class */
class BaseDataIO {
public:
	BaseDataIO();
	virtual ~BaseDataIO();

	/** write x/y data to file */
	template<class T>
	void write_xy(const XYData<T>& data, const char* filename, int precision = 12) const throw(Exception*);
	/** write element data to file */
	template<class T>
	void write_element(const ElementData<T>& data, const char* filename) const { TDKP_GENERAL_EXCEPTION("not implemented yet"); }
	/** write node data to file */
	template<class T>
	void write_nodal(const NodeData<T>& data, const char* filename) const { TDKP_GENERAL_EXCEPTION("not implemented yet"); }
	/** read element data from file */
	template<class T>
	void read_element(ElementData<T>& data, const char* filename) const { TDKP_GENERAL_EXCEPTION("not implemented yet"); }
	/** read node data from file */
	template<class T>
	void read_nodal(NodeData<T>& data, const char* filename) const { TDKP_GENERAL_EXCEPTION("not implemented yet"); }

	/** returns adjusted xydata file name */
	string get_xydata_filename(const char* filename) const;
	/** returns adjusted node data file name */
	string get_nodedata_filename(const char* filename) const;
	/** returns adjusted element data file name */
	string get_elementdata_filename(const char* filename) const;

protected:

	string xy_data_ending;			//!< string which is appended to xydata file names
	string element_data_ending;		//!< string which is appended to node data file names
	string node_data_ending;		//!< string which is appended to element data file names


private:

};

/** write x/y data to file */
template<class T>
void BaseDataIO::write_xy(const XYData<T>& data, const char* cfilename, int precision) const throw(Exception*) {
	string filename = get_xydata_filename(cfilename);
	ofstream fout(filename.c_str());
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
		TDKP_GENERAL_EXCEPTION("could not write to file " << filename);
	}
}


}

#endif /* OUTPUTWRITER_H_ */
