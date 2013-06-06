#ifndef VALIDATOR_H_
#define VALIDATOR_H_


class ElementNumericalValidator {
public:

	static const bool be_chatty = true;
	string message;

	ElementNumericalValidator(const Element& element);


	void test_evaluate_integral(const char* msg) {
		this->message = msg;
		double volume = evaluate_integral(&ElementNumericalValidator::one);
		if(be_chatty) {
			cout << "-----------------------------------------------\n"
			     << " running " << msg << "\n"; 
			cout << "evaluate_integral returns: " 
			     << volume
			     << " while the element volume is: "
			     << element.get_volume() 
			     << "\n";
		}
		TSM_ASSERT_DELTA(message, volume, element.get_volume(), ERROR_THRESHOLD);				     
		test_integrate_solution();
		test_differentiate_solution();
		test_derivative();		     
		test_zero_order();
		test_zero_order_nodal();		 
		test_first_order();   
		test_second_order();
		test_zero_order_single();
		test_first_order_single();		
	}

	void plot_shape_functions(const char* filename);
	
	void test_integrate_solution();
	void test_differentiate_solution();
	void test_derivative();		     
	void test_zero_order();
	void test_zero_order_nodal();		 
	void test_first_order();   
	void test_second_order();
	void test_zero_order_single();
	void test_first_order_single();	
	
	double one(const double& gx, const double& gy, const double& gz) const { return 1.0; }
	
	double zero_order(const double& gx, const double& gy, const double& gz) const;
	double zero_order_nodal(const double& gx, const double& gy, const double& gz) const;
	double first_order(const double& gx, const double& gy, const double& gz) const;
	double second_order(const double& gx, const double& gy, const double& gz) const;
	double zero_order_single(const double& gx, const double& gy, const double& gz) const;
	double first_order_single(const double& gx, const double& gy, const double& gz) const;	
	
	double get_derivative(short elem_shape_func, short diffop, const double& gx, const double& gy, const double& gz) const; 
	
private:
	double evaluate_integral(double (ElementNumericalValidator::*func)(const double&,const double&,const double&) const) const;
	struct PatchCoordinates {
		PatchCoordinates(double gx_, double gy_, double gz_) { gx = gx_; gy = gy_; gz = gz_; }
		double gx; double gy; double gz;	
	};
	
	const Element& element;
	Node low;
	Node high;
	double patch_weight;	
	
	vector<PatchCoordinates> patches;
	
	double increments[3];
	
	short diffops[2];
	short elem_shape_funcs[2];
	short nodal_data_func;
	
};

ElementNumericalValidator::ElementNumericalValidator(const Element& element_) 
: element(element_),
  low(0.0, 0.0, 0.0),
  high(0.0, 0.0, 0.0),
  patch_weight(1.0),
  nodal_data_func(0)
{	
	const int rough_patch_num = 250; // was 400		
	element.get_bounding_box(low, high);
				 			
	// get increments (rough estimate)
	for(unsigned short ii = 0; ii < element.get_dimension(); ii++) {
		increments[ii] = (high.get_coord(ii) - low.get_coord(ii)) / double(rough_patch_num - 1);
		patch_weight *= increments[ii];	
	}
			
	int xx[] = {0, 0, 0};
	xx[1]    = (element.get_dimension() < 2 ? rough_patch_num : 0); 
	xx[2]    = (element.get_dimension() < 3 ? rough_patch_num : 0);
	
	double tmp_coordinates[] = {0.0, 0.0, 0.0};
	
	while(xx[0] < rough_patch_num) {
		// calculate position
		for(short ii = 0; ii < (signed)element.get_dimension(); ii++) {
			tmp_coordinates[ii] = low.get_coord(ii) + double(xx[ii]) * increments[ii];
		}
		// test if its inside and add to valid patches if so
		//cout << "testing " << tmp_coordinates[0] << ", " << tmp_coordinates[1] << ", " << tmp_coordinates[2] << endl;		
		if(element.inside_element(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2])) {
			patches.push_back(PatchCoordinates(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2]));
			//cout << "patch " << tmp_coordinates[0] << ", " << tmp_coordinates[1] << ", " << tmp_coordinates[2] << endl;					
		} 

		
		// -----------------------------------
		// update position (depending on dimension)
		// -----------------------------------
		xx[2]++;
		if(xx[2] >= rough_patch_num) {
			xx[2] = (element.get_dimension() < 3 ? rough_patch_num : 0);
			xx[1]++;
			if(xx[1] >= rough_patch_num) {
				xx[1] = (element.get_dimension() < 2 ? rough_patch_num : 0);
				xx[0]++;
			} 	
		}
		 	
	}	
}

double ElementNumericalValidator::zero_order(const double& gx, const double& gy, const double& gz) const {
	return element.evaluate_form_function_global(elem_shape_funcs[0], gx, gy, gz)
	     * element.evaluate_form_function_global(elem_shape_funcs[1], gx, gy, gz);	
}

double ElementNumericalValidator::zero_order_nodal(const double& gx, const double& gy, const double& gz) const {
	return element.evaluate_form_function_global(nodal_data_func, gx, gy, gz)
	     * element.evaluate_form_function_global(elem_shape_funcs[0], gx, gy, gz)
	     * element.evaluate_form_function_global(elem_shape_funcs[1], gx, gy, gz);	
}

double ElementNumericalValidator::get_derivative(short elem_shape_func, short diffop, const double& gx, const double& gy, const double& gz) const {
	// finite differences meets finite elements
	double off_coordinates[] = {gx, gy, gz};
	double cache             = off_coordinates[diffop];	
	double distance          = increments[diffop] / 2.5;
	int    decrements        = 0;
	bool   swap              = false;
	do {
		distance /= 2.0;
		off_coordinates[diffop] = cache + distance;
		decrements++;
	} while (!element.inside_element(off_coordinates[0], off_coordinates[1], off_coordinates[2]) && decrements < 10);
	if(decrements == 10) {
		decrements = 0;
		distance   = increments[diffop] / 2.5;
		do {
			distance /= 2.0;
			off_coordinates[diffop] = cache - distance;
			decrements++;	
		} while (!element.inside_element(off_coordinates[0], off_coordinates[1], off_coordinates[2]) && decrements < 10);
		swap = true;				
	}
	double left  = element.evaluate_form_function_global(elem_shape_func, off_coordinates[0], off_coordinates[1],off_coordinates[2]);
	double right = element.evaluate_form_function_global(elem_shape_func, gx, gy, gz);
	double deriv = (left - right) / (distance);
	if(swap) {
		deriv *= -1.0;	
	}
	/*
	cout << " --------------------------------------------\n"
	     << "(f(" << off_coordinates[0] << ", " << off_coordinates[1] << ") =  "
	     << left << "\n"
	     << "f(" << gx << ", " << gy << ") = " << right << "\n"
	     << "(" << left << " - " << right << ") / " << (distance) << " = "
	     <<  deriv << "\n";   	              
	}*/
	return deriv;
	     	       	
}

double ElementNumericalValidator::first_order(const double& gx, const double& gy, const double& gz) const {	
	double deriv = this->get_derivative(elem_shape_funcs[0], diffops[0], gx, gy, gz);
	return deriv * element.evaluate_form_function_global(elem_shape_funcs[1], gx, gy, gz);	             			
}

double ElementNumericalValidator::second_order(const double& gx, const double& gy, const double& gz) const {
	return this->get_derivative(elem_shape_funcs[0], diffops[0], gx, gy, gz)
	     * this->get_derivative(elem_shape_funcs[1], diffops[1], gx, gy, gz);	     	
}

double ElementNumericalValidator::zero_order_single(const double& gx, const double& gy, const double& gz) const {
	return element.evaluate_form_function_global(elem_shape_funcs[0], gx, gy, gz); 	
}

double ElementNumericalValidator::first_order_single(const double& gx, const double& gy, const double& gz) const {	
	return this->get_derivative(elem_shape_funcs[0], diffops[0], gx, gy, gz);
}


/** numerically evaluate some integral on the element
 * 
 * @param func callback function that  
 */
double ElementNumericalValidator::evaluate_integral(double (ElementNumericalValidator::*func)(const double&,const double&,const double&) const) const {
	double res = 0.0;	
	int patch_num = patches.size();
	
#pragma omp parallel for default(shared) reduction(+:res)
	for(int ii = 0; ii < patch_num; ii++) {
		const PatchCoordinates& pt = this->patches[ii]; 		
		res += (this->*func)(pt.gx, pt.gy, pt.gz);		
	}		
//	cout << "final = " << res * patch_weight << "\n";
	return res * patch_weight;
}

/** test derivatives of shape functions */
void ElementNumericalValidator::test_derivative() {
	// num tests
	unsigned int num_tests = 6;

	// init rnd	
    timeval tp;
    gettimeofday(&tp,NULL);
    srand48(tp.tv_usec);        	
	
	vector<double> local_coords(3, 0.0);
	vector<double> rand_coords(3, 0.0);
	vector<double> rand_global_coords(3, 0.0);
		
	// test derivatives at random positions inside element  
	for(unsigned int nn = 0; nn < num_tests; nn++) {		
		// -----------------------------------------------
		// create random position inside element (note: element is assumed
		// to have a convex hull ...
		// -----------------------------------------------
		rand_coords.assign(3, 0.0);
		rand_global_coords.assign(3, 0.0);
		double total_rand = 0.0;
		double my_rand    = 0.0;
		for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
			// if its a vertex			
			const Node& nd = element.get_node(ii);
			if(nd.get_index_vertex() != -1) {
				// get local node coords
				element.get_node_local_coords(ii, local_coords);
				my_rand = drand48();
				total_rand += my_rand;
				for(unsigned cc = 0; cc < element.get_dimension(); cc++) {
					rand_coords[cc] += my_rand * local_coords[cc];
					rand_global_coords[cc] += my_rand * nd.get_coord(cc);								
				}
			}	
		}
		TDKP_ASSERT(total_rand > 0.0, "no vertex defined!");
		// normalize interpolation to 1
		for(unsigned cc = 0; cc < element.get_dimension(); cc++) {
			rand_coords[cc] /= total_rand;
			rand_global_coords[cc] /= total_rand;	
		}
		if(be_chatty) {
			cout << "testing derviative at global:";
			for(unsigned int cc = 0; cc < element.get_dimension(); cc++) {
				cout << " " << rand_global_coords[cc];	
			}
			cout << "\n"; 		
		}
		
		// for every shape function	
		for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
			// for all derivatives
			for(unsigned int dd = 0; dd < element.get_dimension(); dd++) {
				double numerical  = get_derivative(ii, dd, rand_global_coords[0], rand_global_coords[1], rand_global_coords[2]);
				double analytical = element.evaluate_form_function_derivative(dd, ii, &rand_coords[0]);			
				TSM_ASSERT_DELTA(message, numerical, analytical, ERROR_THRESHOLD);				 	 
				if(be_chatty) {
					double quot  = 1.0;
					quot = analytical != 0.0 ? analytical : 1.0; 
					quot = numerical  != 0.0 ? numerical  : quot;
					quot = tdkp_math::abs(quot);				
					double rel = tdkp_math::abs(numerical - analytical) / quot;					
					cout << "derivative " << dd << " of shape func " << ii << " "
					     << "numerical: " << numerical << " analytical: " 
					     << analytical << " "
					     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
						 << "\n";  						
				}
			}
		}			
	}	
}

void ElementNumericalValidator::test_integrate_solution() {
	
	double numerical = 0.0;
	vector<double> random_values;
	// integrate numerically
	for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
		double random_value = drand48();
		random_values.push_back(random_value);
		elem_shape_funcs[0] = ii;
		numerical  += random_value * evaluate_integral(&ElementNumericalValidator::zero_order_single);										
	}
	// integrate analytically
	double analytical = element.integrate_solution(&random_values[0]);
	// build output
	double quot  = 1.0;
	quot = analytical != 0.0 ? analytical : 1.0; 
	quot = numerical  != 0.0 ? numerical  : quot;
	quot = tdkp_math::abs(quot);				
	double rel = tdkp_math::abs(numerical - analytical) / quot;
	
	TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);
		
	if(be_chatty) {
		cout << "integrate solution (double): (random solution) "
		     << "numerical: " << numerical << " analytical: " 
		     << analytical << " "
		     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
			 << "\n";  			
	}
		
}

void ElementNumericalValidator::test_differentiate_solution() {
					
	// -----------------------------------------------------
	// differentiate numerically at middle of element (2d rect are implemented to give the derivative at (0,0) reference coordinate system ...)
	// -----------------------------------------------------
	// find middle of element
	double middle_coords[] = {0.0, 0.0, 0.0};	 
	for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
		for(unsigned int jj = 0; jj < element.get_dimension(); jj++) {
			middle_coords[jj] += element.get_node(ii).get_coord(jj);	
		}
	}
	for(unsigned int jj = 0; jj < element.get_dimension(); jj++) {
		middle_coords[jj] /= static_cast<double>(element.get_num_nodes());	
	}
	// for all derivatives				
	for(unsigned int dd = 0; dd < element.get_dimension(); dd++) {
		vector<double> random_values;
		double numerical = 0.0;
		// calculate derivatives
		for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {									
			double random_value = drand48();
			random_values.push_back(random_value);					
			numerical  += random_value * get_derivative(ii, dd, middle_coords[0], middle_coords[1], middle_coords[2]);
		}												
		// calculate analytically
		double analytical = element.differentiate_solution(dd, &random_values[0]);
		// build output
		double quot  = 1.0;
		quot = analytical != 0.0 ? analytical : 1.0; 
		quot = numerical  != 0.0 ? numerical  : quot;
		quot = tdkp_math::abs(quot);				
		double rel = tdkp_math::abs(numerical - analytical) / quot;
		
		TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);
		
		if(be_chatty) {
			cout << "differentiate " << dd << " solution (double): (random solution) "
			     << "numerical: " << numerical << " analytical: " 
			     << analytical << " "
			     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
				 << "\n";  			
		}
	}			
}	

void ElementNumericalValidator::test_zero_order_single() {
	for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
		elem_shape_funcs[0] = ii;
		double numerical  = evaluate_integral(&ElementNumericalValidator::zero_order_single);
		double analytical = element.get_single_integral_0th_order(ii);
		double quot  = 1.0;
		quot = analytical != 0.0 ? analytical : 1.0; 
		quot = numerical  != 0.0 ? numerical  : quot;
		quot = tdkp_math::abs(quot);				
		double rel = tdkp_math::abs(numerical - analytical) / quot;

		TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
		if(be_chatty) {

			cout << "zero order single (" << ii << ") "
			     << "numerical: "  << setw(12) << numerical << " "
			     << "analytical: " << setw(12) << analytical << " "
		  	     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
			     << "\n";
		}
	}	
}


void ElementNumericalValidator::test_first_order_single() {
	for(short ee = 0; ee < (signed)element.get_dimension(); ee++) {
		diffops[0] = ee;
		for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
			elem_shape_funcs[0] = ii;
				
			double numerical  = evaluate_integral(&ElementNumericalValidator::first_order_single);
			double analytical = element.get_single_integral_1st_order(ee,ii);
			double quot  = 1.0;
			quot = analytical != 0.0 ? analytical : 1.0; 
			quot = numerical  != 0.0 ? numerical  : quot;
			quot = tdkp_math::abs(quot);					
			double rel = tdkp_math::abs(numerical - analytical) / quot;

			TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
			if(be_chatty) {					
				
				cout << "first order partial " << ee << " single (" << ii << ") "
				     << "numerical: "  << setw(12) << numerical << " "
				     << "analytical: " << setw(12) << analytical << " "
				     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
				     << "\n"; 							    				     			     				     				
			}
		}	
	}	
}	


void ElementNumericalValidator::test_zero_order() {
	for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
		for(short jj = 0; jj < (signed)element.get_num_nodes(); jj++) {
			elem_shape_funcs[0] = ii;
			elem_shape_funcs[1] = jj;
			double numerical  = evaluate_integral(&ElementNumericalValidator::zero_order);
			double analytical = element.get_element_integral_0th_order(ii,jj);
			double quot  = 1.0;
			quot = analytical != 0.0 ? analytical : 1.0; 
			quot = numerical  != 0.0 ? numerical  : quot;
			quot = tdkp_math::abs(quot);				
			double rel = tdkp_math::abs(numerical - analytical) / quot;
			TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
			if(be_chatty) {
				cout << "zero order (" << ii << ", " << jj << ") "
				     << "numerical: "  << setw(12) << numerical << " "
				     << "analytical: " << setw(12) << analytical << " "
			  	     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
				     << "\n";
			}
		}	
	}	
}

void ElementNumericalValidator::test_zero_order_nodal() {
	for(short cc = 0; cc < (signed)element.get_num_nodes(); cc++) {
		nodal_data_func = cc;
		for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
			for(short jj = 0; jj < (signed)element.get_num_nodes(); jj++) {
				elem_shape_funcs[0] = ii;
				elem_shape_funcs[1] = jj;
				double numerical  = evaluate_integral(&ElementNumericalValidator::zero_order_nodal);
				double analytical = element.get_element_integral_0th_order_nodal_data(cc, ii,jj);
				double quot  = 1.0;
				quot = analytical != 0.0 ? analytical : 1.0; 
				quot = numerical  != 0.0 ? numerical  : quot;
				quot = tdkp_math::abs(quot);				
				double rel = tdkp_math::abs(numerical - analytical) / quot;
				TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
				if(be_chatty) {
					cout << "zero order nodal (nd:" << cc << "; " << ii << ", " << jj << ") "
					     << "numerical: "  << setw(12) << numerical << " "
					     << "analytical: " << setw(12) << analytical << " "
				  	     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
					     << "\n";
				}
			}	
		}
	}	
}

void ElementNumericalValidator::test_first_order() {
	
	for(short ee = 0; ee < (signed)element.get_dimension(); ee++) {
		diffops[0] = ee;
		for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
			elem_shape_funcs[0] = ii;
			for(short jj = 0; jj < (signed)element.get_num_nodes(); jj++) {
				
				elem_shape_funcs[1] = jj;
				
				double numerical  = evaluate_integral(&ElementNumericalValidator::first_order);
				double analytical = element.get_element_integral_1st_order(ee,ii,jj);
				double quot  = 1.0;
				quot = analytical != 0.0 ? analytical : 1.0; 
				quot = numerical  != 0.0 ? numerical  : quot;
				quot = tdkp_math::abs(quot);					
				double rel = tdkp_math::abs(numerical - analytical) / quot;
				TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
				if(be_chatty) {					
					cout << "first order partial " << ee << " (" << ii << ", " << jj << ") "
					     << "numerical: "  << setw(12) << numerical << " "
					     << "analytical: " << setw(12) << analytical << " "
					     << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")   
					     << "\n"; 							    
				}					     			     				     
			}					
		}	
	}	
}
void ElementNumericalValidator::test_second_order() {
	for(short ee = 0; ee < (signed)element.get_dimension(); ee++) {
		diffops[0] = ee;
		for(short ff = 0; ff < (signed)element.get_dimension(); ff++) {
			diffops[1] = ff;
			for(short ii = 0; ii < (signed)element.get_num_nodes(); ii++) {
				elem_shape_funcs[0] = ii;
				for(short jj = 0; jj < (signed)element.get_num_nodes(); jj++) {						
					elem_shape_funcs[1] = jj;
					double numerical  = evaluate_integral(&ElementNumericalValidator::second_order);
					double analytical = element.get_element_integral_2nd_order(ee,ii,ff,jj);
					double quot  = 1.0;
					quot = analytical != 0.0 ? analytical : 1.0; 
					quot = numerical  != 0.0 ? numerical  : quot;
					quot = tdkp_math::abs(quot);
					double rel = tdkp_math::abs(numerical - analytical) / quot;

					TSM_ASSERT_DELTA(message, analytical, numerical, ERROR_THRESHOLD);			
					if(be_chatty) {						
						cout << "second order partial " << ee << "/" << ff << " (" << ii << ", " << jj << ") "
						     << "numerical: "  << setw(12) << numerical << " "
						     << "analytical: " << setw(12) << analytical << " "						     
					     	 << "rel: " << rel << " " << (rel > ERROR_THRESHOLD ? "ERR" : "")						      
						     << "\n";
					}
				}
			}					
		}	
	}	
}	


void ElementNumericalValidator::plot_shape_functions(const char* basename) {
	
	const int rough_patch_num = 250; // was 400
	
	
	// plot element	
	int xx[] = {0, 0, 0};
	xx[1]    = (element.get_dimension() < 2 ? rough_patch_num : 0); 
	xx[2]    = (element.get_dimension() < 3 ? rough_patch_num : 0);
	
	double tmp_coordinates[] = {0.0, 0.0, 0.0};
	double tmp_increments[]  = {0.0, 0.0, 0.0};
	
	for(unsigned int ii = 0; ii < element.get_dimension(); ii++) {
		tmp_increments[ii] = (high.get_coord(ii) - low.get_coord(ii)) / double(rough_patch_num - 1);
	}		
	
	cout << "writing to " << basename << " ";
	
	ofstream fout(basename);
	
	while(xx[0] < rough_patch_num) {
		// calculate position
		for(short ii = 0; ii < (signed)element.get_dimension(); ii++) {
			tmp_coordinates[ii] = low.get_coord(ii) + double(xx[ii]) * tmp_increments[ii];
			fout << tmp_coordinates[ii] << " \t";
		}
		// test if its inside		
		if(element.inside_element(tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2])) {
			for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
				fout << element.evaluate_form_function_global(ii, tmp_coordinates[0], tmp_coordinates[1], tmp_coordinates[2]) 
				     << " \t";
			}						
		} else {
			for(unsigned int ii = 0; ii < element.get_num_nodes(); ii++) {
				fout << "0.0 \t";	
			}					
		} 
		fout << "\n";

		
		// -----------------------------------
		// update position (depending on dimension)
		// -----------------------------------
		xx[2]++;
		if(xx[2] >= rough_patch_num) {
			xx[2] = (element.get_dimension() < 3 ? rough_patch_num : 0);
			xx[1]++;
			if(xx[1] >= rough_patch_num) {
				xx[1] = (element.get_dimension() < 2 ? rough_patch_num : 0);
				xx[0]++;
			} 	
		}
		 	
	}	
	fout.close();
	cout << "done\n";
}

#endif /*VALIDATOR_H_*/
