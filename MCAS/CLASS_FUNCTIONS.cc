/*
 * CLASS_FUNCTIONS.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */
#include <vector>
#include <math.h>
#include <sstream>
 
#include "define.hh"
#include "functions.hh"
#include "class.hh"

statistical_set::statistical_set(std::vector<double> DATA) {
	data = DATA;
}
statistical_set::statistical_set(const statistical_set &dist) {
    PDF = dist.PDF;
    hist = dist.hist;
    data = dist.data;
}
double statistical_set::mean(void) {
	return sum_v(data)/((double)data.size());
}
void statistical_set::calculate_hist(unsigned int n) {
	hist = histogram(&data, n);
	hist_bins = histogram_bins(&data, n);
}
void statistical_set::calculate_PDF(void) {
	PDF = normalize_v(hist);
}

void statistical_set::clear_data(void) {
	data.clear();
}

void statistical_set::convert_to_PDF_and_clear(unsigned int n) {
	calculate_hist(n);
	calculate_PDF();
	clear_data();
}
int statistical_set::save_histogram_to_file(std::string file) {
	std::vector<std::vector<double> > output;
	std::vector<double> temp;
	unsigned int i;
	for(i=0; i < hist.size(); i++) {
		temp.push_back((double)hist[i]);
	}
	output.push_back(hist_bins);
	output.push_back(temp);

	output = matrix_transpose(output);
	return save_mat(file, &output);
}
/*
void OBSERVATION_record::print(void) {

	std::cout << "-- Orbital elements: " << std::endl;
	std::cout << "a=" <<a << std::endl;
	std::cout << "e=" <<e << std::endl;
	std::cout << "i=" <<i << std::endl;
	std::cout << "omega=" <<omega << std::endl;
	std::cout << "Omega=" <<Omega << std::endl;

	std::cout << "-- Observational elements: " << std::endl;
	std::cout << "ra=" <<ra << std::endl;
	std::cout << "dec=" <<dec << std::endl;
	std::cout << "v_g=" <<v_g << std::endl;
	std::cout << "lambda=" <<lambda << std::endl;

	std::cout << "-- Record end" << std::endl;

}*/

void output_stream::endl() {
	switch(type) {
		case 0:
			out1_ << std::endl;
			break;
		case 1:
			out1_ << std::endl;
			out2_ << '\n';
			break;
		default:
			std::cout << "Could not read output type" << std::endl;
	}
}

node::node(std::vector<std::vector<double> > *d_mat,unsigned int ID_TAG) {
	unsigned int j;
    unsigned int objs = ((*d_mat).size());
    int ind_d; 

	if(ID_TAG > 0) {
		ind_d = ID_TAG-1;
        j=0;
        while(ind_d >= 0) {
           d_vec.push_back((*d_mat)[j][ind_d]);
           ind_d--;
           j++;
        }
     }
     d_vec.push_back(0);

    if(ID_TAG < objs) {
    	for(j=0; j < (*d_mat)[ID_TAG].size(); j++) {
     		d_vec.push_back((*d_mat)[ID_TAG][j]);
    	}
    }

    id = ID_TAG;
    member_of_cluster = FALSE;
}

void node::set_node(std::vector<std::vector<double> > *d_mat,unsigned int ID_TAG) {
	d_vec.clear();
	id = ID_TAG;
    member_of_cluster = FALSE;
    
	unsigned int j;
    unsigned int size = ((*d_mat).size());
    int ind_d; 

	if(ID_TAG > 0) {
		ind_d = ID_TAG-1;
        for(j=0; j < ID_TAG; j++) {
           d_vec.push_back((*d_mat)[j][ind_d]);
           ind_d--;
        }
    }
    d_vec.push_back(0);

    if(ID_TAG < size) {
    	for(j=0; j < (*d_mat)[ID_TAG].size(); j++) {
     		d_vec.push_back((*d_mat)[ID_TAG][j]);
    	}
    }

    id = ID_TAG;
    member_of_cluster = FALSE;
}
