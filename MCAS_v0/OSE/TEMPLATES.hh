/*
 * TEMPLATES.hh
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */

#ifndef TEMPLATES_HH_
#define TEMPLATES_HH_

#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <typeinfo>
 
#include "class.hh"
#include "functions.hh"
#include "OVERLOADING.hh"

template <typename T>
std::vector<T> sort_from_ref(std::vector<T> const& in,std::vector<std::pair<size_t, double_vector_iterator> > const& reference) {
    std::vector<T> ret(in.size());

    size_t const size = in.size();
    size_t i;
    for(i = 0; i < size; ++i) {
        ret[i] = in[reference[i].first];
    }
    return ret;
}

template <typename T>
T sum_v(std::vector<T> u) {
	T sum = 0;
	unsigned int i;
	for(i=0;i<u.size();i++){
		sum = sum + u[i];
	}
	return sum;
}

template <typename T>
 std::vector<double> normalize_v(std::vector<T> data) {
 	std::vector<double> ret;
 	unsigned int i;
 	T sum = sum_v(data);
 	for(i=0; i < data.size(); i++) {
 		ret.push_back((double)data[i]/((double)sum));
 	}
 	return ret;
 }

template <typename T>
std::vector<T> sum_mat(std::vector<std::vector<T> > M, std::string mode) {
	std::vector<T> ret;
	T temp_sum;
	unsigned int i,j;
	if(mode.compare("rows")) {
		for(i=0; i < M.size(); i++) {
			temp_sum = 0;
			for(j=0; j < M[i].size(); j++) {
				temp_sum += M[i][j];
			}
			ret.push_back(temp_sum);
		}
	}
	else if(mode.compare("cols")) {
		for(i=0; i < M[0].size(); i++) {
			temp_sum = 0;
			for(j=0; j < M.size(); j++) {
				temp_sum += M[j][i];
			}
			ret.push_back(temp_sum);
		}
	}

	return ret;
}

#ifndef OUT_OVERLOAD
#define OUT_OVERLOAD


template <typename T>
output_stream& operator<<(output_stream& O, T const& S) {
	switch(O.type) {
		case 0:
			O.out1_ << S;
			break;
		case 1:
			O.out1_ << S;
			O.out2_ << S;
			break;
		default:
			std::cout << "Could not read output type" << std::endl;
	}

	return O;
}

#endif

template <typename T>
void print_vector(std::vector<T> v, output_stream out) {
	unsigned int i;
	std::string temp_str1;
	temp_str1 = " ";
	std::cout << std::setprecision(10) << std::scientific;
	for(i=0; i < v.size(); i++) {
		out << v[i];
		if(i == (v.size()-1)) {
			out.endl();
		}
		else{
			out << temp_str1;
		}
	}
	
}

template <typename T>
void print_matrix(std::vector<std::vector<T> > M, output_stream out) {
	unsigned int i;
	std::string temp_str1,temp_str2;
	temp_str1 = "Row ";
	temp_str2 = ": ";
	for(i=0; i < M.size(); i++) {
		out << temp_str1 << i << temp_str2;
		print_vector(M[i],out);
	}
	
}

#endif /* TEMPLATES_HH_ */
