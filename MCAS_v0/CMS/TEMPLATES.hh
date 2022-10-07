/*
 * TEMPLATES.hh
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */

#ifndef TEMPLATES_HH_
#define TEMPLATES_HH_

#include <math.h>
#include <string>
#include <vector>
 
#include "functions.hh"
#include "OVERLOADING.hh"

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


#endif /* TEMPLATES_HH_ */
