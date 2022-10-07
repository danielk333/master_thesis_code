/*
 * math_functions.cc
 *
 *  Created on: Sep 15, 2014
 *      Author: dankas
 */

#include <vector>
#include <math.h>
#include <algorithm>    // std::random_shuffle
 
#include "define.hh"
#include "functions.hh"

std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M) {
	std::vector<std::vector<double> > Mt;
	std::vector<double> temp;
	std::vector<unsigned int> size;
	unsigned int i,j,max_size;
	for(i = 0; i < M.size(); i++) {
		size.push_back(M[i].size());
	}

	std::vector<unsigned int>::iterator TEMP_ID;
	TEMP_ID = std::max_element(size.begin(),size.end());
	max_size = *TEMP_ID;
	for(i = 0; i < max_size; i++) {
		temp.clear();
		for (j = 0; j < M.size(); j++) {
			if(size[j] > i) {
				temp.push_back(M[j][i]);
			}
		}
		Mt.push_back(temp);
	}

	return Mt;
}

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v) {
	std::vector<double> s;
	s.push_back(u[1]*v[2]-u[2]*v[1]);
	s.push_back(u[2]*v[0]-u[0]*v[2]);
	s.push_back(u[0]*v[1]-u[1]*v[0]);

	return s;
}

double dot_product(std::vector<double> u,std::vector<double> v) {
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

double abs_v(std::vector<double> u) {
	return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

double sum_v(std::vector<double> u) {
	double sum = 0;
	unsigned int i;
	for(i=0;i<u.size();i++){
		sum += u[i];
	}
	return sum;
}

std::vector<double> multiply_v(double a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]*a);
	}
	return ret_vec;
}

std::vector<double> add_v(double a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]+a);
	}
	return ret_vec;
}


std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M) {
	std::vector<double> temp_vec;
	std::vector<std::vector<double> > return_vec;
	unsigned int i,j;
	for(i=0; i < M.size(); i++) {
		temp_vec.clear();
		for(j=0; j < M[i].size(); j++) {
			temp_vec.push_back(a*((M[i])[j]));
		}
		return_vec.push_back(temp_vec);
	}

	return return_vec;
}

std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M) {
	std::vector<double> temp_vec;
	std::vector<std::vector<double> > return_vec;
	unsigned int i,j;
	for(i=0; i < M.size(); i++) {
		temp_vec.clear();
		for(j=0; j < M[i].size(); j++) {
			temp_vec.push_back(a[j] + ((M[i])[j]));
		}
		return_vec.push_back(temp_vec);
	}

	return return_vec;
}

std::vector<double> add_v_v(std::vector<double> a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]+a[i]);
	}
	return ret_vec;
}


unsigned int draw_from_dist(std::vector<double> dist) {
	unsigned int i;
	double rng = (double)rand() / RAND_MAX;
	double cum = 0;
	std::vector<double> temp = dist;
	double dist_SUM = sum_v(temp);

	if(!(dist_SUM > 1e-10 && dist_SUM < (1.0-1e-10))) {
		for(i=0; i < temp.size(); i++) {
			temp[i] = temp[i]/dist_SUM;
		}
	}

	for(i=0; i < temp.size(); i++) {
		if((rng < (cum+temp[i])) && (rng > cum)) {
			return i;
		}
		else {
			cum += temp[i];
		}
	}

	return 0;
}


std::vector<double> normal_distribution_density_vector(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(normal_distribution_density_function(mu,sigma,min_x + i*step));
	}

	return ret_vec;
}

std::vector<double> normal_distribution_bin_mid_vector(double mu, double sigma, double sigma_range, unsigned int n) {
	std::vector<double> ret_vec;
	unsigned int i;
	double min_x = mu - sigma_range*sigma;
	double max_x = mu + sigma_range*sigma;
	double step = (max_x - min_x)/n;

	for(i=0;i<n;i++){
		ret_vec.push_back(min_x + i*step);
	}

	return ret_vec;
}

double normal_distribution_density_function(double mu, double sigma, double x) {
	return (1/(sigma*sqrt(2*PI))*exp(-pow(x-mu,2)/(2*pow(sigma,2))));
}