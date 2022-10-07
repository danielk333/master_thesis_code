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


std::vector<std::vector<double> > generate_random_sphere_normals(unsigned int N) {
	//N from 0 - infty, whole number
	//n is the number of normal vectors to be generated
	//a is the longest side partition number of spherical corrdinate angle plane
	unsigned int n = pow(2,N*2+1), a = pow(2,N+1);
	unsigned int i,j;
	double dTheta = (2*PI)/a, dPhi = PI/a;
	std::vector<double> temp_v;
	std::vector<std::vector<double> > return_vec;
	temp_v.push_back(0);temp_v.push_back(0);temp_v.push_back(0);

	for(i=0; i < a; i++) {
		for(j=0; j < a/2; j++) {
			temp_v[0] = cos(i*dTheta)*sin(j*dPhi);
			temp_v[1] = sin(i*dTheta)*sin(j*dPhi);
			temp_v[2] = cos(j*dPhi);
			return_vec.push_back(temp_v);
		}
	}

	return return_vec;
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
	for(i=0; i < dist.size(); i++) {
		if((rng < (cum+dist[i])) && (rng > cum)) {
			return i;
		}
		else {
			cum += dist[i];
		}
	}

	return 0;
}