/*
 * math_functions.cc
 *
 *  Created on: Sep 15, 2014
 *      Author: dankas
 */

#include <vector>
#include <math.h>
#include <algorithm>
#include <sstream>
 
#include "define.hh"
#include "functions.hh"


std::vector<size_t> sort_indexes_2(std::vector<double> x) {
	std::vector<size_t> ind;
	size_t id;
	std::vector<std::pair<size_t, double_vector_iterator> > order_str;
	unsigned int i;

    ind.resize(x.size());
    order_str.resize(x.size());
    id = 0;
    for (double_vector_iterator it = x.begin(); it != x.end(); ++it, ++id) {
        order_str[id] = std::make_pair(id, it);
    }

    std::sort(order_str.begin(), order_str.end(), double_vec_ordering_struct());

    for(i=0; i < x.size(); i++) {
		ind[i] = order_str[i].first;
	}

	return ind;
}

double neville(std::vector<double>* x,std::vector<double>* y,double x0) {
	unsigned int n = (*x).size();
	unsigned int i,j,k;
	std::vector<size_t> idx;
	std::vector<double> x_temp,y_temp;

	std::vector<double> xd;
	std::vector<std::vector<double> > P;

	x_temp.resize(n);
	y_temp.resize(n);
	xd.resize(n);

	for(i=0; i < n; i++) {
		xd[i] = abs_working((*x)[i] - x0);
	}

	idx = sort_indexes_2(xd);

	P = fill_mat(n,n,0.0);

	for(i=0; i < n; i++) {
		x_temp[i] = (*x)[idx[i]];
		y_temp[i] = (*y)[idx[i]];
		P[i][0] = (*y)[idx[i]];
	}

	for(i=1; i <= (n-1); i++) {
		for(j=1; j <= (n-i); j++) {
			P[j-1][i] = ( (x0 - x_temp[j-1])*P[j][i-1] + (x_temp[j+i-1] - x0)*P[j-1][i-1] )/(x_temp[j+i-1] - x_temp[j-1]);
		}
	}

	return P[0][n-1];
}

std::vector<std::vector<double> > fill_mat(unsigned int a, unsigned int b, double X) {
	std::vector<std::vector<double> > ret;
	std::vector<double> temp;
	unsigned int i,j;

	for(i=0; i < a; i++) {
		temp.clear();
		for(j=0; j < b; j++) {
			temp.push_back(X);
		}
		ret.push_back(temp);
	}

	return ret;
}

bool is_finite(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX); 
} 
double abs_working(double x) {
	if(x >= 0) {
		return x;
	}
	else {
		return -x;
	}
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

int factorial(int n) {
	int temp = 1,i;
	for(i=n;i > 0; i--) {
		temp = temp*i;
	}
	return temp;
}

int stumpff(int k, int n, double x) {
	double temp = 0;
	int i;

	for(i=0; i < n; i++) {
		temp = temp + pow(-x,i)/factorial(k+2*i);
	}

	return temp;
}


std::vector<double> rot_z(std::vector<double> u ,double theta) {
	std::vector<double> s;
	
	s.push_back(u[0]*cos(theta) - u[1]*sin(theta));
	s.push_back(u[0]*sin(theta) + u[1]*cos(theta));
	s.push_back(u[2]);

	return s;
}

std::vector<double> rot_y(std::vector<double> u ,double theta) {
	std::vector<double> s;
	
	s.push_back(u[0]*cos(theta) + u[2]*sin(theta));
	s.push_back(u[1]);
	s.push_back(-u[0]*sin(theta) + u[2]*cos(theta));

	return s;
}

std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n) {
	std::vector<double> s;

	double r = abs_v(n);
	double phi = acos(n[2]/r);
	double theta = atan2(n[1],n[0]);

	s = rot_y(rot_z(u,-theta),-phi);

	return s;
}
