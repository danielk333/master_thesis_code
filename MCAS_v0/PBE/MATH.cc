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

std::vector<size_t> sort_indexes(std::vector<double> x) {
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

	idx = sort_indexes(xd);

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
	double ret;
	unsigned int i;
	if(u.size() != v.size()) {
		ret = 0;
	}
	else {
		ret = 0;
		for(i=0; i < u.size(); i++) {
			ret = ret + u[i]*v[i];
		}
	}
	return ret;
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
	//N from 1 - infty, whole number
	//N is the number of normal vectors to be generated
	unsigned int n,i;
	double D;
	std::vector<double> temp_v;
	std::vector<double> temp_norm;
	std::vector<std::vector<double> > return_vec;
	n = 0;
	while(n < N) {
		temp_v.clear();
		D = 0;
		for(i=0; i < 4; i++) {
			temp_v.push_back( ((double)rand() / RAND_MAX)*2.0 - 1.0);
			D = D + temp_v.back()*temp_v.back();
		}
		
		if(D < 1) {
			n++;
			temp_norm.clear();
			temp_norm.push_back(2.0*(temp_v[1]*temp_v[3]+temp_v[0]*temp_v[2])/D);
			temp_norm.push_back(2.0*(temp_v[2]*temp_v[3]-temp_v[0]*temp_v[1])/D);
			temp_norm.push_back((temp_v[0]*temp_v[0] + temp_v[3]*temp_v[3] - temp_v[1]*temp_v[1] - temp_v[2]*temp_v[2])/D);
			//std::cout << "NORM " << abs_v(temp_norm) << std::endl;
			return_vec.push_back(temp_norm);
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
