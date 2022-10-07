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
#include "class.hh"



std::vector<std::vector<double> > uint_mat_to_double(std::vector<std::vector<unsigned int> > M) {
	unsigned int i,j;
	std::vector<std::vector<double> > RET;
	std::vector<double> ROW;

	for(i=0; i < M.size(); i++) {
		ROW.clear();
		for(j=0; j < M[i].size(); j++) {
			ROW.push_back((double)(M[i][j]));
		}
		RET.push_back(ROW);
	}

	return RET;
}

double abs_double(double a) {
	return abs(a);
}

std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n) {
	std::vector<double> s;

	double r = abs_v(n);
	double phi = acos(n[2]/r);
	double theta = atan2(n[1],n[0]);

	s = rot_y(rot_z(u,-theta),-phi);

	return s;
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

double min_v(std::vector<double> u) {
	double M = u[0];
	unsigned int i;

	for(i=1; i < u.size(); i++) {
		if(u[i] < M) {
			M = u[i];
		}
	}

	return M;
}
double max_v(std::vector<double> u) {
	double M = u[0];
	unsigned int i;

	for(i=1; i < u.size(); i++) {
		if(u[i] > M) {
			M = u[i];
		}
	}

	return M;
}
/*
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

std::vector<double> sub_v_v(std::vector<double> a, std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]-a[i]);
	}
	return ret_vec;
}
*/
std::vector<double> abs_element_v(std::vector<double> u) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(abs(u[i]));
	}
	return ret_vec;
}

std::vector<double> extract_col(std::vector<std::vector<double> > M, unsigned int id) {
	std::vector<double> ret_vec;
	unsigned int i;
	for(i=0;i<M.size();i++){
		ret_vec.push_back(M[i][id]);
	}
	return ret_vec;
}


std::string IntToString(unsigned int a) {
    std::ostringstream temp;
    temp<<a;
    return temp.str();
}


std::vector<unsigned int> histogram(const std::vector<double> *data, unsigned int n) {
	if(n == 0) {
		n = round(sqrt((*data).size()));
	}
	unsigned int i,j;
	std::vector<double> temp_data_container = *data;
	std::vector<double>::iterator MAX_ID, MIN_ID;
	MIN_ID = std::min_element(temp_data_container.begin(),temp_data_container.end());
	MAX_ID = std::max_element(temp_data_container.begin(),temp_data_container.end());
	double Delta = (*MAX_ID - *MIN_ID)/n;

	std::vector<unsigned int> index_vec;
	std::vector<unsigned int> ret;
	for(i=0; i < temp_data_container.size(); i++) {
		index_vec.push_back(i);
	}

	for(i=0; i < n; i++) {
		ret.push_back(0);
		for(j=index_vec.size(); j-- > 0; ) {
			if(temp_data_container[index_vec[j]] >= (*MIN_ID + ((double)i)*Delta) && temp_data_container[index_vec[j]] <= (*MIN_ID + ((double)(i+1))*Delta)) {
				ret[i]++;
				index_vec.erase(index_vec.begin()+j);
			}
		}
	}
	std::cout << "Histogram created with: " << std::endl;
	std::cout << "MIN = " << *MIN_ID << ", MAX = " << *MAX_ID  << ", Delta = " << Delta  << ", n = " << n << ", N = " << sum_v(ret) << std::endl;
	return ret;
}


std::vector<double> histogram_bins(const std::vector<double> *data, unsigned int n) {
	if(n == 0) {
		n = round(sqrt((*data).size()));
	}
	double min = *std::min_element((*data).begin(),(*data).end());
	double max = *std::max_element((*data).begin(),(*data).end());
	double Delta = (max - min)/n;

	std::vector<double> ret;

	unsigned int i;
	for(i=0; i < n; i++) {
		ret.push_back(min + Delta*(0.5 + (double)i));
	}
	return ret;
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

double draw_from_uniform(double a, double b) {
	return ((double)rand()/RAND_MAX)*(b - a) + a;
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

std::vector<std::vector<double> > swap_cols(const std::vector<std::vector<double> > M, unsigned int a, unsigned int b) {
	std::vector<std::vector<double> > Mt;
	unsigned int i;
	Mt = M;

	for(i = 0; i < M.size(); i++) {
		Mt[i][b] = M[i][a];
		Mt[i][a] = M[i][b];
	}

	return Mt;
}

bool is_finite(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX); 
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


