#include <vector>
#include <math.h>
#include <algorithm>
#include <sstream>
 
#include "resources.hh"
#include "functions.hh"


inline bool is_finite(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX); 
} 
inline double abs2(double x) {
	return x >= 0 ? x : -x;
}


tensor<double> cross_product(tensor<double> &u, tensor<double> &v) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3 || v.size(0) != 3) throw("Cannot perform cross product on size(0) != 3");
	#endif
	tensor<double> s(3);
	s(0) = u(1)*v(2)-u(2)*v(1);
	s(1) = u(2)*v(0)-u(0)*v(2);
	s(2) = u(0)*v(1)-u(1)*v(0);

	return s;
}
tensor<double> cross_product(const tensor<double> &u, const tensor<double> &v) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3 || v.size(0) != 3) throw("Cannot perform cross product on size(0) != 3");
	#endif
	tensor<double> s(3);
	s(0) = u(1)*v(2)-u(2)*v(1);
	s(1) = u(2)*v(0)-u(0)*v(2);
	s(2) = u(0)*v(1)-u(1)*v(0);

	return s;
}

double dot_product(tensor<double> &u, tensor<double> &v) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3 || v.size(0) != 3) throw("Cannot perform dot product on size(0) != 3");
	#endif
	return (u(0)*v(0) + u(1)*v(1) + u(2)*v(2));
}
double dot_product(const tensor<double> &u, const tensor<double> &v) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3 || v.size(0) != 3) throw("Cannot perform dot product on size(0) != 3");
	#endif
	return (u(0)*v(0) + u(1)*v(1) + u(2)*v(2));
}


double norm(const tensor<double> &u) {
	double R = 0;
	for(int i=0; i < u.size(); i++) R += inline_CUB(u(i));
	return sqrt(R);
}

double abs_v(tensor<double> &u) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3) throw("Cannot perform 3d norm on size(0) != 3");
	#endif
	return sqrt(u(0)*u(0) + u(1)*u(1) + u(2)*u(2));
}

double abs_v(const tensor<double> &u) {
	#ifdef _CHECK_BOUNDS_
		if(u.size(0) != 3) throw("Cannot perform 3d norm on size(0) != 3");
	#endif
	return sqrt(u(0)*u(0) + u(1)*u(1) + u(2)*u(2));
}

inline int factorial(int n) {
	int temp = 1,i;
	for(i=n;i > 0; i--) {
		temp = temp*i;
	}
	return temp;
}

inline void rot_z(tensor<double> &u ,double theta) {
	double uxTemp;
	uxTemp = u(0);
	u(0) = u(0)*cos(theta) - u(1)*sin(theta);
	u(1) = uxTemp*sin(theta) + u(1)*cos(theta);
	//u(2) = u(2);
}

inline void rot_y(tensor<double> &u ,double theta) {
	double uxTemp;
	uxTemp = u(0);
	u(0) = u(0)*cos(theta) + u(2)*sin(theta);
	//u(1) = u(1)
	u(2) = -uxTemp*sin(theta) + u(2)*cos(theta);
}

void rot_cols_z(tensor<double> &u ,double theta) {
	double uxTemp;
	for(int i=0; i < u.size(0); i++) {
		uxTemp=u(i,0);
		u(i,0) = u(i,0)*cos(theta) - u(i,1)*sin(theta);
		u(i,1) = uxTemp*sin(theta) + u(i,1)*cos(theta);
	}
}

void rot_cols_y(tensor<double> &u ,double theta) {
	double uxTemp;
	for(int i=0; i < u.size(0); i++) {
		uxTemp=u(i,0);
		u(i,0) = u(i,0)*cos(theta) + u(i,2)*sin(theta);
		u(i,2) = -uxTemp*sin(theta) + u(i,2)*cos(theta);
	}
}

void rot_rows_z(tensor<double> &u ,double theta) {
	double uxTemp;
	for(int i=0; i < u.size(1); i++) {
		uxTemp = u(0,i);
		u(0,i) = u(0,i)*cos(theta) - u(1,i)*sin(theta);
		u(1,i) = uxTemp*sin(theta) + u(1,i)*cos(theta);
	}
}

void rot_rows_y(tensor<double> &u ,double theta) {
	double uxTemp;
	for(int i=0; i < u.size(1); i++) {
		uxTemp = u(0,i);
		u(0,i) = u(0,i)*cos(theta) + u(2,i)*sin(theta);
		u(2,i) = -uxTemp*sin(theta) + u(2,i)*cos(theta);
	}
}


void rot_to_plane(tensor<double> &u ,const tensor<double> &n) {
	double phi = acos(n(2)/abs_v(n));
	double theta = atan2(n(1),n(0));
	rot_z(u,-theta);
	rot_y(u,-phi);
}



void rot_cols_to_plane(tensor<double> &u ,const tensor<double> &n) {
	double phi = acos(n(2)/abs_v(n));
	double theta = atan2(n(1),n(0));
	rot_cols_z(u ,-theta);
	rot_cols_y(u ,-phi);
}
void rot_rows_to_plane(tensor<double> &u ,const tensor<double> &n) {
	double phi = -acos(n(2)/abs_v(n));
	double theta = -atan2(n(1),n(0));
	rot_rows_z(u ,-theta);
	rot_rows_y(u ,-phi);
}
