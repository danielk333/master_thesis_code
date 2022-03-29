/*
 * OVERLOADING.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */

#ifndef OVERLOADING_HH
#define OVERLOADING_HH

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
 
#include "class.hh"
#include "functions.hh"
/*
output_stream& operator<<(output_stream& O, std::ostream&(*f)(std::ostream&)) {
	switch(O.type) {
		case 0:
			O.out1_ << f;
			break;
		case 1:
			O.out1_ << f;
			O.out2_ << f;
			break;
		default:
			std::cout << "Could not read output type" << std::endl;
	}
   return O;
}
*/

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
std::vector<T> operator/(const std::vector<T> &v, const T &a) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((v[i])/(a));
	}
	return ret;
}
template <typename T>
std::vector<T> operator/(const T &a, const std::vector<T> &v) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((a)/(v[i]));
	}
	return ret;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &v, const T &a) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((v[i])*(a));
	}
	return ret;
}
template <typename T>
std::vector<T> operator*(const T &a, const std::vector<T> &v) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((a)*(v[i]));
	}
	return ret;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &v, const T &a) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((v[i])+(a));
	}
	return ret;
}
template <typename T>
std::vector<T> operator+(const T &a, const std::vector<T> &v) {
	std::vector<T> ret;
	unsigned int i;
	for(i = 0; i < v.size(); i++) {
		ret.push_back(a+v[i]);
	}
	return ret;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &v, const T &a) {
	std::vector<T> ret;
	unsigned int i;
	for(i=0; i < v.size(); i++) {
		ret.push_back((v[i])-(a));
	}
	return ret;
}
template <typename T>
std::vector<T> operator-(const T &a, const std::vector<T> &v) {
	std::vector<T> ret;
	unsigned int i;
	for(i = 0; i < v.size(); i++) {
		ret.push_back(a-v[i]);
	}
	return ret;
}


template <typename T> 
std::vector<std::vector<T> > operator*(const T &a, const std::vector<std::vector<T> > &M) {
	std::vector<T> temp_vec;
	std::vector<std::vector<T> > return_vec;
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

template <typename T> 
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &u) {
	std::vector<T> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(u[i]+a[i]);
	}
	return ret_vec;
}

template <typename T> 
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &u) {
	std::vector<T> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(a[i]-u[i]);
	}
	return ret_vec;
}

template <typename T> 
std::vector<T> operator*(const std::vector<T> &a, const std::vector<T> &u) {
	std::vector<T> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(a[i]*u[i]);
	}
	return ret_vec;
}
template <typename T> 
std::vector<T> operator/(const std::vector<T> &a, const std::vector<T> &u) {
	std::vector<T> ret_vec;
	unsigned int i;
	for(i=0;i<u.size();i++){
		ret_vec.push_back(a[i]/u[i]);
	}
	return ret_vec;
}

#endif
