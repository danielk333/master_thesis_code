/*
 * OVERLOADING.cc
 *
 *  Created on: Jul 2, 2015
 *      Author: dankas
 */

#ifndef OVERLOADING_HH
#define OVERLOADING_HH

#include <vector>



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
