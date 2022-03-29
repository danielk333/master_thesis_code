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


#ifndef OVERLOADING_HH
#define OVERLOADING_HH

//Element scalar division
template <typename T>
tensor<T> operator/(const tensor<T> &v, const T &a) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]/a;
	}
	return X;
}
template <typename T>
tensor<T> operator/(const T &a, const tensor<T> &v) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = a/v[i];
	}
	return X;
}

//Element scalar multiplication
template <typename T>
tensor<T> operator*(const tensor<T> &v, const T &a) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]*a;
	}
	return X;
}
template <typename T>
tensor<T> operator*(const T &a, const tensor<T> &v) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = a*v[i];
	}
	return X;
}

//Element scalar addition
template <typename T>
tensor<T> operator+(const tensor<T> &v, const T &a) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]+a;
	}
	return X;
}
template <typename T>
tensor<T> operator+(const T &a, const tensor<T> &v) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = a+v[i];
	}
	return X;
}

//Element scalar subtraction
template <typename T>
tensor<T> operator-(const tensor<T> &v, const T &a) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]-a;
	}
	return X;
}
template <typename T>
tensor<T> operator-(const T &a, const tensor<T> &v) {
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = a-v[i];
	}
	return X;
}


//Element wise addition
//MUST BE SAME SIZE
template <typename T>
tensor<T> operator+(const tensor<T> &v, const tensor<T> &u) {
	#ifdef _CHECK_BOUNDS_
		if(v.size() != u.size()) {throw("Tensor TOTAL sizes does not agree in element wise addition.");}
	#endif
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]+u[i];
	}
	return X;
}


//Element wise subtraction
//MUST BE SAME SIZE
template <typename T>
tensor<T> operator-(const tensor<T> &v, const tensor<T> &u) {
	#ifdef _CHECK_BOUNDS_
		if(v.size() != u.size()) {throw("Tensor TOTAL sizes does not agree in element wise subtraction.");}
	#endif
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]-u[i];
	}
	return X;
}

//Element wise multiplication
//MUST BE SAME SIZE
template <typename T>
tensor<T> operator*(const tensor<T> &v, const tensor<T> &u) {
	#ifdef _CHECK_BOUNDS_
		if(v.size() != u.size()) {throw("Tensor TOTAL sizes does not agree in element wise multiplication.");}
	#endif
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]*u[i];
	}
	return X;
}

//Element wise division
//MUST BE SAME SIZE
template <typename T>
tensor<T> operator/(const tensor<T> &v, const tensor<T> &u) {
	#ifdef _CHECK_BOUNDS_
		if(v.size() != u.size()) {throw("Tensor TOTAL sizes does not agree in element wise division.");}
	#endif
	tensor<T> X;
	X.resize(v);
	for(long int i=0; i < v.size(); i++) {
		X[i] = v[i]/u[i];
	}
	return X;
}

#endif
