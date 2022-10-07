/*
 * functions.hh
 *
 *  Created on: May 6, 2015
 *      Author: dankas
 */

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
#include "define.hh"


#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(const char *out, std::vector<std::vector<double> > *M);

#endif /* DATA_WRITE_HH */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<std::vector<double> > generate_random_sphere_normals(unsigned int N);
std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M);
double sum_v(std::vector<double> u);
std::vector<double> multiply_v(double a, std::vector<double> u);
std::vector<double> add_v(double a, std::vector<double> u);
unsigned int draw_from_dist(std::vector<double> dist);
std::vector<double> add_v_v(std::vector<double> a, std::vector<double> M);

#endif /* MATH_FUNCTIONS_HH_ */

