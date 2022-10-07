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

#ifndef DATA_READ_HH
#define DATA_READ_HH

int load_data_vector(std::vector<double> *data,std::string file);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);

#endif /* DATA_READ_HH */

#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M);

#endif /* DATA_WRITE_HH */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M);
double sum_v(std::vector<double> u);
std::vector<double> multiply_v(double a, std::vector<double> u);
std::vector<double> add_v(double a, std::vector<double> u);
unsigned int draw_from_dist(std::vector<double> dist);
std::vector<double> add_v_v(std::vector<double> a, std::vector<double> M);
std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M);

std::vector<double> normal_distribution_density_vector(double mu, double sigma, double sigma_range, unsigned int n);
double normal_distribution_density_function(double mu, double sigma, double x);
std::vector<double> normal_distribution_bin_mid_vector(double mu, double sigma, double sigma_range, unsigned int n);

#endif /* MATH_FUNCTIONS_HH_ */

#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

bool file_exists(const std::string& name);
std::string get_selfpath();

#endif /* SYSTEM_FUNCTIONS_HH_ */

#ifndef KEPLER_HH_
#define KEPLER_HH_

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> qp_to_kepler_sun(std::vector<double> qs, std::vector<double> ps, std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> q, std::vector<double> p, double M);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);

#endif /* KEPLER_HH_ */

#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_

#endif /* FUNCTIONS_HH_ */
