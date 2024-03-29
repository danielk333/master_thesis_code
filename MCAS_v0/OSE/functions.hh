/*
 * functions.hh
 *
 *  Created on: Jun 17, 2015
 *      Author: dankas
 */


#include "define.hh"
#include "class.hh"

 #include "ast2body.h"

#ifndef DATA_READ_HH
#define DATA_READ_HH

int load_config(OPTIONS *opt,std::string file);
int load_data_vector(std::vector<double> *data,std::string file);
int load_data_matrix(std::vector<std::vector<double> > *data,std::string file);
int read_mercury6_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);
int read_family_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);
int read_orbdist_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);
int read_orbdist_ss_matrix(std::vector<std::vector<double> > *data,std::string matrix_file);
int fetch_rows_from_matrix_file(std::vector<std::vector<double> > *data, std::string matrix_file, unsigned int n, unsigned int m);

int load_file(const void *data, std::string file, int type);

#endif /* DATA_READ_HH */

#ifndef DATA_WRITE_HH
#define DATA_WRITE_HH

int save_mat(std::string out_file_name, std::vector<std::vector<double> > *M);
int save_batch(std::string file_name, std::string temp_file_name, std::vector<std::vector<double> > *M, unsigned int row_skip);
int print_mercury6_big(std::string mercury6_big, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, double date);
int print_mercury6_param(std::string mercury6_param, OPTIONS *opt);

#endif /* DATA_WRITE_HH */

#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

std::vector<std::string> list_files(const std::string& folder);
void remove_folder(const std::string& name);
bool file_exists(const std::string& name);
bool folder_exists(const std::string& name);
std::string get_selfpath();
bool create_file(std::string name);

#endif /* SYSTEM_FUNCTIONS_HH_ */

#ifndef KEPLER_HH_
#define KEPLER_HH_

std::vector<std::vector<double> > H_kep_qp(std::vector<double> q_vec, std::vector<double> p_vec, double m, double M, double dt);
std::vector<std::vector<double> > H_kep_xv(std::vector<double> q_vec, std::vector<double> v_vec, double M, double dt);
double universal_kepler_equation_newton(double dt, double tol, double a, double r_0, double v_r0, double mu);
double kepler_equation_newton(double M, double e, double E_0, double tol);

double nu_to_E(double e,double nu);

std::vector<double> total_angular_momentum_qv_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec);
int rot_to_invariable_plane(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec);

void heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);
void heliocentric_qv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec);

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu);
std::vector<double> kepler_to_xv(std::vector<double> kepler, double mu);

#endif /* KEPLER_HH_ */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

double abs_double(double a);

int factorial(int n);
int stumpff(int k, int n, double x);

 std::vector<double> rot_z(std::vector<double> u ,double theta);
 std::vector<double> rot_y(std::vector<double> u ,double theta);
 std::vector<double> rot_to_plane(std::vector<double> u ,std::vector<double> n);

std::vector<double> cross_product(std::vector<double> u,std::vector<double> v);
double dot_product(std::vector<double> u,std::vector<double> v);
double abs_v(std::vector<double> u);
std::vector<std::vector<double> > multiply_mat(double a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > add_v_mat(std::vector<double> a, std::vector<std::vector<double> > M);
std::vector<std::vector<double> > matrix_transpose(const std::vector<std::vector<double> > M);
std::vector<std::vector<double> > swap_cols(const std::vector<std::vector<double> > M, unsigned int a, unsigned int b);

std::vector<double> abs_element_v(std::vector<double> u);
std::vector<double> extract_col(std::vector<std::vector<double> > M, unsigned int id);
std::string IntToString(unsigned int a);

unsigned int draw_from_dist(std::vector<double> dist);
double draw_from_uniform(double a, double b);

std::vector<unsigned int> histogram(const std::vector<double> *data, unsigned int n);
std::vector<double> histogram_bins(const std::vector<double> *data, unsigned int n);

bool is_finite(double x);

#endif /* MATH_FUNCTIONS_HH_ */



#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_


#endif /* FUNCTIONS_HH_ */


