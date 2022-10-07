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
int print_mercury6_big(std::string mercury6_big, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, double DATE);
int print_mercury6_param(std::string mercury6_param, OPTIONS *opt, double start_time);
int print_mercury6_small(std::string mercury6_small, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, std::vector<double> *temp_t, double rho, unsigned long long int MAX, unsigned long long int *particle_number_mercury);
int print_mercury6_small_pb(std::string mercury6_small, std::vector<std::vector<double> > *temp_q, std::vector<std::vector<double> > *temp_v, std::vector<double> *temp_m, std::vector<double> *temp_t, double rho);
int print_mercury6_element(std::string mercury6_element, std::vector<std::string> bodies,unsigned int format);

#endif /* DATA_WRITE_HH */


#ifndef SYSTEM_FUNCTIONS_HH_
#define SYSTEM_FUNCTIONS_HH_

std::string exec(const char* cmd);
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

double MOID_bonanno(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<std::vector<double> >* m_vec, unsigned int id_1, unsigned int id_2);
double Tisserand(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<std::vector<double> >* m_vec, unsigned int id_1, unsigned int id_2);

double nu_to_E(double e,double nu);

std::vector<double> total_angular_momentum_qv_heliocentric(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec);
int rot_to_invariable_plane(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, std::vector<double>* m_vec);
int rot_all_y(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta);
int rot_all_z(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec, double theta);

void heliocentric_qp(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* p_vec, std::vector<double>* m_vec);
void heliocentric_qv(std::vector<std::vector<double> >* q_vec, std::vector<std::vector<double> >* v_vec);

std::vector<double> qp_to_kepler(std::vector<double> q, std::vector<double> p, double m, double M);
std::vector<double> kepler_to_qp(std::vector<double> kepler, double m, double M);
std::vector<double> xv_to_kepler(std::vector<double> x, std::vector<double> v, double mu);
std::vector<double> kepler_to_xv(std::vector<double> kepler, double mu);

#endif /* KEPLER_HH_ */

#ifndef MATH_FUNCTIONS_HH_
#define MATH_FUNCTIONS_HH_

 double stdv_mat(std::vector<std::vector<double> > M, double M_mean);
 double mean_mat(std::vector<std::vector<double> > M);

double draw_from_smoth_dist(std::vector<double> dist, std::vector<double> bins);
unsigned int find_closest_XpowN(double X, unsigned int n);
unsigned int find_closest_pow2(unsigned int n);
std::vector<std::vector<double> > fill_mat(unsigned int a, unsigned int b, double X);
std::vector<std::vector<double> > generate_N_bin_list(std::vector<std::vector<double> > sample, unsigned int N_usr);
std::vector<double> generate_N_hist_list(std::vector<std::vector<double> > sample, std::vector<std::vector<double> > bins, std::vector<double> dx);
std::vector<double> draw_from_N_dist(std::vector<double> dist, std::vector<std::vector<double> > bins, std::vector<double> dx);

double abs_double(double a);

double JDtoJ(double JD);
double DATEtoJD(double year, double month, double day, double hour, double minute, double second);

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

std::vector<double> normal_distribution_density_vector(double mu, double sigma, double sigma_range, unsigned int n);
std::vector<double> normal_distribution_bin_mid_vector(double mu, double sigma, double sigma_range, unsigned int n);
std::vector<double> normal_distribution_density_vector_above_zero(double mu, double sigma, double sigma_range, unsigned int n);
std::vector<double> normal_distribution_bin_mid_vector_above_zero(double mu, double sigma, double sigma_range, unsigned int n);
double normal_distribution_density_function(double mu, double sigma, double x);

bool is_finite(double x);

#endif /* MATH_FUNCTIONS_HH_ */

#ifndef D_CRITERION_HH_
#define D_CRITERION_HH_

std::vector<std::vector<double> > calculate_metric_matrix(const std::vector<OBSERVATION_record> *list, std::string function);
std::vector<std::vector<double> > augment_metric_matrix(const std::vector<OBSERVATION_record> *list, const std::vector<OBSERVATION_record> *list_addition, std::string function);

double d_SH(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_D(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_V(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double d_N(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double rho2(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);
double varrho1(const OBSERVATION_record *orb1, const OBSERVATION_record *orb2);

#endif /* D_CRITERION_HH_ */

#ifndef MODELS_HH_
#define MODELS_HH_

double grun_cum_flux(double m);
std::vector<double> grun_flux(double min_mass_detectable, unsigned int grun_model_resolution);

#endif /* MODELS_HH_ */

#ifndef ASSOCIATION_HH_
#define ASSOCIATION_HH_


int node_associate_file(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster, std::string matrix_folder, std::vector<std::string> file_list, const std::vector<std::vector<unsigned int> > &index_matrix);
std::vector<std::vector<unsigned int> > single_linkage_clustering_files(std::string matrix_folder, std::vector<std::string> file_list, std::vector<std::vector<unsigned int> > index_matrix, double d_crit);

std::vector<std::vector<unsigned int> > single_linkage_clustering_matrix(std::vector<std::vector<double> > *d_matrix, double d_crit);
int node_associate(std::vector<node> *nodes, unsigned int *cluster_counter, double d_crit, unsigned int node_id, int mother_cluster);



#endif /* ASSOCIATION_HH_ */



#ifndef FUNCTIONS_HH_
#define FUNCTIONS_HH_


#endif /* FUNCTIONS_HH_ */


